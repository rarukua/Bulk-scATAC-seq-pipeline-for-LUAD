# downstream/R/07_scatac_archr.R
# Single-cell ATAC-seq analysis using ArchR
# ─────────────────────────────────────────────────────────────
# ArchR (Granja et al. 2021 Nat Genet) is the current standard
# for scATAC-seq analysis. Key advantages over Signac:
#   - Disk-based storage (handles 100k+ cells without RAM crash)
#   - Built-in iterative LSI (handles batch effects)
#   - Integrated gene score model (no RNA needed for annotation)
#   - Direct peak-to-gene linkage
#   - Trajectory analysis
#
# This script covers:
#   1. ArchR project creation from fragments.tsv.gz
#   2. QC filtering (TSS enrichment, fragment count)
#   3. Dimensionality reduction (LSI → UMAP)
#   4. Clustering (Seurat/scran graph-based)
#   5. Cell type annotation (gene scores + known markers)
#   6. Peak calling per cluster (pseudobulk)
#   7. Differential accessibility between clusters
#   8. TF deviation scores (chromVAR)

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(dplyr)
  library(optparse)
})

set.seed(42)
addArchRThreads(threads=8)
addArchRGenome("hg38")   # downloads hg38 annotations if not cached

opt_list <- list(
  make_option("--fragments_dir", type="character",
              help="Directory with fragments.tsv.gz files from Nextflow scATAC workflow"),
  make_option("--outdir",        type="character", default="results/scatac_archr"),
  make_option("--sample_name",   type="character", default="LUAD_scATAC")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
setwd(opt$outdir)

# ── Step 1: Create ArchR project ──────────────────────────────────────────
# Reads fragments.tsv.gz, computes per-barcode QC metrics,
# and creates an Arrow file (disk-based HDF5 storage)
message("── Creating ArchR project ──")

fragment_files <- list.files(opt$fragments_dir,
                              pattern="fragments.tsv.gz$",
                              full.names=TRUE, recursive=TRUE)

arrow_files <- createArrowFiles(
  inputFiles     = fragment_files,
  sampleNames    = tools::file_path_sans_ext(basename(dirname(fragment_files))),
  minTSS         = 4,       # minimum TSS enrichment score
  minFrags       = 1000,    # minimum unique fragments per barcode
  addTileMat     = TRUE,    # 500bp tile matrix for LSI
  addGeneScoreMat= TRUE     # gene activity scores (no RNA needed)
)

# ── Doublet detection ─────────────────────────────────────────────────────
# Doublets (two cells captured in one droplet) are a major scATAC artifact.
# ArchR's addDoubletScores uses a simulation approach:
# synthetic doublets are created and cells with similar profiles are flagged.
doublet_scores <- addDoubletScores(
  input       = arrow_files,
  k           = 10,
  knnMethod   = "UMAP",
  LSIMethod   = 1
)

# ── Create ArchR project ───────────────────────────────────────────────────
proj <- ArchRProject(
  ArrowFiles         = arrow_files,
  outputDirectory    = "ArchR_LUAD",
  copyArrows         = TRUE
)

# Filter doublets
proj <- filterDoublets(proj)
message(sprintf("Cells after doublet filtering: %d", nCells(proj)))

# ── Step 2: Dimensionality reduction (LSI) ────────────────────────────────
# LSI (Latent Semantic Indexing) is the standard for scATAC:
#   1. TF-IDF normalization of tile matrix (accounts for variable depth)
#   2. SVD to reduce dimensionality
#   3. L2 normalization of LSI components
# Iterative: first LSI removes depth effects, second is more biological
proj <- addIterativeLSI(
  ArchRProj   = proj,
  useMatrix   = "TileMatrix",
  name        = "IterativeLSI",
  iterations  = 2,
  clusterParams = list(resolution=c(0.2), sampleCells=10000, n.start=10),
  varFeatures = 25000,
  dims        = 1:30
)

# ── Step 3: Clustering ─────────────────────────────────────────────────────
# Graph-based clustering using Seurat's Louvain algorithm
# Resolution 0.4-0.8 recommended for tumor heterogeneity studies
proj <- addClusters(
  input      = proj,
  reducedDims= "IterativeLSI",
  method     = "Seurat",
  name       = "Clusters",
  resolution = 0.5
)
message("Clusters found: ", paste(unique(proj$Clusters), collapse=", "))

# ── Step 4: UMAP embedding ────────────────────────────────────────────────
proj <- addUMAP(
  ArchRProj   = proj,
  reducedDims = "IterativeLSI",
  name        = "UMAP",
  nNeighbors  = 30,
  minDist     = 0.5
)

# Plot UMAP colored by cluster
p_umap_cluster <- plotEmbedding(
  ArchRProj = proj,
  colorBy   = "cellColData",
  name      = "Clusters",
  embedding = "UMAP"
)
plotPDF(p_umap_cluster, name="UMAP_clusters", width=6, height=6,
        addDOC=FALSE, ArchRProj=proj)

# ── Step 5: Cell type annotation via gene scores ──────────────────────────
# ArchR gene activity scores: accessibility around gene body / TSS
# Used to infer expression without RNA-seq.
# LUAD cell type markers:
marker_genes_luad <- list(
  Tumor_cells    = c("NKX2-1","EPCAM","KRT8","KRT18","SFTA3"),
  T_cells        = c("CD3D","CD3E","CD4","CD8A","FOXP3"),
  B_cells        = c("CD79A","MS4A1","CD19"),
  Myeloid        = c("LYZ","CD68","FCGR3A","S100A8"),
  Fibroblasts    = c("COL1A1","FAP","ACTA2","PDGFRA"),
  Endothelial    = c("PECAM1","VWF","CDH5"),
  NK_cells       = c("NCAM1","KLRB1","NKG7")
)

# Plot gene scores for marker genes on UMAP
p_markers <- plotEmbedding(
  ArchRProj    = proj,
  colorBy      = "GeneScoreMatrix",
  name         = unlist(marker_genes_luad),
  embedding    = "UMAP",
  quantCut     = c(0.01, 0.95),
  imputeWeights= NULL
)
plotPDF(p_markers, name="UMAP_marker_gene_scores", width=5, height=5,
        addDOC=FALSE, ArchRProj=proj)

# Manually assign cell types based on marker scores
# (in a real analysis, use automated methods like scType or SingleR)
cell_type_map <- c(
  "C1"="Tumor", "C2"="T_cell", "C3"="Myeloid",
  "C4"="Fibroblast", "C5"="B_cell", "C6"="Endothelial",
  "C7"="Tumor", "C8"="NK_cell"
)
proj$CellType <- recode(proj$Clusters, !!!cell_type_map)

p_umap_celltype <- plotEmbedding(
  ArchRProj = proj, colorBy="cellColData",
  name="CellType", embedding="UMAP"
)
plotPDF(p_umap_celltype, name="UMAP_cell_types", width=6, height=6,
        addDOC=FALSE, ArchRProj=proj)

# ── Step 6: Peak calling per cluster (pseudobulk) ─────────────────────────
# Call peaks for each cell type separately using MACS2.
# Pseudobulk = aggregate all cells of a type into one "sample"
# before peak calling — essential for sparse scATAC data.
# This gives cell-type-specific peak sets.
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy   = "CellType"
)

proj <- addReproduciblePeakSet(
  ArchRProj  = proj,
  groupBy    = "CellType",
  pathToMacs2= findMacs2()
)

# Add peak matrix (cells × peaks) for downstream DA analysis
proj <- addPeakMatrix(proj)

# ── Step 7: chromVAR TF deviation scores ──────────────────────────────────
# chromVAR computes per-cell TF activity scores from motif accessibility.
# Unlike TOBIAS (which requires high coverage per sample), chromVAR
# works at single-cell resolution by aggregating signal across
# all instances of each motif genome-wide.
proj <- addMotifAnnotations(ArchRProj=proj, motifSet="JASPAR2022",
                             name="Motif")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj    = proj,
  peakAnnotation = "Motif",
  force          = TRUE
)

# Plot TF deviation scores
luad_motifs <- c("NKX2-1","FOXA1","AP-1","TP63","CEBPA","SOX2")
p_tf_dev <- plotEmbedding(
  ArchRProj    = proj,
  colorBy      = "MotifMatrix",
  name         = paste0("z:", luad_motifs),
  embedding    = "UMAP",
  imputeWeights= getImputeWeights(proj)
)
plotPDF(p_tf_dev, name="UMAP_TF_deviations", width=5, height=5,
        addDOC=FALSE, ArchRProj=proj)

# ── Save project ──────────────────────────────────────────────────────────
saveArchRProject(ArchRProj=proj, outputDirectory="ArchR_LUAD", load=FALSE)
message("── scATAC ArchR analysis complete ──")
message("  Project saved to: ArchR_LUAD/")
message("  Proceed to 08_scatac_pseudobulk.R for differential accessibility")
