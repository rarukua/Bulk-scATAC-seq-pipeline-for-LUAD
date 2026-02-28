# downstream/R/02_differential_access.R
# Differential chromatin accessibility: LUAD Tumor vs Normal
# Tools: DiffBind (count matrix handling) + DESeq2 (statistical testing)

suppressPackageStartupMessages({
  library(DiffBind)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(dplyr)
  library(optparse)
})

opt_list <- list(
  make_option("--count_matrix", type="character", help="featureCounts matrix from Nextflow"),
  make_option("--samplesheet",  type="character", help="Sample metadata CSV"),
  make_option("--peaks_dir",    type="character", help="Directory of per-sample narrowPeak files"),
  make_option("--outdir",       type="character", default="results/differential")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ── Build DiffBind sample sheet ────────────────────────────────────────────
# DiffBind requires: SampleID, Condition, bamReads, Peaks
message("── Loading samples into DiffBind ──")
meta <- read.csv(opt$samplesheet, stringsAsFactors=FALSE)

db_sheet <- data.frame(
  SampleID  = meta$sample,
  Condition = meta$condition,        # "Tumor" or "Normal"
  bamReads  = meta$bam_path,
  Peaks     = file.path(opt$peaks_dir, paste0(meta$sample, "_peaks.narrowPeak")),
  PeakCaller = "macs"
)

# ── Load count data ────────────────────────────────────────────────────────
dba <- dba(sampleSheet=db_sheet)

# Count reads in consensus peak set
# summits=250: re-center peaks on summit ± 250bp (standardizes peak width)
# This is important for DiffBind because variable-width peaks make
# comparisons difficult — re-centering normalizes the feature space
dba <- dba.count(dba, summits=250, bUseSummarizeOverlaps=TRUE)

# ── Normalize ─────────────────────────────────────────────────────────────
# RLE normalization (DESeq2 style) is recommended for ATAC-seq over
# library-size normalization because ATAC signal is inherently sparse.
# normalize=DBA_NORM_RLE, library=DBA_LIBSIZE_PEAKREADS:
# normalizes only by reads in peaks (not total library size)
# — important because mito reads can inflate total library estimates
dba <- dba.normalize(dba,
                     normalize = DBA_NORM_RLE,
                     library   = DBA_LIBSIZE_PEAKREADS)

# ── Contrast setup ─────────────────────────────────────────────────────────
dba <- dba.contrast(dba,
                    contrast   = c("Condition", "Tumor", "Normal"),
                    minMembers = 2)

# ── Differential analysis ──────────────────────────────────────────────────
# DESeq2 backend: Wald test with Benjamini-Hochberg FDR
dba <- dba.analyze(dba, method=DBA_DESEQ2)

# ── Extract results ────────────────────────────────────────────────────────
da_results <- dba.report(dba,
                          th         = 1,        # return all peaks
                          bCounts    = TRUE,
                          bNormCounts= TRUE)

da_df <- as.data.frame(da_results) %>%
  mutate(
    direction = case_when(
      FDR < 0.05 & Fold > 1  ~ "More open in Tumor",
      FDR < 0.05 & Fold < -1 ~ "More open in Normal",
      TRUE                   ~ "NS"
    )
  )

message(sprintf("  More open in Tumor:  %d peaks", sum(da_df$direction=="More open in Tumor")))
message(sprintf("  More open in Normal: %d peaks", sum(da_df$direction=="More open in Normal")))

# ── Annotate peaks ─────────────────────────────────────────────────────────
# ChIPseeker annotates each peak with nearest gene + genomic feature
# (promoter, intron, exon, intergenic, etc.)
da_gr <- makeGRangesFromDataFrame(da_df, keep.extra.columns=TRUE)
peak_anno <- annotatePeak(da_gr, tssRegion=c(-3000, 3000), TxDb=txdb,
                           annoDb="org.Hs.eg.db")
da_annotated <- as.data.frame(peak_anno)

write.csv(da_annotated,
          file.path(opt$outdir, "differential_peaks_annotated.csv"),
          row.names=FALSE)

# ── Volcano plot ───────────────────────────────────────────────────────────
# Highlight LUAD-relevant genes at differentially accessible peaks
luad_genes <- c("NKX2-1","FOXA1","FOXA2","TP63","SOX2","KRAS",
                "EGFR","ALK","STK11","KEAP1","TP53","RB1")

da_plot <- da_df %>%
  mutate(
    label = ifelse(direction != "NS" &
                   grepl(paste(luad_genes, collapse="|"), da_annotated$SYMBOL),
                   da_annotated$SYMBOL, NA)
  )

p_volcano <- ggplot(da_plot, aes(x=Fold, y=-log10(FDR), color=direction)) +
  geom_point(alpha=0.4, size=0.8) +
  geom_label_repel(aes(label=label), size=3, max.overlaps=15, na.rm=TRUE) +
  scale_color_manual(values=c(
    "More open in Tumor"  = "#D7191C",
    "More open in Normal" = "#2C7BB6",
    "NS"                  = "grey70"
  )) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="grey40") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  labs(title="Differential Chromatin Accessibility: LUAD Tumor vs Normal",
       subtitle=sprintf("FDR<0.05, |Fold|>1 | Tumor-open: %d | Normal-open: %d",
                        sum(da_df$direction=="More open in Tumor"),
                        sum(da_df$direction=="More open in Normal")),
       x="Log2 Fold Change (Tumor / Normal)",
       y="-log10(FDR)",
       color=NULL) +
  theme_bw(base_size=13)

ggsave(file.path(opt$outdir, "volcano_differential_access.pdf"),
       p_volcano, width=8, height=6)

# ── Peak annotation pie chart ──────────────────────────────────────────────
pdf(file.path(opt$outdir, "peak_annotation_pie.pdf"), width=8, height=6)
plotAnnoPie(peak_anno, main="Genomic Annotation of Differential Peaks")
dev.off()

# ── TSS distance plot ──────────────────────────────────────────────────────
pdf(file.path(opt$outdir, "tss_distance_plot.pdf"), width=8, height=5)
plotDistToTSS(peak_anno,
              title="Distance from Differential Peaks to Nearest TSS")
dev.off()

# ── Heatmap of top DA peaks ───────────────────────────────────────────────
library(ComplexHeatmap)
top_peaks <- da_annotated %>%
  filter(FDR < 0.01) %>%
  slice_max(abs(Fold), n=50)

count_mat <- dba.peakset(dba, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
top_mat <- as.matrix(count_mat[rownames(top_peaks),
                                 grep("^Conc", colnames(count_mat))])

pdf(file.path(opt$outdir, "top_peaks_heatmap.pdf"), width=10, height=10)
Heatmap(t(scale(t(top_mat))),
        name             = "Scaled accessibility",
        show_row_names   = FALSE,
        column_labels    = meta$sample,
        top_annotation   = HeatmapAnnotation(
          condition = meta$condition,
          col       = list(condition=c(Tumor="#D7191C", Normal="#2C7BB6"))
        ),
        column_title = "Top 50 Differential Peaks — LUAD Tumor vs Normal")
dev.off()

message("── Differential accessibility analysis complete ──")
