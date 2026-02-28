# Tool Rationale — ATAC-seq LUAD Pipeline

## Aligner: Bowtie2 (not BWA, not HISAT2)

Bowtie2 is the standard for ATAC-seq because:
- Optimised for short-read gapped alignment (ATAC reads are typically 50-75bp)
- `--no-mixed --no-discordant` flags are well-implemented and critical for ATAC
- Lower memory footprint than BWA-MEM for short reads
- nf-core/atacseq, ENCODE ATAC pipeline, and TCGA all use Bowtie2

BWA is preferred for longer reads (>100bp) and variant calling.
HISAT2 is for splicing-aware RNA-seq alignment — not appropriate for ATAC.

---

## Peak caller: MACS2 (not MACS3, not HMMR, not Genrich)

MACS2 with `--nomodel --shift -100 --extsize 200` is the ENCODE-recommended
approach for ATAC-seq. The parameters centre signal on the Tn5 insertion site
rather than the fragment midpoint, which is critical for:
- Accurate motif enrichment (motif must align to accessibility peak centre)
- TF footprinting (requires precise Tn5 cut site coordinates)

MACS3 is the successor but not yet widely adopted in production pipelines.
Genrich handles multi-mapping better but is less validated for ATAC.
HMMR (nucleosome-free region caller) is complementary, not a replacement.

---

## Differential accessibility: DiffBind + DESeq2 (not edgeR, not limma-voom)

DiffBind provides ATAC-specific normalisation options that edgeR/limma lack:
- `DBA_LIBSIZE_PEAKREADS`: normalises by reads in peaks, not total library
  This is important because mito reads inflate raw library size in ATAC
- `DBA_NORM_RLE`: DESeq2-style RLE normalisation, robust for sparse count data
- Handles the consensus peak set construction automatically

DESeq2 backend is recommended over edgeR for ATAC because ATAC count
distributions are more negative-binomial (sparse, overdispersed) than
ChIP-seq, and DESeq2's shrinkage estimators handle this better for small
cohorts (n<10 per group).

---

## TF analysis: HOMER (enrichment) + TOBIAS (footprinting)

These tools answer different questions and are complementary:

| Question | Tool | Method |
|----------|------|--------|
| Which TF sequences are over-represented in open peaks? | HOMER | Sequence composition of peak set vs background |
| Which TFs are actually occupying their binding sites? | TOBIAS | Tn5 insertion dip (footprint) within open regions |
| Which TFs change binding between conditions? | TOBIAS BINDetect | Differential footprint scores |

Use both: HOMER for initial discovery (fast, sensitive), TOBIAS to
distinguish truly bound from accessible-but-unoccupied motifs (more specific).

A key interview point: a TF motif enriched by HOMER with NO footprint in
TOBIAS suggests the chromatin is permissive (open) but the TF is not present —
perhaps due to post-translational regulation or subcellular localisation.

---

## scATAC framework: ArchR (not Signac, not SnapATAC2)

ArchR is preferred for large datasets because:
- Disk-based Arrow file format — processes 100k+ cells without RAM exhaustion
- Iterative LSI handles batch effects and depth variation better than standard PCA
- Built-in gene activity score model (no matched RNA needed for annotation)
- chromVAR integration for TF deviation scores at single-cell resolution
- Active development team, best documentation in the field

Signac is the Seurat-integrated alternative — better if you need tight
scRNA + scATAC integration (Seurat WNN). SnapATAC2 is newer and faster
but less feature-complete for the analyses shown here.

For this portfolio project, ArchR is the stronger choice to demonstrate
because it is more commonly used in production ATAC-seq environments.

---

## Chromatin states: Roadmap Epigenomics (not ChromHMM from scratch)

Running chromHMM de novo requires multiple histone ChIP-seq experiments
(H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3, H3K9me3) that we don't
have for LUAD ATAC-seq data alone.

Instead we use the pre-computed Roadmap 15-state model for lung tissue
(E096) as a reference annotation. This is the standard approach for
interpreting ATAC peaks in the context of known chromatin states.

The key biological insight this enables: distinguishing whether
tumor-specific open chromatin sites are at:
- ACTIVE PROMOTERS (TssA): direct transcriptional activation
- ENHANCERS (Enh/EnhG): distal regulatory rewiring (more common in cancer)
- BIVALENT regions (TssBiv): polycomb derepression of developmental genes

This distinction is important because enhancer rewiring is considered a
hallmark of cancer epigenetic reprogramming — different mechanism from
simple promoter activation.
