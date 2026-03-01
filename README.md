# ATAC-seq Analysis Pipeline — Lung Adenocarcinoma (LUAD)

> **Portfolio project** | Public datasets: TCGA-LUAD (bulk) · 10x Genomics LUAD (scATAC)  
> Covers: Bulk ATAC-seq · Single-cell ATAC-seq · Full downstream integration

---

## Biological Question

Lung adenocarcinoma (LUAD) is characterized by transcription factor (TF)-driven
chromatin remodeling that silences tumor suppressors and activates oncogenic programs.
Key LUAD TFs include NKX2-1 (lineage identity), FOXA1 (pioneer factor), and AP-1
family members (oncogenic activation).

This pipeline asks:

1. Which genomic regions are **differentially accessible** between LUAD tumor and normal lung?
2. Which **transcription factor motifs** are enriched at tumor-specific open chromatin sites?
3. Do TFs show active **footprints** at their binding sites — distinguishing bound from unbound motifs?
4. At single-cell resolution, which **cell populations** drive chromatin remodeling in the tumor microenvironment?
5. Does chromatin accessibility at gene promoters correlate with **gene expression** (TCGA RNA-seq)?
6. How do ATAC-seq accessibility patterns relate to **DNA methylation** (cross-modality with CRC pipeline)?

---

## Why ATAC-seq is Perfect for Nextflow

Unlike the cfDNA downstream analysis, ATAC-seq upstream processing is genuinely
compute-intensive and parallelizable — exactly what Nextflow is built for:

```
Per-sample parallelism:    each sample processed independently
Tool chaining:             FastQC → Trim → Bowtie2 → Samtools → MACS2 → annotation
Resource heterogeneity:    alignment (CPU-heavy) vs peak calling (memory-heavy)
HPC necessity:             bulk cohort (30+ samples) × scATAC (10k+ cells)
```

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  UPSTREAM — Nextflow (this is where it belongs)                  │
│                                                                   │
│  Bulk ATAC workflow          scATAC workflow                      │
│  FASTQ → BAM → peaks         FASTQ → fragments → ArchR/Signac    │
│  per sample, parallelized    cell barcode demux + QC              │
└───────────────────┬─────────────────────────┬───────────────────┘
                    │ consensus peaks          │ cell × peak matrix
                    │ count matrix             │ cell metadata
                    ▼                          ▼
┌─────────────────────────────────────────────────────────────────┐
│  DOWNSTREAM — R (iterative, exploratory)                         │
│                                                                   │
│  Bulk:                         scATAC:                            │
│  DiffBind differential access  ArchR: LSI + UMAP + clusters      │
│  HOMER motif enrichment        Cell type annotation               │
│  TOBIAS TF footprinting        Pseudobulk differential access     │
│  Chromatin state (chromHMM)    TF footprinting per cell type      │
│  RNA-seq integration           Trajectory analysis                │
└─────────────────────────────────────────────────────────────────┘
```

---

## Repository Structure

```
atac-luad-pipeline/
│
├── README.md
│
├── upstream/
│   ├── workflows/
│   │   ├── bulk_atac.nf           FASTQ → peaks → count matrix
│   │   └── scatac.nf              FASTQ → fragments.tsv.gz → ArchR input
│   └── modules/
│       ├── fastqc.nf
│       ├── trim_galore.nf
│       ├── bowtie2_align.nf       ATAC-specific flags (--no-mixed --no-discordant)
│       ├── samtools_filter.nf     Remove mitochondrial, low MAPQ, duplicates
│       ├── picard_dedup.nf
│       ├── macs2_peaks.nf         Peak calling (--nomodel --shift -100 --extsize 200)
│       ├── bedtools_merge.nf      Consensus peak set across samples
│       ├── featurecounts.nf       Count matrix: samples × peaks
│       ├── multiqc.nf
│       └── ataqv.nf               ATAC-specific QC (TSS enrichment, fragment sizes)
│
├── downstream/
│   ├── R/
│   │   ├── 01_bulk_qc.R               Fragment size, TSS enrichment, FRiP score
│   │   ├── 02_differential_access.R   DiffBind + DESeq2 tumor vs normal
│   │   ├── 03_motif_enrichment.R      HOMER findMotifsGenome
│   │   ├── 04_tobias_footprinting.R   Parse + visualize TOBIAS output
│   │   ├── 05_chromatin_states.R      chromHMM state annotation
│   │   ├── 06_rnaseq_integration.R    Accessibility vs expression correlation
│   │   ├── 07_scatac_archr.R          ArchR: clustering, cell type annotation
│   │   ├── 08_scatac_pseudobulk.R     Pseudobulk DA between cell populations
│   │   └── 09_multimodal_summary.R    Integrate bulk + scATAC + RNA-seq
│   └── notebooks/
│       ├── Bulk_ATAC_Report.Rmd
│       └── scATAC_Report.Rmd
│
├── docs/
│   ├── methods.md
│   ├── atac_qc_metrics.md          What each QC metric means
│   ├── bulk_vs_scatac.md           When to use each
│   └── tool_rationale.md
│
├── samplesheets/
│   ├── bulk_samples.csv
│   └── scatac_samples.csv
│
└── assets/
    ├── blacklist_hg38.bed          ENCODE blacklist regions
    ├── multiqc_config.yml
    └── luad_tfbs_sites.bed         LUAD-relevant TF binding sites for footprinting
```

---

## Datasets

| Assay | Dataset | Source | Samples |
|-------|---------|--------|---------|
| Bulk ATAC | TCGA-LUAD ATAC-seq | GDC / Corces et al. 2018 *Nat Genet* | ~110 tumor + paired normal |
| scATAC | 10x Genomics LUAD | 10x public datasets | ~8,000 cells, LUAD tumor |
| RNA-seq (integration) | TCGA-LUAD RNA-seq | GDC | Matched to ATAC samples |

---

## Quickstart

```bash
# Bulk ATAC-seq
nextflow run upstream/workflows/bulk_atac.nf \
  --input samplesheets/bulk_samples.csv \
  --genome hg38 \
  --outdir results/bulk \
  -profile docker

# scATAC-seq
nextflow run upstream/workflows/scatac.nf \
  --input samplesheets/scatac_samples.csv \
  --genome hg38 \
  --outdir results/scatac \
  -profile docker

# After Nextflow: downstream R analysis
Rscript -e "rmarkdown::render('downstream/notebooks/Bulk_ATAC_Report.Rmd')"
Rscript -e "rmarkdown::render('downstream/notebooks/scATAC_Report.Rmd')"
```

---

## Key ATAC-seq Concepts Demonstrated

| Concept | Where shown | Why it matters |
|---------|-------------|----------------|
| Nucleosome-free regions (NFR) | Fragment size QC | <200bp fragments = open chromatin; key ATAC signal |
| TSS enrichment score | ataqv QC module | Gold-standard ATAC quality metric; >6 = good library |
| FRiP score | featureCounts | Fraction of reads in peaks; >20% expected |
| Tn5 insertion bias | bowtie2 module | +4/-5bp offset correction for footprinting accuracy |
| Blacklist filtering | samtools_filter | Remove ENCODE artifact regions before any analysis |
| Consensus peak set | bedtools_merge | Required for comparing peaks across samples |
| TF footprinting vs motif enrichment | TOBIAS vs HOMER | Footprinting distinguishes *bound* from *unbound* motifs |
| LSI dimensionality reduction | ArchR (scATAC) | Standard for scATAC; TF-IDF + SVD on sparse peak matrix |

---

## References

1. Buenrostro JD et al. (2015) ATAC-seq. *Nat Methods* — Original ATAC-seq method
2. Corces MR et al. (2018) TCGA ATAC-seq pan-cancer. *Science* — Dataset used
3. Granja JM et al. (2021) ArchR. *Nat Genet* — scATAC analysis framework
4. Bentsen M et al. (2020) TOBIAS. *Nat Commun* — TF footprinting
5. Zhang Y et al. (2008) MACS. *Genome Biol* — Peak calling
6. Ross-Innes CS et al. (2012) DiffBind. *Nature* — Differential binding/accessibility
