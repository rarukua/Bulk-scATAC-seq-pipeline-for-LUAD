// upstream/workflows/scatac.nf
//
// Single-cell ATAC-seq upstream pipeline
// ─────────────────────────────────────────────────────────────
// scATAC-seq differs from bulk in two fundamental ways:
//
//   1. Cell barcodes: each read carries a barcode identifying which
//      cell it came from. The Nextflow upstream job is to:
//        a) Demultiplex reads by barcode
//        b) Create a fragments.tsv.gz file (the universal scATAC format)
//        c) Generate a per-barcode QC summary
//      The actual single-cell analysis (clustering, cell types, DA)
//      happens DOWNSTREAM in ArchR/Signac R packages.
//
//   2. Sparsity: each cell has far fewer reads than a bulk sample.
//      Standard bulk peak calling per-cell is not appropriate.
//      Instead we create a fragments file and let ArchR/Signac
//      build peak matrices using pseudobulk approaches.
//
// Input: 10x Genomics scATAC FASTQ (3-file format: R1=barcode, R2=insert, R3=insert)
// Output: fragments.tsv.gz → ArchR/Signac input
// ─────────────────────────────────────────────────────────────

nextflow.enable.dsl=2

include { FASTQC           } from '../modules/fastqc'
include { CELLRANGER_ATAC  } from '../modules/cellranger_atac'
include { SINTO_FRAGMENTS  } from '../modules/sinto_fragments'
include { BARCODE_QC       } from '../modules/barcode_qc'
include { MULTIQC          } from '../modules/multiqc'

workflow SCATAC {
    take:
    ch_reads   // [meta, [R1, R2, R3]] — 10x 3-file format

    main:

    // ── Step 1: Raw QC ────────────────────────────────────────────────────
    FASTQC(ch_reads)

    // ── Step 2: Cell Ranger ATAC ──────────────────────────────────────────
    // Cell Ranger ATAC is the standard 10x Genomics preprocessing tool.
    // It handles:
    //   - Barcode correction and whitelisting
    //   - Alignment to reference genome
    //   - Duplicate removal (per cell barcode, not globally)
    //   - Peak calling (per-sample pseudobulk)
    //   - fragments.tsv.gz generation
    //   - Per-barcode QC metrics (TSS enrichment, fragment count)
    //
    // Alternative: if Cell Ranger license is restrictive, use
    // chromap (faster, open source) for alignment + sinto for fragment generation
    CELLRANGER_ATAC(
        ch_reads,
        params.cellranger_ref ?: "${params.genome}_cellranger_ref"
    )

    // ── Step 3: Fragment file generation (if using chromap instead) ───────
    // If using chromap as the aligner (faster alternative to Cell Ranger):
    // sinto fragments converts position-sorted BAM → fragments.tsv.gz
    // maintaining barcode information in the CB tag
    // Uncomment this block and comment out CELLRANGER_ATAC if using chromap
    // SINTO_FRAGMENTS(CHROMAP_ALIGN.out.bam)

    // ── Step 4: Barcode QC ────────────────────────────────────────────────
    // Filter barcodes by:
    //   - Minimum unique fragments per barcode (default: 500)
    //     Below this = empty droplet or dead cell
    //   - TSS enrichment per barcode (>6 = good quality cell)
    //   - Nucleosome signal ratio (NFR/mono-nucleosomal fragments)
    //
    // This step generates the "knee plot" — the inflection point in
    // log(fragments) distribution that separates real cells from noise
    BARCODE_QC(
        CELLRANGER_ATAC.out.fragments,
        CELLRANGER_ATAC.out.barcode_metrics
    )

    // ── MultiQC ───────────────────────────────────────────────────────────
    ch_qc = Channel.empty()
        .mix(FASTQC.out.zip)
        .mix(CELLRANGER_ATAC.out.summary)
        .collect()
    MULTIQC(ch_qc)

    emit:
    fragments        = CELLRANGER_ATAC.out.fragments       // → ArchR downstream
    barcode_qc       = BARCODE_QC.out.filtered_barcodes    // → ArchR downstream
    cellranger_peaks = CELLRANGER_ATAC.out.peaks           // → optional reference
    multiqc_html     = MULTIQC.out.report
}
