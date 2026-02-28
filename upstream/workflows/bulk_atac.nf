// upstream/workflows/bulk_atac.nf
//
// Bulk ATAC-seq upstream pipeline
// ─────────────────────────────────────────────────────────────
// This is Nextflow in its natural habitat:
//   - Per-sample parallelism (each sample is independent)
//   - Sequential tool chaining (output of each step → next)
//   - Resource-heterogeneous steps (alignment needs 8 CPUs,
//     peak calling needs more memory)
//   - Cohort-level merge step (consensus peaks across all samples)
//
// ATAC-seq key biological concepts encoded here:
//   - Tn5 transposase inserts adapters preferentially at open chromatin
//   - NFR fragments (<200bp) mark active regulatory elements
//   - Mono/di-nucleosomal fragments (200-600bp) mark positioned nucleosomes
//   - Mitochondrial DNA is not relevant and must be removed
//   - ENCODE blacklist regions are technical artifacts
// ─────────────────────────────────────────────────────────────

nextflow.enable.dsl=2

include { FASTQC          } from '../modules/fastqc'
include { TRIM_GALORE     } from '../modules/trim_galore'
include { BOWTIE2_ALIGN   } from '../modules/bowtie2_align'
include { SAMTOOLS_FILTER } from '../modules/samtools_filter'
include { PICARD_DEDUP    } from '../modules/picard_dedup'
include { ATAQV           } from '../modules/ataqv'
include { MACS2_PEAKS     } from '../modules/macs2_peaks'
include { BEDTOOLS_MERGE  } from '../modules/bedtools_merge'
include { FEATURECOUNTS   } from '../modules/featurecounts'
include { TOBIAS_ATCORRECT} from '../modules/tobias_atcorrect'
include { TOBIAS_SCORE    } from '../modules/tobias_score'
include { MULTIQC         } from '../modules/multiqc'

workflow BULK_ATAC {
    take:
    ch_reads   // [meta, [R1, R2]]

    main:

    // ── Step 1: Raw read QC ───────────────────────────────────────────────
    FASTQC(ch_reads)

    // ── Step 2: Adapter trimming ──────────────────────────────────────────
    // Nextera adapters used in ATAC-seq library prep
    // Trim Galore auto-detects these
    // --paired: required for proper pair-end trimming
    TRIM_GALORE(ch_reads)

    // ── Step 3: Alignment to hg38 ─────────────────────────────────────────
    // Bowtie2 with ATAC-specific flags:
    //   --no-mixed:      only concordant paired alignments
    //   --no-discordant: discard pairs mapping to different chromosomes
    //   --maxins 2000:   allow large fragments (di/tri-nucleosomal)
    // These flags prevent spurious single-end alignments that inflate
    // open chromatin signal at repetitive elements
    BOWTIE2_ALIGN(
        TRIM_GALORE.out.trimmed_reads,
        params.bowtie2_index
    )

    // ── Step 4: BAM filtering ─────────────────────────────────────────────
    // Sequential filters applied in one samtools pipeline:
    //   a) Remove mitochondrial reads (chrM) — typically 20-80% of ATAC reads
    //      These are real but biologically uninformative for chromatin analysis
    //   b) Remove low MAPQ reads (MAPQ < 30) — multimappers
    //   c) Remove reads in ENCODE blacklist regions — known artifact loci
    //      (centromeres, telomeres, repetitive elements with anomalous signal)
    //   d) Keep only properly paired reads (-f 0x2)
    SAMTOOLS_FILTER(
        BOWTIE2_ALIGN.out.bam,
        params.blacklist
    )

    // ── Step 5: Remove PCR duplicates ─────────────────────────────────────
    // Picard MarkDuplicates: essential for ATAC-seq because library
    // complexity is typically lower than RNA-seq.
    // Track duplication rate in QC — >50% duplication = poor library quality
    PICARD_DEDUP(SAMTOOLS_FILTER.out.bam)

    // ── Step 6: ATAC-specific QC ──────────────────────────────────────────
    // ataqv computes metrics that standard RNA-seq QC tools miss:
    //   - TSS enrichment score: fold enrichment of reads at TSS vs background
    //     >6 = good, >10 = excellent ATAC library
    //   - Fragment size distribution: should show NFR peak (<200bp) and
    //     nucleosomal peaks at ~200bp, ~400bp, ~600bp
    //   - FRiP score: fraction of reads in peaks (computed after peak calling)
    //   - HQAA: high-quality autosomal alignments count
    ATAQV(
        PICARD_DEDUP.out.bam,
        params.fasta
    )

    // ── Step 7: Peak calling with MACS2 ───────────────────────────────────
    // ATAC-specific MACS2 parameters:
    //   --nomodel:    don't model fragment size (ATAC is multi-modal)
    //   --shift -100: shift reads upstream to center on Tn5 cut site
    //   --extsize 200: extend 200bp around cut site (= NFR window)
    //   --nolambda:   don't use local lambda for background (sparse ATAC)
    //   -q 0.05:      FDR threshold
    // Output: narrowPeak file per sample
    MACS2_PEAKS(PICARD_DEDUP.out.bam)

    // ── Step 8: Consensus peak set ────────────────────────────────────────
    // Merge all sample peaks into a single reference set:
    //   1. Concatenate all narrowPeak files
    //   2. Sort by coordinate
    //   3. Merge overlapping peaks
    //   4. Filter: keep only peaks present in ≥N samples
    // This consensus set is used as the feature space for differential analysis
    ch_all_peaks = MACS2_PEAKS.out.narrowpeak.collect()
    BEDTOOLS_MERGE(ch_all_peaks)

    // ── Step 9: Count matrix (samples × consensus peaks) ─────────────────
    // featureCounts quantifies how many reads from each sample fall into
    // each consensus peak → the input for DiffBind/DESeq2 downstream
    ch_all_bams = PICARD_DEDUP.out.bam.collect()
    FEATURECOUNTS(
        ch_all_bams,
        BEDTOOLS_MERGE.out.consensus_peaks
    )

    // ── Step 10: TOBIAS Tn5 bias correction ───────────────────────────────
    // TOBIAS (TF Occupancy and Binding from ATAC-seq Signals):
    //   Step A — ATACorrect: corrects for Tn5 insertion sequence bias
    //     Tn5 has a 9bp sequence preference that creates false peaks
    //     at certain motifs. Correction is done by comparing signal
    //     to expected from genome sequence context.
    //   Step B — ScoreBigwig: computes footprint scores per position
    //     High score = nucleosome-free + TF-bound (protected from Tn5)
    //     Low score  = open but unoccupied
    //   These outputs feed the footprinting analysis in downstream R
    if (params.run_tobias) {
        TOBIAS_ATCORRECT(
            PICARD_DEDUP.out.bam,
            params.fasta,
            BEDTOOLS_MERGE.out.consensus_peaks
        )
        TOBIAS_SCORE(TOBIAS_ATCORRECT.out.corrected_bw)
        ch_tobias = TOBIAS_SCORE.out.score_bw
    } else {
        ch_tobias = Channel.empty()
    }

    // ── MultiQC ───────────────────────────────────────────────────────────
    ch_qc = Channel.empty()
        .mix(FASTQC.out.zip)
        .mix(TRIM_GALORE.out.log)
        .mix(BOWTIE2_ALIGN.out.log)
        .mix(PICARD_DEDUP.out.metrics)
        .mix(ATAQV.out.json)
        .collect()
    MULTIQC(ch_qc)

    emit:
    bam              = PICARD_DEDUP.out.bam
    peaks            = MACS2_PEAKS.out.narrowpeak
    consensus_peaks  = BEDTOOLS_MERGE.out.consensus_peaks
    count_matrix     = FEATURECOUNTS.out.matrix     // → DiffBind downstream
    tobias_scores    = ch_tobias                    // → footprinting downstream
    ataqv_json       = ATAQV.out.json               // → QC report
    multiqc_html     = MULTIQC.out.report
}
