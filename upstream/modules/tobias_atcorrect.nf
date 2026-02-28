// upstream/modules/tobias_atcorrect.nf
//
// TOBIAS ATACorrect — Tn5 sequence bias correction
// ─────────────────────────────────────────────────────────────
// Reference: Bentsen et al. 2020 Nature Communications
//
// Why bias correction matters for footprinting:
//   Tn5 transposase has a strong sequence preference — it inserts
//   more efficiently at certain 9bp sequences. This means some
//   genomic regions appear "open" simply because they match the
//   Tn5 sequence motif, not because chromatin is actually accessible.
//
//   For motif ENRICHMENT (HOMER), this bias is less critical.
//   For TF FOOTPRINTING (TOBIAS), it is essential — a footprint
//   is defined as a DIP in Tn5 insertion within an otherwise open
//   region, so Tn5 bias that creates artificial dips would produce
//   false footprints.
//
//   ATACorrect:
//   1. Models Tn5 insertion frequency at each position from genome sequence
//   2. Subtracts expected from observed → corrected signal
//   3. Output: bias-corrected bigWig (one per strand)
// ─────────────────────────────────────────────────────────────

process TOBIAS_ATCORRECT {
    tag "${meta.id}"
    publishDir "${params.outdir}/tobias/atcorrect/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  peaks_bed

    output:
    tuple val(meta), path("${meta.id}_corrected.bw"), emit: corrected_bw
    tuple val(meta), path("${meta.id}_atcorrect.log")

    script:
    """
    TOBIAS ATACorrect \\
        --bam     ${bam} \\
        --genome  ${fasta} \\
        --peaks   ${peaks_bed} \\
        --outdir  . \\
        --prefix  ${meta.id} \\
        --cores   ${task.cpus} \\
        2>&1 | tee ${meta.id}_atcorrect.log

    # Rename output to consistent name
    mv ${meta.id}_corrected.bw ${meta.id}_corrected.bw 2>/dev/null || true
    """
}

// upstream/modules/tobias_score.nf
//
// TOBIAS ScoreBigwig — compute footprint score per position
// ─────────────────────────────────────────────────────────────
// The footprint score at each position measures whether a TF is
// likely BOUND at that position:
//
//   High score = protected from Tn5 (dip within open region) = TF bound
//   Low score  = accessible (high Tn5 insertion) = no TF bound
//
// The score is: mean(flanking signal) - mean(central signal)
// normalized by the standard deviation of the flanking signal.
//
// This score is then used in TOBIAS BINDetect downstream (in R)
// to classify each motif occurrence as BOUND, UNBOUND, or NO CHANGE
// between conditions (tumor vs normal).
// ─────────────────────────────────────────────────────────────

process TOBIAS_SCORE {
    tag "${meta.id}"
    publishDir "${params.outdir}/tobias/scores/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(corrected_bw)

    output:
    tuple val(meta), path("${meta.id}_footprints.bw"), emit: score_bw

    script:
    """
    TOBIAS ScoreBigwig \\
        --signal  ${corrected_bw} \\
        --regions ${params.consensus_peaks ?: 'peaks.bed'} \\
        --output  ${meta.id}_footprints.bw \\
        --cores   ${task.cpus}
    """
}
