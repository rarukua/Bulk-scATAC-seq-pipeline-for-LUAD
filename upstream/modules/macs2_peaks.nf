// upstream/modules/macs2_peaks.nf
//
// Peak calling for ATAC-seq using MACS2
// Parameters chosen specifically for ATAC-seq (NOT ChIP-seq defaults)

process MACS2_PEAKS {
    tag "${meta.id}"
    publishDir "${params.outdir}/peaks/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_peaks.narrowPeak"), emit: narrowpeak
    tuple val(meta), path("${meta.id}_summits.bed"),      emit: summits
    tuple val(meta), path("${meta.id}_peaks.xls"),        emit: xls
    tuple val(meta), path("${meta.id}_treat_pileup.bdg"), emit: bedgraph

    script:
    def prefix = meta.id

    // MACS2 for ATAC-seq — parameter rationale:
    //
    // --nomodel
    //   Normally MACS2 models fragment size from paired-end data.
    //   For ATAC-seq we skip this because:
    //   a) Fragment size is multi-modal (NFR + mono + di nucleosomal)
    //   b) We want to center on the Tn5 cut site, not the fragment midpoint
    //
    // --shift -100
    //   Shift reads -100bp relative to 5' end. Combined with --extsize 200,
    //   this recenters signal on the Tn5 insertion site.
    //   Why -100? Tn5 creates a 9bp staggered cut; the effective cut site
    //   is at position +4 (forward strand) or -5 (reverse strand).
    //   The -100 shift moves the read to span ±100bp around the cut site.
    //
    // --extsize 200
    //   After shifting, extend each read to 200bp total.
    //   This creates a 200bp window centered on the Tn5 cut site,
    //   which approximates the nucleosome-free region width.
    //
    // --nolambda
    //   Don't use local lambda (local background estimation).
    //   ATAC-seq has sparse background — local lambda often creates
    //   false negatives in genuinely open regions.
    //   Use genome-wide lambda instead (--slocal 0 --llocal 0 is implied)
    //
    // -q 0.05: FDR threshold (q-value, not p-value)
    //
    // -g hs: effective genome size for human hg38
    //
    // --bdg: output bedGraph for visualization
    //
    // -B: save fragment pileup bedGraph (useful for QC inspection)
    """
    macs2 callpeak \\
        -t ${bam} \\
        -f BAMPE \\
        -n ${prefix} \\
        --outdir . \\
        --nomodel \\
        --shift -100 \\
        --extsize 200 \\
        --nolambda \\
        -q ${params.macs2_qvalue} \\
        -g ${params.macs2_gsize} \\
        --bdg \\
        -B \\
        2>&1 | tee ${prefix}_macs2.log

    # Sort narrowPeak by score for easier downstream use
    sort -k8,8nr ${prefix}_peaks.narrowPeak > ${prefix}_peaks.sorted.narrowPeak
    mv ${prefix}_peaks.sorted.narrowPeak ${prefix}_peaks.narrowPeak
    """
}
