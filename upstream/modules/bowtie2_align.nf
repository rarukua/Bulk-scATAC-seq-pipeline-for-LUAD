// upstream/modules/bowtie2_align.nf
//
// Bowtie2 alignment for ATAC-seq
// Key flags are ATAC-specific — different from RNA-seq or ChIP-seq

process BOWTIE2_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/alignments/${meta.id}", mode: 'copy',
               saveAs: { filename -> filename.endsWith('.log') ? filename : null }

    input:
    tuple val(meta), path(reads)
    path  bowtie2_index

    output:
    tuple val(meta), path("${meta.id}.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}_bowtie2.log"), emit: log

    script:
    def prefix = meta.id
    def r1 = reads[0]
    def r2 = reads.size() > 1 ? reads[1] : null
    def paired = r2 ? "-1 ${r1} -2 ${r2}" : "-U ${r1}"

    // ATAC-seq specific Bowtie2 flags:
    //
    // --very-sensitive: maximum sensitivity (-N 1 -L 20 -i S,1,0.50)
    //   ATAC requires this because open chromatin regions may have
    //   non-standard sequence composition
    //
    // --no-mixed: CRITICAL for ATAC
    //   Discard pairs where only one mate aligns. Mixed pairs create
    //   false single-end fragments that inflate apparent open chromatin
    //
    // --no-discordant: CRITICAL for ATAC
    //   Discard pairs aligning to different chromosomes. Discordant
    //   pairs arise from transposase jumping between loci — a real
    //   biology but a major source of false peaks
    //
    // --maxins 2000: allow large insert sizes
    //   ATAC fragments range from ~100bp (NFR) to ~600bp (trinucleosomal)
    //   Default max of 500bp would discard nucleosomal fragments
    //
    // --no-unal: don't output unaligned reads (smaller BAM)
    """
    INDEX=\$(find -L ${bowtie2_index} -name "*.1.bt2" | sed 's/.1.bt2//')

    bowtie2 \\
        -x \$INDEX \\
        ${paired} \\
        --very-sensitive \\
        --no-mixed \\
        --no-discordant \\
        --maxins 2000 \\
        --no-unal \\
        -p ${task.cpus} \\
        2> ${prefix}_bowtie2.log \\
    | samtools sort -@ ${task.cpus} -o ${prefix}.bam

    samtools index ${prefix}.bam
    """
}
