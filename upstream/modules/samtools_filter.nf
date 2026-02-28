// upstream/modules/samtools_filter.nf
//
// Multi-step BAM filtering for ATAC-seq
// Each filter targets a specific source of artifact or noise

process SAMTOOLS_FILTER {
    tag "${meta.id}"
    publishDir "${params.outdir}/filtered_bams/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  blacklist_bed

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.filtered.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}_filter_stats.txt"), emit: stats

    script:
    def prefix = meta.id
    """
    # ── Record pre-filter read counts ─────────────────────────────────────
    echo "Pre-filter reads:" > ${prefix}_filter_stats.txt
    samtools flagstat ${bam} >> ${prefix}_filter_stats.txt

    # ── Filter 1: Remove mitochondrial reads ──────────────────────────────
    # Tn5 transposase preferentially inserts into mitochondrial DNA because
    # it is not chromatinized (lacks nucleosomes). This creates an enormous
    # fraction of reads (20-80%) that are uninformative for chromatin state.
    # We exclude by chromosome name.
    samtools view -h ${bam} \\
        | grep -v "\\b${params.mito_name}\\b" \\
        | samtools view -b \\
        > ${prefix}.no_mito.bam

    # ── Filter 2: Remove low mapping quality reads ────────────────────────
    # MAPQ < 30 indicates multi-mapping or ambiguous alignment.
    # For ATAC-seq at open chromatin, multi-mappers are particularly
    # problematic because repeat elements often appear accessible
    # when in reality it's just multi-mapping artifact.
    # -f 0x2: keep only properly paired reads
    # -F 0x4: remove unmapped
    # -F 0x8: remove mate unmapped
    # -F 0x100: remove secondary alignments
    # -F 0x400: remove optical duplicates (before Picard step)
    samtools view -b \\
        -q ${params.min_mapq} \\
        -f 0x2 \\
        -F 0x4 -F 0x8 -F 0x100 \\
        ${prefix}.no_mito.bam \\
        > ${prefix}.mapq_filtered.bam

    # ── Filter 3: Remove ENCODE blacklist regions ──────────────────────────
    # Blacklist regions are genomic loci with anomalously high signal in
    # any ENCODE experiment regardless of experiment type. They arise from:
    //   - Satellite repeats (centromeres, telomeres)
    //   - Ribosomal DNA repeats
    //   - Segmental duplications
    //   - Virus integration hotspots
    # bedtools intersect -v: keep reads NOT overlapping blacklist
    # -abam: input is BAM format
    bedtools intersect \\
        -abam ${prefix}.mapq_filtered.bam \\
        -b ${blacklist_bed} \\
        -v \\
        > ${prefix}.filtered.bam

    samtools sort -@ ${task.cpus} -o ${prefix}.tmp.bam ${prefix}.filtered.bam
    mv ${prefix}.tmp.bam ${prefix}.filtered.bam
    samtools index ${prefix}.filtered.bam

    # ── Record post-filter stats ───────────────────────────────────────────
    echo "" >> ${prefix}_filter_stats.txt
    echo "Post-filter reads:" >> ${prefix}_filter_stats.txt
    samtools flagstat ${prefix}.filtered.bam >> ${prefix}_filter_stats.txt

    # ── Cleanup intermediates ──────────────────────────────────────────────
    rm -f ${prefix}.no_mito.bam ${prefix}.mapq_filtered.bam
    """
}
