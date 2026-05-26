process GATK_MARK_DUPLICATES {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(mergedBam)

    output:
        tuple val(meta), path("*.dedup_reads.bam"),    emit: bam_tuple
        tuple val(meta), path("*.MarkDupMetrics.txt"), emit: metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gatk MarkDuplicates --java-options "-Xmx${task.memory.giga}G" \\
        --METRICS_FILE ${meta.id}.MarkDupMetrics.txt \\
        -I ${mergedBam} \\
        -O ${meta.id}.dedup_reads.bam
    """

    stub:
    """
    touch ${meta.id}.MarkDupMetrics.txt
    touch ${meta.id}.dedup_reads.bam
    """
}
