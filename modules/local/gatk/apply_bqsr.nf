process GATK_APPLY_BQSR {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(recalibrationTable), path(dedupedBam)
        path(ref_fasta)
        path("*")

    output:
        tuple val(meta), path("*.baserecalibrated.bam")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gatk ApplyBQSR --java-options "-Xmx${task.memory.giga}G" \\
        -R ${ref_fasta} \\
        -I ${dedupedBam} \\
        --bqsr ${recalibrationTable} \\
        -O ${meta.id}.baserecalibrated.bam
    """

    stub:
    """
    touch ${meta.id}.baserecalibrated.bam
    """
}
