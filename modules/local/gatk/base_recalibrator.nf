process GATK_BASE_RECALIBRATOR {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(dedupedBam)
        path(dbsnp)
        path(ref_fasta)
        path("*")

    output:
        tuple val(meta), path("*.recal_data.table"), path(dedupedBam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gatk BaseRecalibrator --java-options "-Xmx${task.memory.giga}G" \\
        --known-sites ${dbsnp} \\
        -R ${ref_fasta} \\
        -I ${dedupedBam} \\
        -O ${meta.id}.recal_data.table
    """

    stub:
    """
    touch ${meta.id}.recal_data.table
    """
}
