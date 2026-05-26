process GATK_INDEX_FEATURE_FILE {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*.tbi"), path(vcf), emit: sample_vcf_tuple
        tuple path("*.tbi"), path(vcf),             emit: vcf_tuple

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gatk IndexFeatureFile --java-options "-Xmx${task.memory.giga}G" \\
        -I ${vcf}
    """

    stub:
    """
    touch ${vcf}.tbi
    """
}
