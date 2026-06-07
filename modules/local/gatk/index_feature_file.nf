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
    // torch-magma's GATK_INDEX_FEATURE_FILE emits an explicit `-O <vcf>.tbi` when
    // an output prefix is set (minor-variants / LOFREQ path) and omits it for the
    // cohort path. We reproduce that per-alias via ext.args in conf/modules.config
    // so the emitted command matches standard MAGMA exactly (the .tbi produced is
    // identical either way — gatk defaults to <input>.tbi).
    def args = task.ext.args ?: ''
    """
    gatk IndexFeatureFile --java-options "-Xmx${task.memory.giga}G" \\
        ${args} \\
        -I ${vcf}
    """

    stub:
    """
    touch ${vcf}.tbi
    """
}
