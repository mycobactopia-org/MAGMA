process GATK_APPLY_VQSR {
    tag "${meta.id}"
    label 'process_medium'

    input:
        val(analysisMode)
        tuple val(meta), path(variantsVcfIndex), path(variantsVcf), path(recalVcfIndex), path(recalVcf), path(tranches)
        path(reference)
        path("*")

    output:
        tuple val(meta), path("*.filtered_${analysisMode}_inc-rRNA.vcf.gz.tbi"), path("*.filtered_${analysisMode}_inc-rRNA.vcf.gz"), emit: filteredVcfTuple

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gatk ApplyVQSR --java-options "-Xmx${task.memory.giga}G" \\
        -R ${reference} \\
        -V ${variantsVcf} \\
        --tranches-file ${tranches} \\
        --recal-file ${recalVcf} \\
        ${args} \\
        -mode ${analysisMode} \\
        -O ${meta.id}.filtered_${analysisMode}_inc-rRNA.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.filtered_${analysisMode}_inc-rRNA.vcf.gz
    touch ${meta.id}.filtered_${analysisMode}_inc-rRNA.vcf.gz.tbi
    """
}
