process GATK_GENOTYPE_GVCFS {
    tag "${meta.id}"
    label 'process_high'

    input:
        tuple val(meta), path(combinedVcfIndex), path(combinedVcf)
        path(ref_fasta)
        path("*")

    output:
        tuple val(meta), path("*.raw_variants.vcf.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gatk GenotypeGVCFs --java-options "-Xmx${task.memory.giga}G" \\
        -R ${ref_fasta} \\
        -V ${combinedVcf} \\
        ${args} \\
        -O ${meta.id}.raw_variants.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.raw_variants.vcf.gz
    """
}
