process GATK_HAPLOTYPE_CALLER {
    tag "${meta.id}"
    label 'process_medium'

    input:
        tuple val(meta), path(bai), path(bam)
        path(ref_fasta)
        path("*")

    output:
        tuple val(meta), path("*.g.vcf.gz.tbi"), path("*.g.vcf.gz"), emit: gvcf_ch
        tuple val(meta), path("*haplotype_caller.bam"),               emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gatk HaplotypeCaller --java-options "-Xmx${task.memory.giga}G" \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -ERC GVCF \\
        ${args} \\
        -bamout ${meta.id}.haplotype_caller.bam \\
        -O ${meta.id}.g.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.g.vcf.gz
    touch ${meta.id}.g.vcf.gz.tbi
    touch ${meta.id}.haplotype_caller.bam
    """
}
