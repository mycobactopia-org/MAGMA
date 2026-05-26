process LOFREQ_FILTER {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(vcf)
        path(ref_fasta)

    output:
        tuple val(meta), path("*.Filtered_AF.vcf")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    lofreq filter \\
        ${args} \\
        -i ${vcf} \\
    > ${meta.id}.Filtered_AF.vcf
    """

    stub:
    """
    touch ${meta.id}.Filtered_AF.vcf
    """
}
