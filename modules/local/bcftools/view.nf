process BCFTOOLS_VIEW {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*.filtered.bcf.gz"), path("*.filtered.bcf.gz.csi")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view ${args} ${vcf} -o ${prefix}.filtered.bcf
    bgzip ${prefix}.filtered.bcf
    bcftools index ${prefix}.filtered.bcf.gz
    """

    stub:
    """
    touch ${meta.id}.filtered.bcf.gz
    touch ${meta.id}.filtered.bcf.gz.csi
    """
}
