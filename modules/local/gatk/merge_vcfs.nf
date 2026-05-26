process GATK_MERGE_VCFS {
    tag "${meta.id}"
    label 'process_medium'

    input:
        tuple val(meta), path(snpVcfIndex), path(snpVcf), path(indelVcfIndex), path(indelVcf)

    output:
        tuple val(meta), path("*.filtered_SNP.RawIndels.vcf.gz.tbi"), path("*.filtered_SNP.RawIndels.vcf.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gatk MergeVcfs --java-options "-Xmx${task.memory.giga}G" \\
        -I ${snpVcf} \\
        -I ${indelVcf} \\
        -O ${meta.id}.filtered_SNP.RawIndels.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.filtered_SNP.RawIndels.vcf.gz
    touch ${meta.id}.filtered_SNP.RawIndels.vcf.gz.tbi
    """
}
