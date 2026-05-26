process BCFTOOLS_MERGE {
    label 'process_medium'

    input:
        val(joint_name)
        val(file_format)
        path(vcfs_file)      // a file listing VCFs (lofreq) or use path("*") glob approach
        path("vcf_files/*")  // staged VCF files

    output:
        tuple val(joint_name), path("*.vcf.gz.csi"), path("*.${file_format}.vcf.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools merge -o ${joint_name}.${file_format}.vcf -l ${vcfs_file}
    bgzip ${joint_name}.${file_format}.vcf
    bcftools index ${joint_name}.${file_format}.vcf.gz
    """

    stub:
    """
    touch ${joint_name}.${file_format}.vcf.gz
    touch ${joint_name}.${file_format}.vcf.gz.csi
    """
}
