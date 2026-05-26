process SNPDISTS {
    tag "${meta.id}"
    label 'process_low'

    input:
        val(prefix)
        tuple val(meta), path(alignmentFasta)

    output:
        tuple val(meta), path("*.snp_dists.tsv")
        path("*snp_dists.tsv"), emit: snp_dists_file

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    snp-dists ${alignmentFasta} -b \\
    > ${meta.id}.${prefix}.snp_dists.tsv
    """

    stub:
    """
    touch ${meta.id}.${prefix}.snp_dists.tsv
    """
}
