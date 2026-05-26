process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(bam)
        path(reference)

    output:
        tuple val(meta), path("*.SamtoolStats.txt")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools stats \\
        ${args} \\
        ${bam} \\
        -r ${reference} \\
    > ${meta.id}.SamtoolStats.txt
    """

    stub:
    """
    touch ${meta.id}.SamtoolStats.txt
    """
}
