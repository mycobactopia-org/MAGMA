process SAMTOOLS_MERGE {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path("bams/*")

    output:
        tuple val(meta), path("*.sorted_reads.bam")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools merge \\
        -f \\
        ${meta.id}.sorted_reads.bam \\
        bams/* \\
        -@ ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}.sorted_reads.bam
    """
}
