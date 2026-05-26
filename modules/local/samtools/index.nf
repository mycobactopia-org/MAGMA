process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.bai"), path(bam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools index ${bam}
    """

    stub:
    """
    touch ${meta.id}.bam.bai
    """
}
