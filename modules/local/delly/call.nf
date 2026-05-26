process DELLY_CALL {
    tag "${meta.id}"
    label 'process_medium'

    input:
        tuple val(meta), path(bai), path(recalibratedBam)
        path(reference)

    output:
        tuple val(meta), path("*.delly.bcf")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    delly call \\
        -g ${reference} \\
        ${recalibratedBam} \\
        ${args} \\
        -o ${meta.id}.delly.bcf
    """

    stub:
    """
    touch ${meta.id}.delly.bcf
    """
}
