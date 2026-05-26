process LOFREQ_INDELQUAL {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(recalibratedBam)
        path(ref_fasta)

    output:
        tuple val(meta), path("*.dindel.bam")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    lofreq indelqual \\
        -f ${ref_fasta} \\
        --dindel \\
        ${args} \\
        -o ${meta.id}.dindel.bam \\
        ${recalibratedBam}
    """

    stub:
    """
    touch ${meta.id}.dindel.bam
    """
}
