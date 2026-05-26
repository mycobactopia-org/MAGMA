process TBPROFILER_VCF_PROFILE {
    tag "${meta.id}"
    label 'process_medium'

    input:
        tuple val(meta), path(mergedVcfIndex), path(mergedVcf)
        path(resistanceDb)

    output:
        path("results/*")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def optionalDb = resistanceDb.name != 'NO_FILE' ? "--db ${resistanceDb.name}" : ""
    """
    bcftools view ${mergedVcf} | sed 's/${params.ref_fasta_basename}/Chromosome/g' > intermediate.vcf

    cat intermediate.vcf | bcftools view -Oz -o intermediate.vcf.gz

    tb-profiler profile \\
        ${optionalDb} \\
        --threads ${task.cpus} \\
        --vcf intermediate.vcf.gz \\
        ${args}
    """

    stub:
    """
    mkdir results
    touch results/${meta.id}.results.json
    """
}
