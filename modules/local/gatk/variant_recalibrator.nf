process GATK_VARIANT_RECALIBRATOR {
    tag "mode: ${analysisMode} ann: ${annotations}"
    label 'process_high'

    input:
        val(analysisMode)
        val(annotations)
        tuple val(meta), path(variantsIndex), path(variantsVcf)
        val(resourceFilesArg)
        path(resourceFiles)
        path(resourceFileIndexes)
        path(reference)
        path("*")

    output:
        tuple val(meta), path("*.tbi"), path("*.recal.vcf.gz"),              emit: recalVcfTuple
        tuple val(meta), path("*.tranches"),                                  emit: tranchesFile
        path("*.R")
        path("*.model")
        path("*.pdf"), optional: true
        tuple val(meta), path("*${analysisMode}*.command.log"),              emit: annotationsLog

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                   = task.ext.args ?: ''
    def finalResourceFilesArg  = resourceFilesArg ? "--resource:${resourceFilesArg}" : ""
    def annotationSuffix       = task.ext.prefix ? ".${task.ext.prefix}" : ""
    """
    gatk VariantRecalibrator --java-options "-Xmx${task.memory.giga}G" \\
        -R ${reference} \\
        -V ${variantsVcf} \\
        ${finalResourceFilesArg} \\
        ${annotations} \\
        ${args} \\
        -mode ${analysisMode} \\
        --tranches-file ${meta.id}.${analysisMode}${annotationSuffix}.tranches \\
        --rscript-file ${meta.id}.${analysisMode}${annotationSuffix}.R \\
        --output ${meta.id}.${analysisMode}${annotationSuffix}.recal.vcf.gz \\
        --output-model ${meta.id}.${analysisMode}${annotationSuffix}.model \\
        2>${meta.id}.${analysisMode}${annotationSuffix}.command.log

    cp ${meta.id}.${analysisMode}${annotationSuffix}.command.log .command.log
    """

    stub:
    """
    touch ${meta.id}.${analysisMode}.tranches
    touch ${meta.id}.${analysisMode}.R
    touch ${meta.id}.${analysisMode}.recal.vcf.gz
    touch ${meta.id}.${analysisMode}.model
    touch ${meta.id}.${analysisMode}.command.log
    """
}
