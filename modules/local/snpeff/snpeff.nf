process SNPEFF {
    tag "${meta.id}"
    label 'process_medium'

    input:
        tuple val(meta), path(rawJointVariantsFile)
        path(ref_fasta)
        path(snpeff_config)   // staged snpEff.config — or NO_FILE sentinel for Mtb (container DB)
        path(snpeff_db)       // staged DB directory  — or NO_FILE sentinel for Mtb (container DB)

    output:
        tuple val(meta), path("*.annotated.vcf"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def genome_id  = params.snpeff_genome_id
    def ref_base   = params.ref_fasta_basename
    def use_custom = (snpeff_config.name != 'NO_FILE')

    if (use_custom) {
        // M. bovis (or any custom genome): build the DB from staged GFF3 + sequences,
        // then annotate.  We copy the staged directory under data/<genome_id>/ so that
        // snpEff's data.dir = ./data/ setting resolves correctly and the build writes
        // the .bin file locally (not back into the project resources directory).
        """
        mkdir -p data/${genome_id}
        cp -r ${snpeff_db}/* data/${genome_id}/

        snpEff build -c ${snpeff_config} -dataDir data -gff3 -v ${genome_id}

        rename_vcf_chrom.py --vcf ${rawJointVariantsFile} --source ${ref_base} --target Chromosome \\
            | snpEff ${args} -c ${snpeff_config} -dataDir data \\
            | rename_vcf_chrom.py --target ${ref_base} --source Chromosome \\
         > ${meta.id}.raw_variants.annotated.vcf
        """
    } else {
        // Mtb default: DB is pre-installed in the container — no config staging needed.
        """
        rename_vcf_chrom.py --vcf ${rawJointVariantsFile} --source ${ref_base} --target Chromosome \\
            | snpEff ${args} \\
            | rename_vcf_chrom.py --target ${ref_base} --source Chromosome \\
         > ${meta.id}.raw_variants.annotated.vcf
        """
    }

    stub:
    """
    touch ${meta.id}.raw_variants.annotated.vcf
    """
}
