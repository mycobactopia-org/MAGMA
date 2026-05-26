include { GATK_SELECT_VARIANTS as GATK_SELECT_VARIANTS_SNP             } from '../../modules/local/gatk/select_variants'
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR_SNP   } from '../../modules/local/gatk/variant_recalibrator'
include { GATK_APPLY_VQSR as GATK_APPLY_VQSR_SNP                       } from '../../modules/local/gatk/apply_vqsr'
include { GATK_SELECT_VARIANTS_EXCLUSION as GATK_SELECT_VARIANTS_EXCLUSION_SNP } from '../../modules/local/gatk/select_variants_exclusion'
include { OPTIMIZE_VARIANT_RECALIBRATION } from './optimize_variant_recalibration'


workflow SNP_ANALYSIS {

    take:
    cohort_vcf_and_index_ch   // channel: [ meta, tbi, vcf ]

    main:

    // Select SNP variants from the joint VCF
    GATK_SELECT_VARIANTS_SNP(
        'SNP',
        'raw',
        cohort_vcf_and_index_ch,
        '',
        [],
        [],
        params.ref_fasta,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    // Build resource file channels for VQSR
    arg_files_ch = Channel.of(
            ["coll2014,known=false,training=true,truth=true,prior=15.0",  file(params.coll2014_vcf),   file(params.coll2014_vcf_tbi)],
            ["coll2018,known=false,training=true,truth=true,prior=15.0",  file(params.coll2018_vcf),   file(params.coll2018_vcf_tbi)],
            ["Napier2020,known=false,training=true,truth=true,prior=15.0",file(params.napier2020_vcf), file(params.napier2020_vcf_tbi)],
            ["Benavente2015,known=true,training=false,truth=false,prior=5.0",file(params.benavente2015_vcf),file(params.benavente2015_vcf_tbi)]
        )
        .ifEmpty([])
        .map { row -> row != [] ? [ "${row[0]} ${row[1].getName()}", row[1], row[2] ] : [] }
        .flatten()

    args_ch = arg_files_ch
        .filter { it instanceof String || it instanceof org.codehaus.groovy.runtime.GStringImpl }
        .reduce { a, b -> "$a --resource:$b" }
        .ifEmpty('')

    resources_files_ch = arg_files_ch
        .filter { it instanceof java.nio.file.Path && it.getExtension() == "gz" }
        .collect()
        .ifEmpty([])

    resources_file_indexes_ch = arg_files_ch
        .filter { it instanceof java.nio.file.Path && it.getExtension() == "tbi" }
        .collect()
        .ifEmpty([])

    if (!params.skip_variant_recalibration) {

        // Optimized VQSR: try 6 annotation combinations and pick the best
        OPTIMIZE_VARIANT_RECALIBRATION(
            'SNP',
            GATK_SELECT_VARIANTS_SNP.out.variantsVcfTuple,
            args_ch,
            resources_files_ch,
            resources_file_indexes_ch
        )

        vqsr_ch = OPTIMIZE_VARIANT_RECALIBRATION.out.optimized_vqsr_ch

    } else {

        // Simple VQSR with default annotations
        GATK_VARIANT_RECALIBRATOR_SNP(
            'SNP',
            '-an DP -an AS_MQ',
            GATK_SELECT_VARIANTS_SNP.out.variantsVcfTuple,
            args_ch,
            resources_files_ch,
            resources_file_indexes_ch,
            params.ref_fasta,
            [params.ref_fasta_fai, params.ref_fasta_dict]
        )

        vqsr_ch = GATK_SELECT_VARIANTS_SNP.out.variantsVcfTuple
            .join(GATK_VARIANT_RECALIBRATOR_SNP.out.recalVcfTuple)
            .join(GATK_VARIANT_RECALIBRATOR_SNP.out.tranchesFile)
    }

    GATK_APPLY_VQSR_SNP(
        'SNP',
        vqsr_ch,
        params.ref_fasta,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    // Exclude rRNA loci
    GATK_SELECT_VARIANTS_EXCLUSION_SNP(
        'SNP',
        GATK_APPLY_VQSR_SNP.out.filteredVcfTuple,
        params.rrna_list,
        params.ref_fasta,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    emit:
    snp_exc_vcf_ch = GATK_SELECT_VARIANTS_EXCLUSION_SNP.out   // rRNA-excluded SNP VCF
    snp_inc_vcf_ch = GATK_APPLY_VQSR_SNP.out.filteredVcfTuple // rRNA-included SNP VCF
}
