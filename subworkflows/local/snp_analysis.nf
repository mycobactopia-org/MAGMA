include { GATK_SELECT_VARIANTS as GATK_SELECT_VARIANTS_SNP             } from '../../modules/local/gatk/select_variants'
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR_SNP   } from '../../modules/local/gatk/variant_recalibrator'
include { GATK4_APPLYVQSR as GATK_APPLY_VQSR_SNP                       } from '../../modules/nf-core/gatk4/applyvqsr/main'
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

    // Migration: GATK_APPLY_VQSR is now the standard nf-core gatk4/applyvqsr.
    //   Local input: 'SNP', [meta, tbi, vcf, recalTbi, recalVcf, tranches],
    //                ref_fasta, [ref_fasta_fai, ref_fasta_dict]
    //   nf-core    : [meta, vcf, vcf_tbi, recal, recal_index, tranches],
    //                fasta, fai, dict  (note tuple element order swapped:
    //                vcf is index 1 not 2, recal at index 3, recal_index at 4)
    // The 'SNP' analysisMode val is moved into ext.args (-mode SNP) in
    // conf/modules.config; the output filename suffix is preserved via ext.prefix.
    def apply_vqsr_input_ch = vqsr_ch.map { meta, tbi, vcf, recalTbi, recalVcf, tranches ->
        [ meta, vcf, tbi, recalVcf, recalTbi, tranches ]
    }
    GATK_APPLY_VQSR_SNP(
        apply_vqsr_input_ch,
        file(params.ref_fasta),
        file(params.ref_fasta_fai),
        file(params.ref_fasta_dict)
    )

    // Downstream consumers used the local `filteredVcfTuple` emit shape [meta, tbi, vcf].
    // nf-core gatk4/applyvqsr emits vcf / tbi separately; recompose with .join.
    def apply_vqsr_snp_filtered_ch = GATK_APPLY_VQSR_SNP.out.tbi.join(GATK_APPLY_VQSR_SNP.out.vcf)

    // Exclude rRNA loci
    GATK_SELECT_VARIANTS_EXCLUSION_SNP(
        'SNP',
        apply_vqsr_snp_filtered_ch,
        params.rrna_list,
        params.ref_fasta,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    emit:
    snp_exc_vcf_ch = GATK_SELECT_VARIANTS_EXCLUSION_SNP.out   // rRNA-excluded SNP VCF
    snp_inc_vcf_ch = apply_vqsr_snp_filtered_ch // rRNA-included SNP VCF — [meta, tbi, vcf]
}
