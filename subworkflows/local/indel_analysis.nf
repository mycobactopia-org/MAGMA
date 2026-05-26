include { GATK_SELECT_VARIANTS as GATK_SELECT_VARIANTS_INDEL             } from '../../modules/local/gatk/select_variants'
include { GATK_SELECT_VARIANTS_EXCLUSION as GATK_SELECT_VARIANTS_EXCLUSION_INDEL } from '../../modules/local/gatk/select_variants_exclusion'


// NOTE: Full INDEL VQSR is experimental (XBS_merge#L164) and disabled by default.
//       The indel_vcf_ch emits the raw selected INDELs (no VQSR applied).
workflow INDEL_ANALYSIS {

    take:
    cohort_vcf_and_index_ch   // channel: [ meta, tbi, vcf ]

    main:

    // Select INDEL/MNP/MIXED variants from the joint VCF
    GATK_SELECT_VARIANTS_INDEL(
        'INDEL',
        'raw',
        cohort_vcf_and_index_ch,
        '',
        [],
        [],
        params.ref_fasta,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    // Exclude rRNA loci from indel VCF
    GATK_SELECT_VARIANTS_EXCLUSION_INDEL(
        'INDEL',
        GATK_SELECT_VARIANTS_INDEL.out.variantsVcfTuple,
        params.rrna_list,
        params.ref_fasta,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    emit:
    // NOTE: returns raw selected INDELs (VQSR deferred to future work)
    indel_vcf_ch     = GATK_SELECT_VARIANTS_INDEL.out.variantsVcfTuple
    indel_exc_vcf_ch = GATK_SELECT_VARIANTS_EXCLUSION_INDEL.out
}
