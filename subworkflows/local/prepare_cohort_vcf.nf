include { GATK_COMBINE_GVCFS              } from '../../modules/local/gatk/combine_gvcfs'
include { GATK_GENOTYPE_GVCFS             } from '../../modules/local/gatk/genotype_gvcfs'
include { SNPEFF                          } from '../../modules/local/snpeff/snpeff'
include { BGZIP                           } from '../../modules/local/bgzip/bgzip'
include { GATK_INDEX_FEATURE_FILE as GATK_INDEX_FEATURE_FILE_COHORT } from '../../modules/local/gatk/index_feature_file'


workflow PREPARE_COHORT_VCF {

    take:
    cohort_gvcfs_ch   // channel: [ [meta], path(gvcf), path(tbi) ]

    main:

    // Build the --variant string by collecting all .gz filenames
    gvcfs_string_ch = cohort_gvcfs_ch
        .flatten()
        .filter { it instanceof java.nio.file.Path && it.getExtension() == "gz" }
        .map { it -> file(it).name }
        .reduce { a, b -> "$a --variant $b" }

    // Collect all GVCFs (paths) for staging
    gvcfs_paths_ch = cohort_gvcfs_ch
        .flatten()
        .filter { it instanceof java.nio.file.Path }
        .collect()

    def refExitRifGvcf    = params.use_ref_gvcf ? file(params.ref_gvcf,     checkIfExists: true) : []
    def refExitRifGvcfTbi = params.use_ref_gvcf ? file(params.ref_gvcf_tbi, checkIfExists: true) : []

    GATK_COMBINE_GVCFS(
        params.vcf_name,
        gvcfs_string_ch,
        gvcfs_paths_ch,
        params.ref_fasta,
        refExitRifGvcf,
        refExitRifGvcfTbi,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    // Wrap the joint_name output in a meta map for downstream consistency
    combined_ch = GATK_COMBINE_GVCFS.out
        .map { joint_name, tbi, vcf -> [ [id: joint_name], tbi, vcf ] }

    GATK_GENOTYPE_GVCFS(combined_ch, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict])

    SNPEFF(GATK_GENOTYPE_GVCFS.out, params.ref_fasta)

    BGZIP(SNPEFF.out)

    GATK_INDEX_FEATURE_FILE_COHORT(BGZIP.out)

    emit:
    cohort_vcf_and_index_ch = GATK_INDEX_FEATURE_FILE_COHORT.out.sample_vcf_tuple // [ meta, tbi, vcf ]
}
