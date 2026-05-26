include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR_ANN7 } from '../../modules/local/gatk/variant_recalibrator'
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR_ANN6 } from '../../modules/local/gatk/variant_recalibrator'
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR_ANN5 } from '../../modules/local/gatk/variant_recalibrator'
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR_ANN4 } from '../../modules/local/gatk/variant_recalibrator'
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR_ANN3 } from '../../modules/local/gatk/variant_recalibrator'
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR_ANN2 } from '../../modules/local/gatk/variant_recalibrator'
include { UTILS_ELIMINATE_ANNOTATION as UTILS_ELIMINATE_ANNOTATION_ANN7 } from '../../modules/local/utils/eliminate_annotation'
include { UTILS_ELIMINATE_ANNOTATION as UTILS_ELIMINATE_ANNOTATION_ANN6 } from '../../modules/local/utils/eliminate_annotation'
include { UTILS_ELIMINATE_ANNOTATION as UTILS_ELIMINATE_ANNOTATION_ANN5 } from '../../modules/local/utils/eliminate_annotation'
include { UTILS_ELIMINATE_ANNOTATION as UTILS_ELIMINATE_ANNOTATION_ANN4 } from '../../modules/local/utils/eliminate_annotation'
include { UTILS_ELIMINATE_ANNOTATION as UTILS_ELIMINATE_ANNOTATION_ANN3 } from '../../modules/local/utils/eliminate_annotation'
include { UTILS_ELIMINATE_ANNOTATION as UTILS_ELIMINATE_ANNOTATION_ANN2 } from '../../modules/local/utils/eliminate_annotation'
include { UTILS_SELECT_BEST_ANNOTATIONS } from '../../modules/local/utils/select_best_annotations'


workflow OPTIMIZE_VARIANT_RECALIBRATION {

    take:
    analysisType                  // val: "SNP" or "INDEL"
    select_variants_vcftuple_ch   // channel: [ meta, tbi, vcf ]
    args_ch                       // channel: val(resourceFilesArg string)
    resources_files_ch            // channel: path(vcf files)
    resources_file_indexes_ch     // channel: path(tbi files)

    main:

    // ANN7 — all 7 annotations
    GATK_VARIANT_RECALIBRATOR_ANN7(
        analysisType,
        '-an DP -an AS_QD -an AS_MQ -an AS_FS -an AS_SOR -an AS_MQRankSum -an AS_ReadPosRankSum',
        select_variants_vcftuple_ch,
        args_ch,
        resources_files_ch,
        resources_file_indexes_ch,
        params.ref_fasta,
        [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    UTILS_ELIMINATE_ANNOTATION_ANN7(
        analysisType,
        GATK_VARIANT_RECALIBRATOR_ANN7.out.annotationsLog
            .join(GATK_VARIANT_RECALIBRATOR_ANN7.out.tranchesFile)
    )

    // ANN6 — remove lowest-scoring annotation iteratively
    ann6_ch = UTILS_ELIMINATE_ANNOTATION_ANN7.out.reducedAnnotationsFile.map { it.text }

    GATK_VARIANT_RECALIBRATOR_ANN6(
        analysisType, ann6_ch, select_variants_vcftuple_ch,
        args_ch, resources_files_ch, resources_file_indexes_ch,
        params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    UTILS_ELIMINATE_ANNOTATION_ANN6(
        analysisType,
        GATK_VARIANT_RECALIBRATOR_ANN6.out.annotationsLog
            .join(GATK_VARIANT_RECALIBRATOR_ANN6.out.tranchesFile)
    )

    // ANN5
    ann5_ch = UTILS_ELIMINATE_ANNOTATION_ANN6.out.reducedAnnotationsFile.map { it.text }

    GATK_VARIANT_RECALIBRATOR_ANN5(
        analysisType, ann5_ch, select_variants_vcftuple_ch,
        args_ch, resources_files_ch, resources_file_indexes_ch,
        params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    UTILS_ELIMINATE_ANNOTATION_ANN5(
        analysisType,
        GATK_VARIANT_RECALIBRATOR_ANN5.out.annotationsLog
            .join(GATK_VARIANT_RECALIBRATOR_ANN5.out.tranchesFile)
    )

    // ANN4
    ann4_ch = UTILS_ELIMINATE_ANNOTATION_ANN5.out.reducedAnnotationsFile.map { it.text }

    GATK_VARIANT_RECALIBRATOR_ANN4(
        analysisType, ann4_ch, select_variants_vcftuple_ch,
        args_ch, resources_files_ch, resources_file_indexes_ch,
        params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    UTILS_ELIMINATE_ANNOTATION_ANN4(
        analysisType,
        GATK_VARIANT_RECALIBRATOR_ANN4.out.annotationsLog
            .join(GATK_VARIANT_RECALIBRATOR_ANN4.out.tranchesFile)
    )

    // ANN3
    ann3_ch = UTILS_ELIMINATE_ANNOTATION_ANN4.out.reducedAnnotationsFile.map { it.text }

    GATK_VARIANT_RECALIBRATOR_ANN3(
        analysisType, ann3_ch, select_variants_vcftuple_ch,
        args_ch, resources_files_ch, resources_file_indexes_ch,
        params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    UTILS_ELIMINATE_ANNOTATION_ANN3(
        analysisType,
        GATK_VARIANT_RECALIBRATOR_ANN3.out.annotationsLog
            .join(GATK_VARIANT_RECALIBRATOR_ANN3.out.tranchesFile)
    )

    // ANN2
    ann2_ch = UTILS_ELIMINATE_ANNOTATION_ANN3.out.reducedAnnotationsFile.map { it.text }

    GATK_VARIANT_RECALIBRATOR_ANN2(
        analysisType, ann2_ch, select_variants_vcftuple_ch,
        args_ch, resources_files_ch, resources_file_indexes_ch,
        params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict]
    )

    UTILS_ELIMINATE_ANNOTATION_ANN2(
        analysisType,
        GATK_VARIANT_RECALIBRATOR_ANN2.out.annotationsLog
            .join(GATK_VARIANT_RECALIBRATOR_ANN2.out.tranchesFile)
    )

    // Collect all recal VCFs, tranches and annotation JSON files
    recal_vcf_files_ch = GATK_VARIANT_RECALIBRATOR_ANN7.out.recalVcfTuple
        .join(GATK_VARIANT_RECALIBRATOR_ANN6.out.recalVcfTuple)
        .join(GATK_VARIANT_RECALIBRATOR_ANN5.out.recalVcfTuple)
        .join(GATK_VARIANT_RECALIBRATOR_ANN4.out.recalVcfTuple)
        .join(GATK_VARIANT_RECALIBRATOR_ANN3.out.recalVcfTuple)
        .join(GATK_VARIANT_RECALIBRATOR_ANN2.out.recalVcfTuple)
        .flatten()
        .filter { it instanceof java.nio.file.Path }
        .collect()

    tranches_files_ch = GATK_VARIANT_RECALIBRATOR_ANN7.out.tranchesFile
        .join(GATK_VARIANT_RECALIBRATOR_ANN6.out.tranchesFile)
        .join(GATK_VARIANT_RECALIBRATOR_ANN5.out.tranchesFile)
        .join(GATK_VARIANT_RECALIBRATOR_ANN4.out.tranchesFile)
        .join(GATK_VARIANT_RECALIBRATOR_ANN3.out.tranchesFile)
        .join(GATK_VARIANT_RECALIBRATOR_ANN2.out.tranchesFile)
        .flatten()
        .filter { it instanceof java.nio.file.Path }
        .collect()

    annotations_and_tranches_json_files_ch = UTILS_ELIMINATE_ANNOTATION_ANN7.out.annotationsTranchesFile
        .join(UTILS_ELIMINATE_ANNOTATION_ANN6.out.annotationsTranchesFile)
        .join(UTILS_ELIMINATE_ANNOTATION_ANN5.out.annotationsTranchesFile)
        .join(UTILS_ELIMINATE_ANNOTATION_ANN4.out.annotationsTranchesFile)
        .join(UTILS_ELIMINATE_ANNOTATION_ANN3.out.annotationsTranchesFile)
        .join(UTILS_ELIMINATE_ANNOTATION_ANN2.out.annotationsTranchesFile)
        .flatten()
        .filter { it instanceof java.nio.file.Path }
        .collect()

    UTILS_SELECT_BEST_ANNOTATIONS(
        params.vcf_name,
        annotations_and_tranches_json_files_ch,
        recal_vcf_files_ch,
        tranches_files_ch
    )

    emit:
    // optimized_vqsr_ch is the tuple needed by GATK_APPLY_VQSR: [ meta, variants_tbi, variants_vcf, recal_tbi, recal_vcf, tranches ]
    optimized_vqsr_ch = select_variants_vcftuple_ch
        .join(UTILS_SELECT_BEST_ANNOTATIONS.out.bestRecalVcfTuple)
        .join(UTILS_SELECT_BEST_ANNOTATIONS.out.bestTranchesFile)
}
