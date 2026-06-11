/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                 } from '../modules/local/multiqc/main'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_magma_pipeline'

// ── Validation ────────────────────────────────────────────────────────────────
include { SAMPLESHEET_VALIDATION   } from '../modules/local/utils/samplesheet_validation'
include { VALIDATE_FASTQS_WF       } from '../subworkflows/local/validate_fastqs_wf'

// ── Mapping ────────────────────────────────────────────────────────────────────
include { MAP_WF                   } from '../subworkflows/local/map_wf'

// ── Per-sample calling + minor variants ────────────────────────────────────────
include { CALL_WF                       } from '../subworkflows/local/call_wf'
include { MINOR_VARIANTS_ANALYSIS_WF    } from '../subworkflows/local/minor_variants_analysis_wf'
include { UTILS_MERGE_COHORT_STATS      } from '../modules/local/utils/merge_cohort_stats'

// ── QC (FastQC, NTMprofiler, TBprofiler fastq, SpoTyping) ─────────────────────
include { QUALITY_CHECK_WF         } from '../subworkflows/local/quality_check_wf'

// ── Structural variants (DELLY) ────────────────────────────────────────────────
include { STRUCTURAL_VARIANTS_ANALYSIS_WF } from '../subworkflows/local/structural_variants_analysis_wf'

// ── Cohort joint-genotyping subworkflows ──────────────────────────────────────
include { PREPARE_COHORT_VCF     } from '../subworkflows/local/prepare_cohort_vcf'
include { SNP_ANALYSIS           } from '../subworkflows/local/snp_analysis'
include { INDEL_ANALYSIS         } from '../subworkflows/local/indel_analysis'
include { MAJOR_VARIANT_ANALYSIS } from '../subworkflows/local/major_variant_analysis'
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS_INCCOMPLEX } from '../subworkflows/local/phylogeny_analysis'
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS_EXCOMPLEX  } from '../subworkflows/local/phylogeny_analysis'
include { CLUSTER_ANALYSIS as CLUSTER_ANALYSIS_INCCOMPLEX    } from '../subworkflows/local/cluster_analysis'
include { CLUSTER_ANALYSIS as CLUSTER_ANALYSIS_EXCOMPLEX     } from '../subworkflows/local/cluster_analysis'

// ── Final merge + reports ──────────────────────────────────────────────────────
include { GATK4_MERGEVCFS as GATK_MERGE_VCFS_INC      } from '../modules/nf-core/gatk4/mergevcfs/main'
include { UTILS_SUMMARIZE_RESISTANCE_RESULTS           } from '../modules/local/utils/summarize_resistance_results'
include { UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION } from '../modules/local/utils/summarize_resistance_results_mixed_infection'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAGMA {

    take:
    ch_samplesheet             // channel: samplesheet read in from --input (nf-core interface)
    multiqc_config             //   param: path to multiqc_config or null
    multiqc_logo               //   param: path to multiqc_logo or null
    multiqc_methods_description //  param: path to methods_description or null
    outdir                     //   param: output directory

    main:

    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    // =========================================================================
    // SAMPLESHEET VALIDATION + FASTQ VALIDATION
    // NOTE: MAGMA uses its own CSV format (study, sample, library, …, r1, r2).
    //       We run SAMPLESHEET_VALIDATION directly on params.input_samplesheet
    //       rather than relying on the nf-schema-parsed ch_samplesheet, which
    //       expects a different schema.
    // =========================================================================

    SAMPLESHEET_VALIDATION(file(params.input_samplesheet ?: params.input, checkIfExists: true))

    VALIDATE_FASTQS_WF(
        SAMPLESHEET_VALIDATION.out.validated_samplesheet,
        SAMPLESHEET_VALIDATION.out.status
    )

    def approved_fastqs_ch = VALIDATE_FASTQS_WF.out.approved_fastqs_ch

    // =========================================================================
    // QUALITY CHECK (FastQC, NTMprofiler, TBprofiler FASTQ, SpoTyping)
    // =========================================================================

    QUALITY_CHECK_WF(approved_fastqs_ch)

    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CHECK_WF.out.reports_fastqc_ch)

    // =========================================================================
    // EARLY EXIT: only_validate_fastqs mode
    // =========================================================================

    if (!params.only_validate_fastqs) {

        // =====================================================================
        // MAPPING — MAP_WF (BWA-MEM)
        // =====================================================================

        MAP_WF(approved_fastqs_ch)

        // =====================================================================
        // CALL_WF — per-sample variant calling + QC stats
        // =====================================================================

        CALL_WF(MAP_WF.out.sorted_reads_ch)

        // gvcf channel reconstructed from CALL_WF's separate tbi+vcf emits.
        // Match torch-magma's gvcf_ch shape: collected stream of [meta, tbi, vcf].
        def gvcf_ch = CALL_WF.out.gvcf_tbi_ch.join(CALL_WF.out.gvcf_vcf_ch).collect()

        def lofreq_vcf_tuple_ch = CALL_WF.out.reformatted_lofreq_vcfs_tuple_ch

        // =====================================================================
        // MINOR VARIANTS ANALYSIS (LoFreq cohort merge + TBprofiler)
        // =====================================================================

        MINOR_VARIANTS_ANALYSIS_WF(lofreq_vcf_tuple_ch, outdir)

        UTILS_MERGE_COHORT_STATS(
            MINOR_VARIANTS_ANALYSIS_WF.out.approved_samples_ch,
            MINOR_VARIANTS_ANALYSIS_WF.out.rejected_samples_ch,
            CALL_WF.out.cohort_stats_tsv
        )

        // MultiQC (MAGMA): the merged cohort stats feed preprocess_multiqc_input.py
        ch_multiqc_files = ch_multiqc_files.mix(UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch)

        // Parse merged cohort stats to identify all samples and approved samples
        all_samples_ch = UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch
            .splitCsv(header: false, skip: 1, sep: '\t')
            .map { row -> [ row.first() ] }

        approved_samples_ch = UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch
            .splitCsv(header: false, skip: 1, sep: '\t')
            .map { row -> [ row.first(), row.last().toInteger() ] }
            .filter { it[1] == 1 }
            .map { [ it[0] ] }

        // =====================================================================
        // STRUCTURAL VARIANTS ANALYSIS (DELLY)
        // =====================================================================

        STRUCTURAL_VARIANTS_ANALYSIS_WF(approved_fastqs_ch, all_samples_ch, outdir)

        // =====================================================================
        // MERGE_WF — cohort joint-genotyping + phylogeny
        // =====================================================================

        if (!params.skip_merge_analysis) {

            // Build gvcf channel keyed on the sample name so it can join the
            // approved-samples channel (which carries sample-name Strings).
            // Port note: GATK_HAPLOTYPE_CALLER.out.gvcf_ch is [meta_MAP, tbi, gvcf];
            // torch-magma's equivalent was [sampleName_STR, tbi, gvcf] (which is
            // why a 3-element collate worked there with a String key). After
            // .collect().flatten() we re-collate to 3 and project the meta map
            // to its id so the .join key matches approved_samples_ch [sampleName].
            //
            // Without this fix the join silently produced an empty channel
            // (meta map ≠ string), the entire merge_wf (PREPARE_COHORT_VCF,
            // SNP_ANALYSIS, INDEL_ANALYSIS, MAJOR_VARIANT_ANALYSIS,
            // PHYLOGENY_*, CLUSTER_*) never fired, and MULTIQC failed at
            // preprocess_multiqc_input.py for missing joint.ExDR.ExComplex.snp_dists.tsv.
            // After the nf-core migration GATK_HAPLOTYPE_CALLER emits .vcf and
            // .tbi as separate channels; rebuild the [meta, tbi, vcf] 3-tuple
            // the existing .flatten/.collate(3) pipeline expects.
            collated_gvcfs_ch = CALL_WF.out.gvcf_tbi_ch
                .join(CALL_WF.out.gvcf_vcf_ch)
                .collect()
                .flatten()
                .collate(3)
                .map { meta, tbi, vcf -> [ meta.id, tbi, vcf ] }

            // Filter to approved samples only
            selected_gvcfs_ch = collated_gvcfs_ch
                .join(approved_samples_ch.map { [it[0], it[0]] })
                .flatten()
                .filter { it instanceof java.nio.file.Path }
                .collect()

            PREPARE_COHORT_VCF(selected_gvcfs_ch)

            SNP_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)
            INDEL_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

            // Build the SNP+INDEL VCF channel for MergeVcfs. The join produces
            // [meta, snp_tbi, snp_vcf, indel_tbi, indel_vcf]; nf-core wants
            // [meta, [vcfs...]] (a list of VCFs to merge).
            def merge_inc_vcf_ch = SNP_ANALYSIS.out.snp_inc_vcf_ch
                .join(INDEL_ANALYSIS.out.indel_vcf_ch)
                .map { meta, _snp_tbi, snp_vcf, _indel_tbi, indel_vcf -> [ meta, [ snp_vcf, indel_vcf ] ] }

            // nf-core gatk4/mergevcfs also takes an (optional) dict value tuple
            // for --SEQUENCE_DICTIONARY; pass an empty tuple — MAGMA's merge
            // doesn't use a dict (torch-magma's command never had one).
            GATK_MERGE_VCFS_INC(merge_inc_vcf_ch, Channel.value([ [id: 'none'], [] ]))

            // Rebuild the downstream [meta, tbi, vcf] tuple that MAJOR_VARIANT_ANALYSIS
            // expects from the two separate nf-core emits (vcf / tbi).
            def merge_inc_vcf_out_ch = GATK_MERGE_VCFS_INC.out.tbi.join(GATK_MERGE_VCFS_INC.out.vcf)

            MAJOR_VARIANT_ANALYSIS(merge_inc_vcf_out_ch, lofreq_vcf_tuple_ch)

            // ExComplex phylogeny (excludes DR loci + complex/repetitive regions)
            excomplex_exclude_interval_ref_ch = Channel.of(
                file(params.coll2018_vcf),
                file(params.coll2018_vcf_tbi),
                file(params.excluded_loci_list)
            ).flatten()

            if (!params.skip_phylogeny_and_clustering) {
                PHYLOGENY_ANALYSIS_EXCOMPLEX(
                    Channel.value('ExDR.ExComplex'),
                    excomplex_exclude_interval_ref_ch,
                    SNP_ANALYSIS.out.snp_exc_vcf_ch
                )

                CLUSTER_ANALYSIS_EXCOMPLEX(
                    PHYLOGENY_ANALYSIS_EXCOMPLEX.out.snpsites_tree_tuple,
                    Channel.value('ExDR.ExComplex')
                )

                // MultiQC (MAGMA): the ExComplex SNP-distance matrix feeds
                // preprocess_multiqc_input.py (joint.ExDR.ExComplex.snp_dists.tsv)
                ch_multiqc_files = ch_multiqc_files.mix(PHYLOGENY_ANALYSIS_EXCOMPLEX.out.snp_dists_ch)
            }

            // IncComplex phylogeny (excludes DR loci only)
            if (!params.skip_complex_regions && !params.skip_phylogeny_and_clustering) {
                inccomplex_exclude_interval_ref_ch = Channel.of(
                    file(params.coll2018_vcf),
                    file(params.coll2018_vcf_tbi)
                ).flatten()

                PHYLOGENY_ANALYSIS_INCCOMPLEX(
                    Channel.value('ExDR.IncComplex'),
                    inccomplex_exclude_interval_ref_ch,
                    SNP_ANALYSIS.out.snp_exc_vcf_ch
                )

                CLUSTER_ANALYSIS_INCCOMPLEX(
                    PHYLOGENY_ANALYSIS_INCCOMPLEX.out.snpsites_tree_tuple,
                    Channel.value('ExDR.IncComplex')
                )
            }

            // Final resistance summary reports
            UTILS_SUMMARIZE_RESISTANCE_RESULTS(
                UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch,
                MAJOR_VARIANT_ANALYSIS.out.major_variants_results_ch,
                MINOR_VARIANTS_ANALYSIS_WF.out.minor_variants_results_ch,
                STRUCTURAL_VARIANTS_ANALYSIS_WF.out.structural_variants_results_ch
            )

            UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION(
                UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch,
                MINOR_VARIANTS_ANALYSIS_WF.out.minor_variants_results_ch,
                STRUCTURAL_VARIANTS_ANALYSIS_WF.out.structural_variants_results_ch
            )

        } // end !skip_merge_analysis
    } // end !only_validate_fastqs

    // =========================================================================
    // Collate and save software versions
    // =========================================================================

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'magma_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )

    // =========================================================================
    // MultiQC — faithful to standard MAGMA (preprocess_multiqc_input.py + multiqc)
    // =========================================================================
    // Note: ch_collated_versions above is still published to pipeline_info. We do
    // NOT mix the nf-core boilerplate (versions / workflow-summary / methods) into
    // the MULTIQC input — standard MAGMA feeds MULTIQC only the FASTQC reports,
    // the merged cohort stats and the ExComplex SNP-distance matrix (mixed in
    // above), matching the reference .command.sh exactly.

    MULTIQC(
        multiqc_config
            ? file(multiqc_config, checkIfExists: true)
            : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
        ch_multiqc_files.flatten().collect()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
