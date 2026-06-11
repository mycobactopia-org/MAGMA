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
include { BWA_MEM as BWA_MEM_DELLY               } from '../modules/local/bwa/mem'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_DELLY } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DELLY  } from '../modules/nf-core/samtools/index/main'
include { GATK4_MARKDUPLICATES as GATK_MARK_DUPLICATES_DELLY   } from '../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR as GATK_BASE_RECALIBRATOR_DELLY } from '../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR as GATK_APPLY_BQSR_DELLY              } from '../modules/nf-core/gatk4/applybqsr/main'
include { DELLY_CALL               } from '../modules/nf-core/delly/call/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_DELLY   } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_DELLY } from '../modules/local/bcftools/merge'
include { BGZIP as BGZIP_MINOR_VARIANTS          } from '../modules/local/bgzip/bgzip'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_LOFREQ } from '../modules/local/bcftools/merge'

// ── TB drug-resistance & structural variant profiles ──────────────────────────
include { TBPROFILER_VCF_PROFILE as TBPROFILER_VCF_PROFILE_DELLY  } from '../modules/local/tbprofiler/vcf_profile'
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE_DELLY           } from '../modules/local/tbprofiler/collate'
include { TBPROFILER_VCF_PROFILE as TBPROFILER_VCF_PROFILE_LOFREQ } from '../modules/local/tbprofiler/vcf_profile'
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE_LOFREQ          } from '../modules/local/tbprofiler/collate'
include { UTILS_MULTIPLE_INFECTION_FILTER } from '../modules/local/utils/multiple_infection_filter'

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

        BWA_MEM_DELLY(
            approved_fastqs_ch,
            params.ref_fasta,
            [params.ref_fasta_dict, params.ref_fasta_amb, params.ref_fasta_ann,
             params.ref_fasta_bwt, params.ref_fasta_fai, params.ref_fasta_pac,
             params.ref_fasta_sa]
        )

        delly_normalize_ch = BWA_MEM_DELLY.out
            .map { meta, bam ->
                def splitted = meta.id.split("\\.")
                def sampleId = splitted[0] + '.' + splitted[1]
                [ meta + [id: sampleId], bam ]
            }
            .groupTuple()

        // Join with approved samples to filter.
        // Shapes after join:
        //   all_samples_ch.map { [it[0], it[0]] }                       → [id, id]
        //   delly_normalize_ch.map { meta, bam -> [meta.id, meta, bam] } → [id, meta, [bams]]
        //   join (key=id)                                                → [id, id, meta, [bams]]
        // The previous 3-arg destructure dropped the second id, causing
        // 'Invalid method invocation `call` with arguments: [id, id, meta, [bam]]'.
        delly_filtered_ch = all_samples_ch
            .map { [ it[0], it[0] ] }
            .join(delly_normalize_ch.map { meta, bam -> [ meta.id, meta, bam ] })
            .map { _id, _id2, meta, bam -> [ meta, bam ] }

        def samtools_merge_delly_input_ch = delly_filtered_ch.map { meta, bams -> [meta, bams, []] }
        SAMTOOLS_MERGE_DELLY(samtools_merge_delly_input_ch, Channel.value([ [id: 'none'], [], [], [] ]))
        GATK_MARK_DUPLICATES_DELLY(SAMTOOLS_MERGE_DELLY.out.bam, [], [])

        if (!params.skip_base_recalibration) {
            def base_recal_delly_input_ch = GATK_MARK_DUPLICATES_DELLY.out.bam.map { meta, bam -> [ meta, bam, [], [] ] }
            GATK_BASE_RECALIBRATOR_DELLY(
                base_recal_delly_input_ch,
                Channel.value([ [id:'ref'],   file(params.ref_fasta)      ]),
                Channel.value([ [id:'ref'],   file(params.ref_fasta_fai)  ]),
                Channel.value([ [id:'ref'],   file(params.ref_fasta_dict) ]),
                Channel.value([ [id:'dbsnp'], file(params.dbsnp_vcf)      ]),
                Channel.value([ [id:'dbsnp'], file(params.dbsnp_vcf_tbi)  ])
            )
            def apply_bqsr_delly_input_ch = GATK_BASE_RECALIBRATOR_DELLY.out.table
                .join(GATK_MARK_DUPLICATES_DELLY.out.bam)
                .map { meta, table, bam -> [ meta, bam, [], table, [] ] }
            GATK_APPLY_BQSR_DELLY(
                apply_bqsr_delly_input_ch,
                file(params.ref_fasta),
                file(params.ref_fasta_fai),
                file(params.ref_fasta_dict)
            )
            recal_delly_ch = GATK_APPLY_BQSR_DELLY.out.bam
        } else {
            recal_delly_ch = GATK_MARK_DUPLICATES_DELLY.out.bam
        }

        SAMTOOLS_INDEX_DELLY(recal_delly_ch)
        // Adapt to the nf-core DELLY_CALL signature, which takes a wide tuple
        //   (meta, input_bam, input_bai, vcf, vcf_index, exclude_bed)
        // with the last three optional, plus separate fasta + fai value tuples
        // and a 'bcf'/'vcf' suffix selector. MAGMA does discovery-mode DELLY
        // calls (no genotyping VCF, no exclude bed) and wants .bcf output.
        def delly_call_input_ch = SAMTOOLS_INDEX_DELLY.out.index
            .join(recal_delly_ch)
            .map { meta, bai, bam -> [ meta, bam, bai, [], [], [] ] }
        def delly_call_fasta_ch = Channel.value([ [id: 'ref'], file(params.ref_fasta)     ])
        def delly_call_fai_ch   = Channel.value([ [id: 'ref'], file(params.ref_fasta_fai) ])
        DELLY_CALL(delly_call_input_ch, delly_call_fasta_ch, delly_call_fai_ch, 'bcf')
        // nf-core BCFTOOLS_VIEW takes [meta, vcf, index] (the index can be passed
        // empty []); we have DELLY_CALL.out.bcf=[meta, bcf] and .csi=[meta, csi],
        // so join them. Optional regions/targets/samples are all empty.
        def bcftools_view_delly_input_ch = DELLY_CALL.out.bcf.join(DELLY_CALL.out.csi)
        BCFTOOLS_VIEW_DELLY(bcftools_view_delly_input_ch, [], [], [])

        // The local module emitted a 3-tuple [meta, bcf.gz, csi] which the
        // downstream pipeline flattened and filtered. The nf-core module emits
        // .vcf (the filtered file) and .index (the auto-generated csi from
        // --write-index=csi) separately; join+flatten the same way.
        def bcftools_view_delly_out_ch = BCFTOOLS_VIEW_DELLY.out.vcf.join(BCFTOOLS_VIEW_DELLY.out.index)

        delly_vcfs_ch = bcftools_view_delly_out_ch
            .collect()
            .flatten()
            .filter { it instanceof java.nio.file.Path }
            .unique()
            .collect(sort: true)

        delly_vcfs_file = bcftools_view_delly_out_ch
            .collect()
            .flatten()
            .filter { it instanceof java.nio.file.Path && it.getExtension() == "gz" }
            .map { it.name }
            .collectFile(name: "${outdir}/structural_variant_vcfs.txt", newLine: true)

        BCFTOOLS_MERGE_DELLY(
            params.vcf_name,
            'delly',
            delly_vcfs_file,
            delly_vcfs_ch
        )

        def delly_resistanceDb = []
        TBPROFILER_VCF_PROFILE_DELLY(BCFTOOLS_MERGE_DELLY.out, delly_resistanceDb)
        TBPROFILER_COLLATE_DELLY(
            params.vcf_name,
            TBPROFILER_VCF_PROFILE_DELLY.out.collect(),
            delly_resistanceDb
        )

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
                TBPROFILER_COLLATE_DELLY.out.per_sample_results
            )

            UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION(
                UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch,
                MINOR_VARIANTS_ANALYSIS_WF.out.minor_variants_results_ch,
                TBPROFILER_COLLATE_DELLY.out.per_sample_results
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
