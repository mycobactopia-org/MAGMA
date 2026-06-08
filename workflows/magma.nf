/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                  } from '../modules/nf-core/fastqc/main'
include { MULTIQC                 } from '../modules/local/multiqc/main'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_magma_pipeline'

// ── Validation ────────────────────────────────────────────────────────────────
include { SAMPLESHEET_VALIDATION   } from '../modules/local/utils/samplesheet_validation'
include { FASTQ_VALIDATOR          } from '../modules/local/fastqvalidator/main'
include { UTILS_FASTQ_COHORT_VALIDATION } from '../modules/local/utils/fastq_cohort_validation'

// ── Mapping ────────────────────────────────────────────────────────────────────
include { BWA_MEM                  } from '../modules/local/bwa/mem'

// ── Per-sample calling ─────────────────────────────────────────────────────────
include { SAMTOOLS_MERGE           } from '../modules/local/samtools/merge'
include { SAMTOOLS_INDEX           } from '../modules/local/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_LOFREQ } from '../modules/local/samtools/index'
include { SAMTOOLS_STATS           } from '../modules/local/samtools/stats'
include { GATK_MARK_DUPLICATES     } from '../modules/local/gatk/mark_duplicates'
include { GATK_BASE_RECALIBRATOR   } from '../modules/local/gatk/base_recalibrator'
include { GATK_APPLY_BQSR          } from '../modules/local/gatk/apply_bqsr'
include { GATK_HAPLOTYPE_CALLER    } from '../modules/local/gatk/haplotype_caller'
include { GATK_COLLECT_WGS_METRICS } from '../modules/local/gatk/collect_wgs_metrics'
include { GATK_FLAG_STAT           } from '../modules/local/gatk/flag_stat'
include { LOFREQ_INDELQUAL         } from '../modules/local/lofreq/indelqual'
include { LOFREQ_CALL              } from '../modules/local/lofreq/call'
include { LOFREQ_CALL_NTM          } from '../modules/local/lofreq/call_ntm'
include { LOFREQ_FILTER            } from '../modules/local/lofreq/filter'
include { UTILS_SAMPLE_STATS       } from '../modules/local/utils/sample_stats'
include { UTILS_COHORT_STATS       } from '../modules/local/utils/cohort_stats'
include { UTILS_REFORMAT_LOFREQ    } from '../modules/local/utils/reformat_lofreq'
include { GATK_INDEX_FEATURE_FILE as GATK_INDEX_FEATURE_FILE_LOFREQ } from '../modules/local/gatk/index_feature_file'
include { BGZIP as BGZIP_LOFREQ    } from '../modules/local/bgzip/bgzip'
include { UTILS_MERGE_COHORT_STATS } from '../modules/local/utils/merge_cohort_stats'

// ── QC (FastQC, NTMprofiler, TBprofiler fastq, SpoTyping) ─────────────────────
include { NTMPROFILER_PROFILE      } from '../modules/local/ntmprofiler/profile'
include { NTMPROFILER_COLLATE      } from '../modules/local/ntmprofiler/collate'
include { TBPROFILER_FASTQ_PROFILE } from '../modules/local/tbprofiler/fastq_profile'
include { TBPROFILER_COLLATE as TBPROFILER_FASTQ_COLLATE } from '../modules/local/tbprofiler/collate'
include { SPOTYPING                } from '../modules/local/spotyping/main'
include { UTILS_CAT_SPOTYPING      } from '../modules/local/utils/cat_spotyping'

// ── Structural variants (DELLY) ────────────────────────────────────────────────
include { BWA_MEM as BWA_MEM_DELLY               } from '../modules/local/bwa/mem'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_DELLY } from '../modules/local/samtools/merge'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DELLY  } from '../modules/local/samtools/index'
include { GATK_MARK_DUPLICATES as GATK_MARK_DUPLICATES_DELLY   } from '../modules/local/gatk/mark_duplicates'
include { GATK_BASE_RECALIBRATOR as GATK_BASE_RECALIBRATOR_DELLY } from '../modules/local/gatk/base_recalibrator'
include { GATK_APPLY_BQSR as GATK_APPLY_BQSR_DELLY              } from '../modules/local/gatk/apply_bqsr'
include { DELLY_CALL               } from '../modules/local/delly/call'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_DELLY   } from '../modules/local/bcftools/view'
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
include { GATK_MERGE_VCFS as GATK_MERGE_VCFS_INC      } from '../modules/local/gatk/merge_vcfs'
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

    // Build a per-FASTQ-file channel from the validated JSON
    fastqs_ch = SAMPLESHEET_VALIDATION.out.validated_samplesheet
        .splitJson()
        .map { row ->
            def meta = [
                id:             row.magma_sample_name,
                bam_rg_string:  row.magma_bam_rg_string,
                paired:         row.R2 ? true : false,
                single_end:     row.R2 ? false : true
            ]
            row.R2 ? [ meta, [row.R1, row.R2] ] : [ meta, [row.R1] ]
        }
        .transpose()   // emit one record per fastq file

    // Validate each FASTQ individually
    FASTQ_VALIDATOR(fastqs_ch, SAMPLESHEET_VALIDATION.out.status)

    // Cohort-level validation (QC cutoffs, contamination check)
    UTILS_FASTQ_COHORT_VALIDATION(
        FASTQ_VALIDATOR.out.fastq_report.collect(),
        SAMPLESHEET_VALIDATION.out.validated_samplesheet
    )

    // Build approved-samples channel from magma_analysis.json
    approved_fastqs_ch = UTILS_FASTQ_COHORT_VALIDATION.out.magma_analysis_json
        .splitJson()
        .filter { it.value.fastqs_approved }
        .map { row ->
            def meta = [
                id:            row.value.magma_sample_name,
                bam_rg_string: row.value.magma_bam_rg_string,
                paired:        row.value.R2 ? true : false,
                single_end:    row.value.R2 ? false : true
            ]
            row.value.R2 ? [ meta, [row.value.R1, row.value.R2] ] : [ meta, [row.value.R1] ]
        }

    // =========================================================================
    // QUALITY CHECK (FastQC, NTMprofiler, TBprofiler FASTQ, SpoTyping)
    // =========================================================================

    FASTQC(approved_fastqs_ch)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, file -> file })

    if (!params.skip_ntmprofiler) {
        NTMPROFILER_PROFILE(approved_fastqs_ch)
        NTMPROFILER_COLLATE(params.vcf_name, NTMPROFILER_PROFILE.out.profile_json.collect())
    }

    if (!params.skip_tbprofiler_fastq) {
        TBPROFILER_FASTQ_PROFILE(approved_fastqs_ch)
        TBPROFILER_FASTQ_COLLATE(
            params.vcf_name,
            TBPROFILER_FASTQ_PROFILE.out.json.map { _meta, j -> j }.collect(),
            []
        )
    }

    if (!params.skip_spotyping) {
        SPOTYPING(approved_fastqs_ch)
        UTILS_CAT_SPOTYPING(SPOTYPING.out.txt.collect())
    }

    // =========================================================================
    // EARLY EXIT: only_validate_fastqs mode
    // =========================================================================

    if (!params.only_validate_fastqs) {

        // =====================================================================
        // MAPPING — BWA (main pipeline)
        // =====================================================================

        BWA_MEM(
            approved_fastqs_ch,
            params.ref_fasta,
            [params.ref_fasta_dict, params.ref_fasta_amb, params.ref_fasta_ann,
             params.ref_fasta_bwt, params.ref_fasta_fai, params.ref_fasta_pac,
             params.ref_fasta_sa]
        )

        // =====================================================================
        // CALL_WF — per-sample variant calling
        // =====================================================================

        // Group by sample (collapse per-library BAMs into one tuple per sample)
        normalize_libraries_ch = BWA_MEM.out
            .map { meta, bam ->
                def splitted = meta.id.split("\\.")
                def sampleId = splitted[0] + '.' + splitted[1]
                [ meta + [id: sampleId], bam ]
            }
            .groupTuple()

        SAMTOOLS_MERGE(normalize_libraries_ch)

        GATK_MARK_DUPLICATES(SAMTOOLS_MERGE.out)

        // BQSR (disabled by default for Mtb)
        if (!params.skip_base_recalibration) {
            GATK_BASE_RECALIBRATOR(
                GATK_MARK_DUPLICATES.out.bam_tuple,
                params.dbsnp_vcf,
                params.ref_fasta,
                [params.ref_fasta_fai, params.ref_fasta_dict, params.dbsnp_vcf_tbi]
            )
            GATK_APPLY_BQSR(
                GATK_BASE_RECALIBRATOR.out,
                params.ref_fasta,
                [params.ref_fasta_fai, params.ref_fasta_dict]
            )
            recalibrated_bam_ch = GATK_APPLY_BQSR.out
        } else {
            recalibrated_bam_ch = GATK_MARK_DUPLICATES.out.bam_tuple
        }

        SAMTOOLS_INDEX(recalibrated_bam_ch)

        // HaplotypeCaller (major variants / GVCF)
        GATK_HAPLOTYPE_CALLER(
            SAMTOOLS_INDEX.out,
            params.ref_fasta,
            [params.ref_fasta_fai, params.ref_fasta_dict]
        )

        // NTM contamination estimate via LoFreq call at 16S locus
        LOFREQ_CALL_NTM(
            SAMTOOLS_INDEX.out,
            params.ref_fasta,
            [params.ref_fasta_fai]
        )

        // LoFreq minor-variant calling
        LOFREQ_INDELQUAL(recalibrated_bam_ch, params.ref_fasta)
        SAMTOOLS_INDEX_LOFREQ(LOFREQ_INDELQUAL.out)
        LOFREQ_CALL(SAMTOOLS_INDEX_LOFREQ.out, params.ref_fasta, [params.ref_fasta_fai])
        LOFREQ_FILTER(LOFREQ_CALL.out, params.ref_fasta)

        // Reformat LoFreq VCFs for downstream merging
        UTILS_REFORMAT_LOFREQ(LOFREQ_CALL.out)
        BGZIP_LOFREQ(UTILS_REFORMAT_LOFREQ.out)
        GATK_INDEX_FEATURE_FILE_LOFREQ(BGZIP_LOFREQ.out)

        // Per-sample QC stats
        SAMTOOLS_STATS(recalibrated_bam_ch, params.ref_fasta)
        GATK_COLLECT_WGS_METRICS(recalibrated_bam_ch, params.ref_fasta)
        GATK_FLAG_STAT(recalibrated_bam_ch, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict])

        sample_stats_ch = SAMTOOLS_STATS.out
            .join(GATK_COLLECT_WGS_METRICS.out)
            .join(GATK_FLAG_STAT.out)
            .join(LOFREQ_CALL_NTM.out)

        UTILS_SAMPLE_STATS(sample_stats_ch)
        UTILS_COHORT_STATS(UTILS_SAMPLE_STATS.out.collect())

        // =====================================================================
        // MINOR VARIANTS ANALYSIS (LoFreq cohort merge + TBprofiler)
        // =====================================================================

        lofreq_vcf_tuple_ch = GATK_INDEX_FEATURE_FILE_LOFREQ.out.vcf_tuple.collect(sort: true)

        vcfs_file = lofreq_vcf_tuple_ch
            .flatten()
            .filter { it instanceof java.nio.file.Path && it.getExtension() == "gz" }
            .map { it.name }
            .collectFile(name: "${outdir}/minor_variant_vcfs.txt", newLine: true)

        BCFTOOLS_MERGE_LOFREQ(
            params.vcf_name,
            'lofreq',
            vcfs_file,
            lofreq_vcf_tuple_ch
                .flatten()
                .filter { it instanceof java.nio.file.Path }
                .collect()
        )

        def resistanceDb = []
        TBPROFILER_VCF_PROFILE_LOFREQ(BCFTOOLS_MERGE_LOFREQ.out, resistanceDb)
        TBPROFILER_COLLATE_LOFREQ(
            params.vcf_name,
            TBPROFILER_VCF_PROFILE_LOFREQ.out.collect(),
            resistanceDb
        )
        UTILS_MULTIPLE_INFECTION_FILTER(TBPROFILER_COLLATE_LOFREQ.out.per_sample_results)

        UTILS_MERGE_COHORT_STATS(
            UTILS_MULTIPLE_INFECTION_FILTER.out.approved_samples,
            UTILS_MULTIPLE_INFECTION_FILTER.out.rejected_samples,
            UTILS_COHORT_STATS.out
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

        SAMTOOLS_MERGE_DELLY(delly_filtered_ch)
        GATK_MARK_DUPLICATES_DELLY(SAMTOOLS_MERGE_DELLY.out)

        if (!params.skip_base_recalibration) {
            GATK_BASE_RECALIBRATOR_DELLY(
                GATK_MARK_DUPLICATES_DELLY.out.bam_tuple,
                params.dbsnp_vcf,
                params.ref_fasta,
                [params.ref_fasta_fai, params.ref_fasta_dict, params.dbsnp_vcf_tbi]
            )
            GATK_APPLY_BQSR_DELLY(
                GATK_BASE_RECALIBRATOR_DELLY.out,
                params.ref_fasta,
                [params.ref_fasta_fai, params.ref_fasta_dict]
            )
            recal_delly_ch = GATK_APPLY_BQSR_DELLY.out
        } else {
            recal_delly_ch = GATK_MARK_DUPLICATES_DELLY.out.bam_tuple
        }

        SAMTOOLS_INDEX_DELLY(recal_delly_ch)
        DELLY_CALL(SAMTOOLS_INDEX_DELLY.out, params.ref_fasta)
        BCFTOOLS_VIEW_DELLY(DELLY_CALL.out)

        delly_vcfs_ch = BCFTOOLS_VIEW_DELLY.out
            .collect()
            .flatten()
            .filter { it instanceof java.nio.file.Path }
            .unique()
            .collect(sort: true)

        delly_vcfs_file = BCFTOOLS_VIEW_DELLY.out
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

        TBPROFILER_VCF_PROFILE_DELLY(BCFTOOLS_MERGE_DELLY.out, resistanceDb)
        TBPROFILER_COLLATE_DELLY(
            params.vcf_name,
            TBPROFILER_VCF_PROFILE_DELLY.out.collect(),
            resistanceDb
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
            collated_gvcfs_ch = GATK_HAPLOTYPE_CALLER.out.gvcf_ch
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

            merge_inc_vcf_ch = SNP_ANALYSIS.out.snp_inc_vcf_ch
                .join(INDEL_ANALYSIS.out.indel_vcf_ch)

            GATK_MERGE_VCFS_INC(merge_inc_vcf_ch)

            MAJOR_VARIANT_ANALYSIS(GATK_MERGE_VCFS_INC.out, lofreq_vcf_tuple_ch)

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
                TBPROFILER_COLLATE_LOFREQ.out.per_sample_results,
                TBPROFILER_COLLATE_DELLY.out.per_sample_results
            )

            UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION(
                UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch,
                TBPROFILER_COLLATE_LOFREQ.out.per_sample_results,
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
