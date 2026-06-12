# Plan 3 — Extract XBS Variant Calling as a Generic-Bacterial Subworkflow

**Status:** approved, Phase 0 about to start.
**Goal:** lift the canonical XBS variant-calling spine out of MAGMA into a standalone, reusable Nextflow pipeline + subworkflow that works on any haploid bacterium with a reference + truth sets.
**New repo:** `mycobactopia-org/xbs-variant-calling` (to be created).
**Consumer:** MAGMA on `nf-core-port` (and, later, any other bacterial pipeline — NF_CORE_TBANALYZER, plan-2 *M. bovis*, etc.).
**Reference paper:** Goig et al., *Microbial Genomics*, 2022 — *Comprehensive and accurate genetic variant identification from contaminated and low-coverage Mycobacterium tuberculosis whole genome sequencing data* ([PMC8743552](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743552/)).
**Reference implementation:** [TORCH-Consortium/XBS-variant-calling-core](https://github.com/TORCH-Consortium/XBS-variant-calling-core) — 2 bash scripts.

---

## 1. Why this is the right factoring

XBS-core was *already* extracted from MAGMA once (as 2 bash scripts in the TORCH-Consortium/XBS-variant-calling-core repo, last touched 2022). The current MAGMA codebase duplicates that logic in Nextflow form. Honoring the original architectural split — and modernizing it to Nextflow — gives:

- **Reuse**: any future bacterial pipeline (Plan 2 *M. bovis*, NF_CORE_TBANALYZER, hypothetical *S. aureus* / *E. coli* pipelines) imports one subworkflow rather than re-implementing GATK joint calling.
- **Fidelity**: the standalone repo can be unit-tested against the Goig 2021 bash core for byte-level equivalence, then MAGMA inherits that fidelity for free.
- **Cleaner MAGMA**: MAGMA's `CALL_WF`, `PREPARE_COHORT_VCF`, `SNP_ANALYSIS`, `INDEL_ANALYSIS` shrink down to "compose the XBS subworkflow + add MAGMA-specific layers (LoFreq, DELLY, profilers, phylogeny)".

## 2. Architectural boundaries

| Concern | In XBS subworkflow | Out (stays in MAGMA / other consumer) |
|---|---|---|
| FASTQ → BAM (BWA-MEM) | ✅ | — |
| Library merge (samtools) | ✅ | — |
| MarkDuplicates | ✅ | — |
| BQSR | ✅ (opt-in, `skip_bqsr = true` by default) | — |
| HaplotypeCaller (per-sample GVCF) | ✅ | — |
| CombineGVCFs + GenotypeGVCFs (cohort) | ✅ | — |
| VQSR (SNP + INDEL) | ✅ | — |
| Per-sample QC stats (SamtoolStats, WgsMetrics, FlagStat) | ✅ (emit only, no gating) | — |
| Region exclusion (rRNA, PE/PPE, complex) | ❌ | MAGMA + consumers |
| SnpEff annotation | ❌ | MAGMA + consumers |
| LoFreq minor variant calling | ❌ | MAGMA |
| DELLY structural variants | ❌ | MAGMA |
| TBprofiler / NTMprofiler / SpoTyping | ❌ | MAGMA |
| Phylogeny + clustering | ❌ | MAGMA |
| Multi-infection filter / sample QC scoring | ❌ | MAGMA's `UTILS_COHORT_STATS` + `UTILS_MULTIPLE_INFECTION_FILTER` |

## 3. Genericness contract — what differs per organism

Only 5 things change:

| Parameter | M. tuberculosis (MAGMA default) | What other organisms supply |
|---|---|---|
| `reference` bundle | NC_000962.3 H37Rv + BWA index | their own reference + index |
| `snp_truth_vcf` | Coll2018 (1,500+ lineage SNPs) | a high-confidence species-level callset |
| `indel_truth_vcf` | (MAGMA uses a curated set) | their own |
| `target_titv` | 1.85 | organism-specific or `null` to disable |
| `dbsnp_vcf` (only if BQSR enabled) | Benavente2015 | their own known-sites VCF |

Ploidy stays at 1 for all bacteria.

## 4. Subworkflow interface

```groovy
workflow XBS_VARIANT_CALLING {
    take:
      ch_reads          // [meta(study, sample, library), [r1, r2]]
      ch_reference      // value: [meta, fasta, fai, dict, bwa_index_files…]
      ch_snp_truth      // value: [meta, vcf, tbi]   (required)
      ch_indel_truth    // value: [meta, vcf, tbi]   (required)
      ch_dbsnp          // value: [meta, vcf, tbi]   (empty unless !skip_bqsr)

    main: …

    emit:
      gvcfs              // per-sample GVCFs                              [meta, vcf, tbi]
      sample_qc          // SamtoolStats + WgsMetrics + FlagStat          [meta, paths…]
      raw_variants       // cohort raw VCF post-GenotypeGVCFs             [meta, vcf, tbi]
      snp_filtered       // FilteredSNPs.vcf.gz + tbi                     [meta, vcf, tbi]
      indel_filtered     // FilteredINDELs.vcf.gz + tbi                   [meta, vcf, tbi]
      vqsr_diagnostics   // tranches, recal, .R plots                     [meta, paths…]
      versions           // standard nf-core versions topic               channel
}
```

Display: `XBS_VARIANT_CALLING:<process>` standalone, `MAGMA:XBS_VARIANT_CALLING:<process>` when MAGMA imports it.

## 5. The canonical 12-stage spine

Taken from the 2 XBS-core bash scripts and the Goig 2021 paper, byte-faithful to the published flags:

| # | Stage | Tool | Key flags |
|---|---|---|---|
| 1 | Map | BWA-MEM | `-M -R <rg>`  (`bwa_extra_args = ""` by default; MAGMA passes `-k 100`) |
| 2 | Sort + index | samtools sort, samtools index | `-O BAM` |
| 3 | Library merge | samtools merge | default |
| 4 | Mark duplicates | GATK MarkDuplicates | default |
| 5 | (Optional) BQSR | GATK BaseRecalibrator + ApplyBQSR | `skip_bqsr = true` by default; opt-in |
| 6 | Per-sample QC | samtools stats, GATK CollectWgsMetrics, GATK FlagStat | `--READ_LENGTH 0 --COVERAGE_CAP 10000 --COUNT_UNPAIRED` |
| 7 | Per-sample HC | GATK HaplotypeCaller | `-ploidy 1 -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation --read-filter MappingQualityNotZeroReadFilter` |
| 8 | Cohort combine | GATK CombineGVCFs | `-G StandardAnnotation -G AS_StandardAnnotation` |
| 9 | Cohort joint genotype | GATK GenotypeGVCFs | `--sample-ploidy 1 -G StandardAnnotation -G AS_StandardAnnotation` |
| 10 | Split SNP/INDEL | GATK SelectVariants (×2) | `--select-type-to-include {SNP,INDEL} --remove-unused-alternates --exclude-non-variants` |
| 11 | SNP VQSR | GATK VariantRecalibrator + ApplyVQSR (SNP) | `-AS --resource:5000SNP,known=false,training=true,truth=true,prior=20.0 $snp_truth -an AS_MQRankSum -an AS_QD -an AS_MQ -an DP --max-gaussians 4 -mq-cap 60 --target-titv ${params.target_titv}` → ApplyVQSR: `--ts-filter-level 99.9 -mode SNP` |
| 12 | INDEL VQSR | GATK VariantRecalibrator + ApplyVQSR (INDEL) | `--resource:500INDEL,known=false,training=true,truth=true,prior=20.0 $indel_truth -an MQRankSum -an QD -an DP --max-gaussians 2 -mq-cap 60` → ApplyVQSR: `--lod-score-cutoff 0.0000 -mode INDEL` |

## 6. Phase plan

| Phase | Output | Validation |
|---|---|---|
| **0. Bootstrap repo** | `mycobactopia-org/xbs-variant-calling` initialized from nf-core template v4.0.2 — schema, samplesheet, test profile, CI green on `nextflow lint`. | nf-core lint |
| **1. Lift modules** | `nf-core modules install` for all 12 tools (bwa/mem, samtools/{sort,index,merge,stats}, gatk4/{markduplicates,baserecalibrator,applybqsr,haplotypecaller,combinegvcfs,genotypegvcfs,selectvariants,variantrecalibrator,applyvqsr,collectwgsmetrics,flagstat}). Apply the same patches MAGMA uses (combinegvcfs: add `tbi` emit; selectvariants: add `--exclude-intervals` mode). | `nextflow lint` per module |
| **2. Build `XBS_VARIANT_CALLING`** | The 12-stage spine wired through `ext.args` / `modules.config` to match Goig 2021. Internal split: `XBS_PER_SAMPLE` (map → merge → MarkDup → optional BQSR → HC → QC) and `XBS_COHORT_VQSR` (CombineGVCFs → GenotypeGVCFs → 2× Select+VQSR+Apply). | Unit-test per module's `.command.sh` against XBS-core bash |
| **3. End-to-end MTB validation** | Standalone test profile using MAGMA's 600k subsampled MTB dataset + Coll2018 truth + H37Rv refs. SciVer-style diff vs MAGMA v0.3.0 outputs for `joint.raw_variants`, `joint.filtered_SNP*`, `joint.filtered_INDEL*`. | byte-identical record counts |
| **4. MAGMA integration** | Install xbs-variant-calling as nf-core subworkflow in MAGMA. Replace MAGMA's `CALL_WF` per-sample calling + `PREPARE_COHORT_VCF` + `SNP_ANALYSIS` + `INDEL_ANALYSIS` core with the subworkflow's emits. MAGMA's surrounding layers (LoFreq, DELLY, profilers, region exclusion, GATK_MERGE_VCFS_INC, phylogeny, clustering, reports) untouched. Tag `v0.4.0` on `nf-core-port`. | **MAGMA SciVer stays green vs `v0.2.4-sciver`** |
| **5. Non-MTB demo** | Add a `test_mbovis` profile (riding on Plan 2 assets) or a synthetic non-MTB profile, proving "5 config values swap → works on another organism." | run completes, sensible VCF |

## 7. Repo layout (nf-core conventional)

```
mycobactopia-org/xbs-variant-calling/
├── main.nf                              # entry → XBS_VARIANT_CALLING
├── nextflow.config                      # default params + manifest
├── nextflow_schema.json                 # nf-schema-validated params
├── modules/nf-core/                     # 12 registry modules
├── subworkflows/local/
│   ├── xbs_variant_calling_wf.nf        # composer
│   ├── xbs_per_sample.nf                # stages 1–7
│   └── xbs_cohort_vqsr.nf               # stages 8–12
├── conf/
│   ├── modules.config                   # ext.args reproducing Goig 2021 flags
│   ├── test_mtb.config                  # H37Rv + Coll2018
│   └── test_mbovis.config               # AF2122 + mbovis truth set (Phase 5)
├── assets/
│   └── samplesheet_schema.json          # study, sample, library, r1, r2 (MAGMA-compatible)
├── CITATIONS.md                         # Goig 2021, XBS-core, GATK
└── README.md
```

## 8. Design decisions

1. **Modules come from the nf-core registry**, not copied from MAGMA. Same modules; avoids drift.
2. **`bwa_extra_args` param** so MAGMA can pass `-k 100` (current tuning) without forcing it on others. Default empty → byte-identical to XBS-core bash.
3. **Per-sample QC inside the subworkflow but no gating.** Emit it; let consumers gate.
4. **Versions via the standard nf-core `versions` topic.** MAGMA's `softwareVersionsToYAML` picks them up automatically.
5. **MAGMA samplesheet shape** (`study, sample, library, r1, r2`) — flows through unchanged.
6. **License: MIT.** Attribute Goig 2021 + XBS-core in `CITATIONS.md` and `README.md`.
7. **No SnpEff in v1.** MAGMA does its own annotation downstream; XBS subworkflow stops at unannotated VCFs. Keeps the boundary clean.
8. **VQSR is required, not optional, in v1.** Pipelines without a truth set fall outside the XBS contract. Future v2 could add a hard-filter fallback path.

## 9. Risk register

- **VQSR is finicky at small N.** Goig 2021 used 12+ samples. Test datasets may need a special "small-N" config (lower `--max-gaussians`, skip Ti/Tv targeting). Defer to Phase 3.
- **Truth set availability for non-MTB.** Real genericness depends on the user having one. Phase 5 will either find/build one for a demo organism or document the "build your own truth set" recipe.
- **Phase-4 MAGMA refactor invalidates Nextflow caches.** Same one-time cost as v0.2.3 and v0.2.4 SciVer runs (~29 min on abc-cluster).
- **GATK module patches must travel.** The two patches MAGMA carries (`combinegvcfs` tbi emit, `selectvariants` exclude-intervals mode) need to be re-applied in the new repo. Documented in `modules/nf-core/<tool>/<tool>.diff`.

## 10. Relationship to other plans

- **Plan 1 (nf-core port)**: ✅ done. Plan 3 builds on Plan 1's nf-core module migration — the same registry modules get installed in the new repo.
- **Plan 2 (M. bovis)**: see `PLAN_2_MBOVIS_INTEGRATION.md`. Plan 2 can be implemented standalone *or* via Plan 3's Path B (M. bovis as the worked example of Plan 3's genericness contract). Path B is recommended once Plan 3 Phase 4 lands.
- **NF_CORE_TBANALYZER**: future. When MAGMA is included as a subworkflow of TBANALYZER, processes will display `NF_CORE_TBANALYZER:MAGMA:XBS_VARIANT_CALLING:<process>` — clean naming all the way down the include stack.
