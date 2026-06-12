# Plan 3 — Extract XBS Variant Calling as a Generic-Bacterial Subworkflow

**Status:** approved, Phase 0 about to start.
**Goal:** lift the canonical XBS variant-calling spine out of MAGMA into a standalone, reusable Nextflow pipeline + subworkflow that works on any haploid bacterium with a reference + truth sets.
**New repo:** `mycobactopia-org/xbs-variant-calling` (to be created).
**Consumer:** MAGMA on `nf-core-port` (and, later, any other bacterial pipeline — NF_CORE_TBANALYZER, plan-2 *M. bovis*, etc.).
**Reference paper:** Heupink et al., *Microbial Genomics*, 2021 — *Comprehensive and accurate genetic variant identification from contaminated and low-coverage Mycobacterium tuberculosis whole genome sequencing data* ([PMC8743552](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743552/)).
**Reference implementation:** [TORCH-Consortium/XBS-variant-calling-core](https://github.com/TORCH-Consortium/XBS-variant-calling-core) — 2 bash scripts.

---

## 1. Why this is the right factoring

XBS-core was *already* extracted from MAGMA once (as 2 bash scripts in the TORCH-Consortium/XBS-variant-calling-core repo, last touched 2022). The current MAGMA codebase duplicates that logic in Nextflow form. Honoring the original architectural split — and modernizing it to Nextflow — gives:

- **Reuse**: any future bacterial pipeline (Plan 2 *M. bovis*, NF_CORE_TBANALYZER, hypothetical *S. aureus* / *E. coli* pipelines) imports one subworkflow rather than re-implementing GATK joint calling.
- **Fidelity**: the standalone repo can be unit-tested against the Heupink 2021 bash core for byte-level equivalence, then MAGMA inherits that fidelity for free.
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

Taken from the 2 XBS-core bash scripts and the Heupink 2021 paper, byte-faithful to the published flags:

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
| **2. Build `XBS_VARIANT_CALLING`** | The 12-stage spine wired through `ext.args` / `modules.config` to match Heupink 2021. Internal split: `XBS_PER_SAMPLE` (map → merge → MarkDup → optional BQSR → HC → QC) and `XBS_COHORT_VQSR` (CombineGVCFs → GenotypeGVCFs → 2× Select+VQSR+Apply). | Unit-test per module's `.command.sh` against XBS-core bash |
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
│   ├── modules.config                   # ext.args reproducing Heupink 2021 flags
│   ├── test_mtb.config                  # H37Rv + Coll2018
│   └── test_mbovis.config               # AF2122 + mbovis truth set (Phase 5)
├── assets/
│   └── samplesheet_schema.json          # study, sample, library, r1, r2 (MAGMA-compatible)
├── CITATIONS.md                         # Heupink 2021, XBS-core, GATK
└── README.md
```

## 8. Design decisions

1. **Modules come from the nf-core registry**, not copied from MAGMA. Same modules; avoids drift.
2. **`bwa_extra_args` param** so MAGMA can pass `-k 100` (current tuning) without forcing it on others. Default empty → byte-identical to XBS-core bash.
3. **Per-sample QC inside the subworkflow but no gating.** Emit it; let consumers gate.
4. **Versions via the standard nf-core `versions` topic.** MAGMA's `softwareVersionsToYAML` picks them up automatically.
5. **MAGMA samplesheet shape** (`study, sample, library, r1, r2`) — flows through unchanged.
6. **License: MIT.** Attribute Heupink 2021 + XBS-core in `CITATIONS.md` and `README.md`.
7. **No SnpEff in v1.** MAGMA does its own annotation downstream; XBS subworkflow stops at unannotated VCFs. Keeps the boundary clean.
8. **VQSR is required, not optional, in v1.** Pipelines without a truth set fall outside the XBS contract. Future v2 could add a hard-filter fallback path.

## 9. Risk register

- **VQSR is finicky at small N.** Heupink 2021 used 12+ samples. Test datasets may need a special "small-N" config (lower `--max-gaussians`, skip Ti/Tv targeting). Defer to Phase 3.
- **Truth set availability for non-MTB.** Real genericness depends on the user having one. Phase 5 will either find/build one for a demo organism or document the "build your own truth set" recipe.
- **Phase-4 MAGMA refactor invalidates Nextflow caches.** Same one-time cost as v0.2.3 and v0.2.4 SciVer runs (~29 min on abc-cluster).
- **GATK module patches must travel.** The two patches MAGMA carries (`combinegvcfs` tbi emit, `selectvariants` exclude-intervals mode) need to be re-applied in the new repo. Documented in `modules/nf-core/<tool>/<tool>.diff`.

## 10. Relationship to other plans

- **Plan 1 (nf-core port)**: ✅ done. Plan 3 builds on Plan 1's nf-core module migration — the same registry modules get installed in the new repo.
- **Plan 2 (M. bovis)**: see `PLAN_2_MBOVIS_INTEGRATION.md`. Plan 2 can be implemented standalone *or* via Plan 3's Path B (M. bovis as the worked example of Plan 3's genericness contract). Path B is recommended once Plan 3 Phase 4 lands.
- **NF_CORE_TBANALYZER**: future. When MAGMA is included as a subworkflow of TBANALYZER, processes will display `NF_CORE_TBANALYZER:MAGMA:XBS_VARIANT_CALLING:<process>` — clean naming all the way down the include stack.

---

## 11. Phase 4 — MAGMA integration (design, awaiting Phase 3 validation)

**Status:** plan only. Implementation blocked on Phase 3 SciVer-style validation of `xbs-variant-calling` standalone (currently waiting on abc-cluster connectivity). When standalone XBS is green, the table below is the playbook.

**Internal subworkflow naming update.** During Phase 2 build, the cohort-side subworkflow was renamed `XBS_COHORT_VQSR` → `XBS_COHORT` (the workflow does CombineGVCFs + GenotypeGVCFs + SelectVariants split + VQSR — the longer name overemphasized one stage). All process display names below reflect this.

### 11.1 Runtime process names after integration

When MAGMA imports the XBS subworkflow, MAGMA's current `MAP_WF`, `CALL_WF` mapping-through-HC, and `MERGE_WF` CombineGVCFs-through-VQSR chains collapse into the call graph below:

```
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:BWA_MEM
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_INDEX_LIB
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_MERGE
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:GATK4_MARKDUPLICATES
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_INDEX_MARKDUP
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:GATK4_BASERECALIBRATOR        ← MAGMA enables BQSR
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:GATK4_APPLYBQSR
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_INDEX_RECAL
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:GATK4_HAPLOTYPECALLER
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_STATS
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_FLAGSTAT
MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:PICARD_COLLECTWGSMETRICS

MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_COMBINEGVCFS
MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_GENOTYPEGVCFS
MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_SELECTVARIANTS_SNP
MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_SELECTVARIANTS_INDEL
MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_VARIANTRECALIBRATOR_SNP
MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_APPLYVQSR_SNP
```

MAGMA's `params.input` integration sets:

| Param | Value | Why |
|---|---|---|
| `bwa_extra_args` | `"-k 100"` | MAGMA's current BWA-MEM tuning |
| `skip_bqsr` | `false` | MAGMA enables BQSR by default (overriding XBS paper-true default) |
| `skip_snp_vqsr` | `false` | MAGMA runs SNP VQSR |
| `skip_indel_vqsr` | `true` | MAGMA uses INDEL hard filters, not VQSR |
| `target_titv` | `1.85` | MTB-specific Ti/Tv from MAGMA |
| `joint_name` | `"joint"` | matches MAGMA's `params.vcf_name` |
| `sample_ploidy` | `1` | bacterial default (unchanged) |

When MAGMA is later wrapped by `NF_CORE_TBANALYZER`, every process gets one extra prefix level: `NF_CORE_TBANALYZER:MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:BWA_MEM`.

### 11.2 Name comparison — old MAGMA process → new MAGMA process

| MAGMA v0.3.0 / nf-core-port | After Phase 4 |
|---|---|
| `MAGMA:MAP_WF:BWA_MEM` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:BWA_MEM` |
| `MAGMA:CALL_WF:SAMTOOLS_MERGE` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_MERGE` |
| `MAGMA:CALL_WF:GATK_MARK_DUPLICATES` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:GATK4_MARKDUPLICATES` |
| `MAGMA:CALL_WF:GATK_BASE_RECALIBRATOR` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:GATK4_BASERECALIBRATOR` |
| `MAGMA:CALL_WF:GATK_APPLY_BQSR` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:GATK4_APPLYBQSR` |
| `MAGMA:CALL_WF:SAMTOOLS_INDEX` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_INDEX_RECAL` |
| `MAGMA:CALL_WF:SAMTOOLS_STATS` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_STATS` |
| `MAGMA:CALL_WF:GATK_COLLECT_WGS_METRICS` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:PICARD_COLLECTWGSMETRICS` ⚠️ |
| `MAGMA:CALL_WF:GATK_FLAG_STAT` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:SAMTOOLS_FLAGSTAT` ⚠️ |
| `MAGMA:CALL_WF:GATK_HAPLOTYPE_CALLER` | `MAGMA:XBS_VARIANT_CALLING:XBS_PER_SAMPLE:GATK4_HAPLOTYPECALLER` |
| `MAGMA:MERGE_WF:PREPARE_COHORT_VCF:GATK_COMBINE_GVCFS` | `MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_COMBINEGVCFS` |
| `MAGMA:MERGE_WF:PREPARE_COHORT_VCF:GATK_GENOTYPE_GVCFS` | `MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_GENOTYPEGVCFS` |
| `MAGMA:MERGE_WF:SNP_ANALYSIS:GATK_SELECT_VARIANTS_SNP` | `MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_SELECTVARIANTS_SNP` |
| `MAGMA:MERGE_WF:SNP_ANALYSIS:GATK_VARIANT_RECALIBRATOR_*` | `MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_VARIANTRECALIBRATOR_SNP` |
| `MAGMA:MERGE_WF:SNP_ANALYSIS:GATK_APPLY_VQSR_SNP` | `MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_APPLYVQSR_SNP` |
| `MAGMA:MERGE_WF:INDEL_ANALYSIS:GATK_SELECT_VARIANTS_INDEL` | `MAGMA:XBS_VARIANT_CALLING:XBS_COHORT:GATK4_SELECTVARIANTS_INDEL` |
| `MAGMA:CALL_WF:LOFREQ_*` | **stays MAGMA-side** — needs new `LOFREQ_PER_SAMPLE_WF` (XBS doesn't do LoFreq) |
| `MAGMA:CALL_WF:UTILS_COHORT_STATS` | **stays MAGMA-side** — MAGMA's cohort QC scoring; consumes XBS's per-sample QC emits |
| `MAGMA:MERGE_WF:SNP_ANALYSIS:GATK_SELECT_VARIANTS_EXCLUSION_SNP` | **stays MAGMA-side** — rRNA exclusion is downstream of XBS |
| `MAGMA:MERGE_WF:INDEL_ANALYSIS:GATK_SELECT_VARIANTS_EXCLUSION_INDEL` | **stays MAGMA-side** — same |
| `MAGMA:MERGE_WF:SNP_ANALYSIS:OPTIMIZE_VARIANT_RECALIBRATION:*` | **DELETED** — XBS's paper-true single-pass VQSR supersedes MAGMA's grid-search variant |

Two module-name renames (⚠️) deserve special attention — they appear in process display *and* in publishDir paths:

- `GATK_COLLECT_WGS_METRICS` → `PICARD_COLLECTWGSMETRICS` (nf-core registers Picard tools under `picard/`)
- `GATK_FLAG_STAT` → `SAMTOOLS_FLAGSTAT` (nf-core has no `gatk4/flagstat`)

Output content is unchanged but file paths shift; MAGMA's `UTILS_COHORT_STATS` parser may grep by filename — needs auditing before Phase 4 lands.

### 11.3 Files in MAGMA that change

| File | Change |
|---|---|
| `workflows/magma.nf` | Replace `MAP_WF` call + `CALL_WF` per-sample chain + `MERGE_WF` CombineGVCFs-through-VQSR chain with one `XBS_VARIANT_CALLING(...)` call. Wire MAGMA-specific channels (`ref_fasta`, `coll2018_vcf`, …) into XBS's input shape. |
| `subworkflows/local/call_wf.nf` | Slim to LoFreq chain + `UTILS_COHORT_STATS`. Or split into a new `lofreq_per_sample_wf.nf` and delete CALL_WF as a unit. |
| `subworkflows/local/map_wf.nf` | **DELETE** — XBS does mapping with `bwa_extra_args="-k 100"`. |
| `subworkflows/local/merge_wf.nf` | `PREPARE_COHORT_VCF` (CombineGVCFs+GenotypeGVCFs+SnpEff+bgzip+index) → SnpEff + bgzip + index of XBS's `raw_variants` emit. Drop SNP_ANALYSIS's VQSR portion (XBS does it). Keep rRNA exclusion in both SNP_ANALYSIS / INDEL_ANALYSIS. Keep `GATK_MERGE_VCFS_INC`. Keep downstream phylogeny / cluster. |
| `subworkflows/local/prepare_cohort_vcf.nf` | Drop CombineGVCFs + GenotypeGVCFs (delegated to XBS); now just SnpEff + bgzip + index. |
| `subworkflows/local/snp_analysis.nf` | Drop the VQSR processes; keep rRNA exclusion. |
| `subworkflows/local/indel_analysis.nf` | Drop the VQSR processes; keep rRNA exclusion. |
| `subworkflows/local/optimize_variant_recalibration.nf` | **DELETE** — MAGMA's grid-search VQSR replaced by XBS's single-pass. |
| `modules/local/*` + `modules/nf-core/gatk4/*` | 16 modules become unused (XBS provides them). Recommend `git rm` at integration time to avoid drift. |
| `modules.json` | Drop entries for the 16 modules. Add an XBS subworkflow entry (mechanism TBD — see §11.4). |
| `conf/modules.config` | Drop `withName` blocks for the 16 XBS-covered modules. |
| `conf/abc_cluster.config` | DELLY / BCFTOOLS overrides stay (unrelated to XBS). If any XBS-covered modules had abc overrides, they move to xbs-variant-calling's `conf/abc_cluster.config`. |
| `NFCORE_MIGRATION_TRACKER.md` | Append "Plan 3 integration" section. |

Net diff estimate: **~3000 lines deleted, ~50 lines added** on the MAGMA side.

### 11.4 How MAGMA pulls in XBS — three mechanism options

| Mechanism | Pros | Cons |
|---|---|---|
| **`nf-core subworkflows install` from custom remote** | Standard tooling; `modules.json` pinning | Requires xbs-variant-calling to expose its subworkflow at `subworkflows/<org>/<name>/main.nf`, not the current `subworkflows/local/`. **Structural change to xbs-variant-calling.** |
| **Git submodule** | No structural change to xbs-variant-calling; commit-SHA pinning; one command to bump version | Submodule UX (`git submodule update --init`); extra step on cloning |
| **Copy-and-pin** | No external infra; pin via `XBS_VERSION` constant + sync script | Manual sync; drift risk; non-idiomatic for nf-core |

**Recommendation for Phase 4 v1: git submodule.** Path becomes:

```
mycobactopia-org/MAGMA/
├── submodules/
│   └── xbs-variant-calling/         # pinned at SHA in .gitmodules
│       └── subworkflows/local/xbs_variant_calling_wf.nf
└── workflows/magma.nf
    include { XBS_VARIANT_CALLING } from '../submodules/xbs-variant-calling/subworkflows/local/xbs_variant_calling_wf'
```

Defer the `nf-core subworkflows install` route to xbs `v0.3.0+` once we restructure for nf-core/subworkflows registry conventions (out of scope for v1 integration).

### 11.5 SciVer impact — what to expect vs `v0.2.4-sciver` baseline

| Output | Risk | Reasoning |
|---|---|---|
| `analyses/snp_distances/joint.ExDR.{Ex,Inc}Complex.snp_dists.tsv` | **LOW** | Downstream of SNP cohort VCF. Byte-identical if SNP VCF unchanged. |
| Phylogeny treefiles + 4 cluster picks | **LOW** | Topology-equivalent if SNP VCF unchanged. |
| `vcf_files/cohort/snp_variant_files/joint.filtered_SNP_exc-rRNA.vcf.gz` record count | **MEDIUM** | MAGMA's `OPTIMIZE_VARIANT_RECALIBRATION` (grid-search over `--max-gaussians`) vs XBS's single-pass at `--max-gaussians 4` may produce different filtered record counts (likely close, not byte-identical). |
| `joint.filtered_SNP_inc-rRNA.vcf.gz` record count | **MEDIUM** | Same. |
| `joint.lofreq.vcf.gz` | **LOW** | LoFreq chain stays MAGMA-side, untouched. |
| `joint.filtered_SNP.RawIndels.vcf.gz` | **MEDIUM** | Goes through XBS's CombineGVCFs + GenotypeGVCFs; should match if ext.args parity holds. |
| `QC_statistics/cohort/joint.merged_cohort_stats.tsv` | **HIGH** | `UTILS_COHORT_STATS` parses per-sample stats. The samtools_stats + flagstat + wgs_metrics filenames change (XBS uses `${meta.id}.SamtoolStats`, `${meta.id}.FlagStat`, `${meta.id}.WgsMetrics`). If the parser greps by filename pattern, it breaks. |
| 3 major_variants TBprofiler itol summaries | **LOW** | Downstream of SNP VCF. |
| 3 minor_variants TBprofiler itol summaries | **LOW** | Unchanged (LoFreq chain). |
| 3 structural_variants TBprofiler itol summaries | **LOW** | Unchanged (DELLY chain). Known divergence stays at 41 vs 39. |

### 11.6 Flag-parity checks needed before Phase 4 lands

These are the most likely sources of unexpected SciVer divergence. Audit each before submitting the Phase 4 PR:

1. **`BWA_MEM`** — MAGMA's local module uses `-k 100`. XBS exposes `bwa_extra_args` for this. MAGMA-side integration **must** set `bwa_extra_args = "-k 100"`.
2. **`GATK_HAPLOTYPE_CALLER`** — MAGMA's local version emits `--bam-output ${meta.id}.haplotype_caller.bam` so downstream MAGMA processes consume the realigned BAM. XBS does not emit the realigned BAM. **Decide:** does MAGMA actually consume it? If yes, add a `haplotypecaller_emit_bam` param to XBS.
3. **`GATK_MARK_DUPLICATES`** — MAGMA's local has specific `ext.prefix` + metrics filename. XBS uses `${meta.id}.markdup.bam`. Filename change may cascade into MAGMA's downstream regex'es.
4. **`GATK_VARIANT_RECALIBRATOR_*`** — MAGMA's `OPTIMIZE_VARIANT_RECALIBRATION` runs VQSR multiple times with different `--max-gaussians` values and picks the best. XBS does single-pass at `--max-gaussians 4`. **This is the biggest behavioral divergence.** Decide: (a) recreate optimize-grid in XBS as a `vqsr_optimize_gaussians = [4, 3, 2, 1]` param, or (b) accept the divergence and update SciVer baselines.
5. **`UTILS_COHORT_STATS`** — MAGMA's parser may pattern-match on the old `GATK_FLAG_STAT` / `GATK_COLLECT_WGS_METRICS` filenames. Update parsing to match XBS's `*.FlagStat`, `*.WgsMetrics`, `*.SamtoolStats`.

### 11.7 Phase 4 deliverables (one PR + one SciVer + one tag)

1. **One PR on `mycobactopia-org/MAGMA` `nf-core-port`** that:
   - Adds `mycobactopia-org/xbs-variant-calling` as a git submodule pinned at the xbs commit that passed Phase 3 standalone validation
   - Refactors `magma.nf`, `call_wf.nf`, `merge_wf.nf`, `prepare_cohort_vcf.nf`, `snp_analysis.nf`, `indel_analysis.nf` per §11.3
   - Deletes `map_wf.nf`, `optimize_variant_recalibration.nf`, 16 modules
   - Updates `UTILS_COHORT_STATS` filename parsing
   - Adds MAGMA-side params (§11.1)
2. **One follow-up SciVer run** against the same MAGMA dataset → produces `v0.4.0-sciver` numbers. Document accepted shifts.
3. **Tag `v0.4.0` on `nf-core-port`** if SciVer is acceptable. (Doesn't need to be byte-identical to `v0.2.4-sciver`; some shifts are expected, especially in the VQSR-derived filtered SNP VCFs from §11.6 item 4.)

### 11.8 Prerequisites — what must be true before Phase 4 starts

- ✅ `xbs-variant-calling` Phase 2 complete (subworkflows written, lint clean, scaffold runnable).
- ⏳ `xbs-variant-calling` Phase 3 complete: end-to-end run on abc-cluster against the 3-sample EXIT-RIF test profile. Expected outputs:
  - per-sample GVCFs landed
  - cohort `raw_variants.vcf.gz` produced
  - `FilteredSNPs.vcf.gz` produced (VQSR converged with `max-gaussians=1`)
  - no fatal errors
  - SciVer-style sanity diff against MAGMA's `v0.3.0` outputs for the equivalent stages (record counts close; full byte-identity not required at this stage)
- ⏳ Decision on the OPTIMIZE_VARIANT_RECALIBRATION question (§11.6 item 4). If recreating the grid in XBS, that's a small follow-up PR on `xbs-variant-calling` first.
- ⏳ Decision on the `--bam-output` question (§11.6 item 2). If MAGMA needs the realigned BAM, add the param to XBS first.

### 11.9 Resumption checklist (for whoever picks Phase 4 up after a break)

1. `git fetch origin && git checkout nf-core-port && git pull`
2. Read this `PLAN_3_XBS_EXTRACTION.md` §11.
3. Check the standalone xbs-variant-calling validation status — most recent tag should be `v0.x.x-validated` or similar.
4. Resolve the two prerequisites in §11.8 with the user before touching MAGMA code.
5. Create a working branch off `nf-core-port`: `git checkout -b feat/xbs-integration`.
6. Walk §11.3 and §11.6 in order, smallest-change-first. Commit after each subworkflow file.
7. Run nextflow lint after each commit; run a local-conda smoke test after the full refactor.
8. Submit to abc-cluster for the v0.4.0-sciver run.
9. Open the PR with the link to this section as the implementation tracker.
