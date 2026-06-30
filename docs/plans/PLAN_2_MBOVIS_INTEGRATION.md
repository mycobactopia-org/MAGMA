# Plan 2 — *Mycobacterium bovis* Integration

**Status:** 🟡 rebase complete — wiring done, end-to-end test run pending.
**Base branch:** `nf-core-port` on `mycobactopia-org/MAGMA` (post-Plan-1 / `v0.3.0`).
**Working branch:** `feat/plan2-mbovis-v2` — 2 commits ahead of `nf-core-port` (rebased, replaces old `feat/plan2-mbovis`).
**Depends on:** Plan 1 ✅ merged.
**Relationship to Plan 3:** can be implemented standalone (Path A) **or** can ride on top of Plan 3's XBS extraction (Path B — recommended once Plan 3 lands). See "Two implementation paths" below.

---

## 1. Goal

Enable MAGMA to run end-to-end on *M. bovis* sequencing data using the *M. bovis* AF2122/97 reference (NC_002945.4) instead of *M. tuberculosis* H37Rv (NC_000962.3), with organism-appropriate region exclusions, SnpEff annotation, and (eventually) resistance profiling.

Out of scope for v1: *M. bovis* drug resistance database curation (would require collaboration with the TBprofiler maintainers or a dedicated TBprofiler-bovis database). v1 produces variants; resistance interpretation stays an open question.

## 2. Current assets on `feat/plan2-mbovis`

These commits already exist on `origin/feat/plan2-mbovis` (mycobactopia-org/MAGMA):

```
551444f  test: add subsampled test dataset infrastructure
7e90a04  feat: add M. bovis AF2122/97A organism profile (Plan 2)
```

Files added / changed (~352K lines, mostly the FASTA + GenBank reference):

| Path | Purpose |
|---|---|
| `resources/genome/NC_002945v4.{fa,fai,dict,amb,ann,bwt,pac,sa,gb}` | M. bovis AF2122/97 reference + BWA index + GenBank record |
| `resources/regions/mbovis/excluded_loci.list` (129 entries) | M. bovis equivalent of the H37Rv `UVP_List_of_Excluded_loci.list` — PE/PPE, complex, repetitive regions in *M. bovis* coords |
| `resources/regions/mbovis/rRNA.list` (3 entries) | rRNA exclusion intervals in *M. bovis* coords |
| `resources/snpeff/data/mbovis_AF2122/{genes.gff,sequences.fa}` | Custom SnpEff database for *M. bovis* AF2122/97 |
| `resources/snpeff/snpEff.config` | SnpEff config registering the `mbovis_AF2122` database |
| `nextflow.config` (11 lines) | Adds *M. bovis* profile params |
| `nextflow_schema.json` (16 lines) | Schema entries for the new params |
| `subworkflows/local/prepare_cohort_vcf.nf` (18 lines) | Threads through the organism-specific SnpEff genome ID |
| `workflows/magma.nf` (25 lines) | Wires the new profile into the top-level workflow |
| `scripts/create_test_data.sh` (157 lines) | Generates the small subsampled *M. bovis* test dataset |

## 3. Work log — done vs. remaining

### ✅ Done (as of 2026-06-12, on `feat/plan2-mbovis-v2`)

1. **Rebase complete** — `feat/plan2-mbovis-v2` branches from `nf-core-port@eb4abef` (post-Plan-1, v0.3.0). All pre-Plan-1 scaffolding discarded; only the M. bovis-specific assets and wiring carried forward.

2. **`mbovis` profile** — `conf/organisms/mbovis.config` sets `ref_fasta_basename`, `snpeff_genome_id`, `snpeff_config_file`, `snpeff_db_dir`, clears all Mtb truth VCFs, and sets all relevant skip flags (`skip_variant_recalibration`, `skip_spotyping`, `skip_ntmprofiler`, `skip_ismapper`). `conf/organisms/mtb.config` documents Mtb defaults for reference.

3. **SnpEff dual-mode** — `modules/local/snpeff/snpeff.nf` now takes 4 inputs (VCF, ref, `snpeff_config`, `snpeff_db`). For Mtb (`snpeff_config_file = ""`), the NO_FILE sentinel is passed and the container's pre-installed DB is used. For M. bovis, the config/db are staged, the DB is built with `snpEff build -dataDir data`, and annotation follows.

4. **SnpEff channel staging** — `subworkflows/local/prepare_cohort_vcf.nf` resolves `snpeff_config_ch` and `snpeff_db_ch` via NO_FILE sentinel logic and passes them through to `SNPEFF`.

5. **ISMAPPER wired** — `workflows/magma.nf` now imports `ISMAPPER` + `BCFTOOLS_VIEW_ISMAPPER` and calls them under an `if (!params.skip_ismapper)` guard (skipped by default in mbovis.config).

6. **`nextflow.config` cleaned up** — default `snpeff_config_file` and `snpeff_db_dir` set to `""` (use container DB for Mtb). Profile entries `mtb` and `mbovis` added.

7. **`.gitignore` fixed** — `!resources/snpeff/data/` exception allows the M. bovis SnpEff source files to be tracked.

### 🔲 Still to do

1. **Run the subsampled *M. bovis* test dataset end-to-end** (`-profile mbovis,abc_cluster` or `mbovis,test`). An *M. bovis* SRA accession is needed for `scripts/create_test_data.sh`. Confirm variants, SnpEff annotations, and region exclusions.
2. **Document gaps** in the README — no DR profiling, VQSR disabled until a *M. bovis* truth VCF is curated, ISMapper pending IS element query sequences.
3. **Open PR** from `feat/plan2-mbovis-v2` → `nf-core-port` once the end-to-end test passes.

## 4. Two implementation paths

### Path A — Drop-in profile on `nf-core-port` (the original Plan 2)

Add `conf/mbovis.config`, rebase `feat/plan2-mbovis`, run the test dataset, ship as MAGMA v0.5.0. Pros: fast, ships in one PR, no cross-repo coordination. Cons: every future organism (Plan 3+) repeats this same dance; *M. bovis*-specific assets live inside MAGMA.

### Path B — Ride on Plan 3's XBS extraction (recommended after Plan 3 lands)

Plan 3 extracts the GATK joint-calling spine into `mycobactopia-org/xbs-variant-calling` as a generic-bacterial subworkflow that takes ref + truth sets + ploidy as inputs. With that in place, Plan 2 becomes "pass *M. bovis* config values to the XBS subworkflow + provide custom SnpEff data + region exclusions" — meaningfully smaller. Pros: clean separation; the *M. bovis* assets become the worked example proving genericness. Cons: depends on Plan 3 reaching Phase 4 (MAGMA integration) first.

**Recommended sequencing:**
- If Plan 3 is on track for landing in a few weeks → wait, do Plan 2 via Path B.
- If Plan 3 is uncertain or paused → do Plan 2 via Path A now, then refactor Plan 2 to Path B once Plan 3 lands. The assets in `resources/regions/mbovis/`, `resources/snpeff/data/mbovis_AF2122/` carry over verbatim.

## 5. How to start a new session on this

```
cd /Users/abhi/projects/PHD-pub-magma-pipeline/analysis/pipelines/TORCH-consortium-magma
git fetch origin
git checkout -b plan2-rebase origin/feat/plan2-mbovis
git rebase origin/nf-core-port    # expect conflicts in workflows/magma.nf and subworkflows/local/prepare_cohort_vcf.nf
```

The new chat should be briefed with: this file, `docs/plans/PLAN_3_XBS_EXTRACTION.md` (so it knows about the Path-B option), and `NFCORE_MIGRATION_TRACKER.md` (so it knows what Plan 1 changed in the codebase the rebase is landing onto).

## 6. Decisions already locked in (from earlier planning sessions)

These were decided before Plan 1 merged and should carry forward into the rebased Plan 2 work unless explicitly revisited:

- **VQSR: skip for *M. bovis* (`skip_variant_recalibration = true`)** — no curated truth VCF for *M. bovis* yet. Fall back to GATK hard filters (or, under Path B, run XBS subworkflow without VQSR via a `skip_vqsr` flag we'd add in Plan 3).
- **ISMapper: skip by default (`skip_ismapper = true`)** — IS900 query infrastructure not built; revisit when query targets are agreed.
- **TBprofiler: use default WHO database** — it already includes *M. bovis* references; the existing chromosome-rename step (`sed 's/${params.ref_fasta_basename}/Chromosome/g'`) is param-driven from Plan 1, so no module changes needed.
- **NTMprofiler: skip** — not meaningful for an MTBC member.
- **SpoTyping: skip** — *M. tuberculosis*-specific.
- **SnpEff strategy: dual-mode.** Default ("") = use the container's pre-installed Mtb database (zero build cost for MAGMA's MTB path). For *M. bovis*, `mbovis.config` sets `snpeff_config_file` → build DB at runtime from the staged `genes.gff` + `sequences.fa` in `resources/snpeff/data/mbovis_AF2122/`.
- **`NO_FILE` sentinel** for optional path inputs — `assets/NO_FILE` is the convention used for "no file supplied"; modules check the basename instead of pathnull.

## 7. Open questions and research notes

### Truth set for GATK VQSR on *M. bovis* (answered 2026-06-30)

**What MAGMA needs for Mtb VQSR (the analogues needed for M. bovis):**

| Mtb input | Purpose | M. bovis analogue |
|---|---|---|
| `coll2018_vcf` | SNP truth/training (1,500+ phylogenetic lineage markers) | Not yet built — see options below |
| `walker2015_vcf` / `napier2020_vcf` | INDEL truth | Not yet built |
| `benavente2015_vcf` | Known sites for BQSR | Not yet built |
| `ref_gvcf` (EXIT-RIF) | Reference GVCF to add diversity for small/clonal cohorts | Not yet built — arguably more important for M. bovis than the truth VCF |

**Plan 2 v1 decision: hard filters (`skip_variant_recalibration = true`)**

GATK's own guidance recommends hard filters when N < 30 samples or no truth VCF is available. For M. bovis both apply. The phylogenetic signal is very strong (highly clonal, low SNP density), so the VQSR Gaussian mixture model adds little over standard hard filters (`QD < 2.0`, `MQ < 40.0`, `FS > 60.0`). This is defensible in any publication.

**Plan 2 v2 (publishable quality): build truth VCFs from public data**

The best public sources for constructing M. bovis truth VCFs (analogous to Coll 2018):

| Source | What it provides | Access |
|---|---|---|
| **APHA/Defra UK cattle TB WGS** | Thousands of UK isolates across multiple clonal complexes; deepest M. bovis WGS resource available | ENA/SRA; accession lists in Crispell et al., Smith et al. |
| **Loiseau et al. 2020** (*Nature Microbiology*) — "An African origin for *M. bovis*" | Global phylogeny across all major clonal complexes; published SNP matrix | SRA + Supplementary |
| **Crispell et al. 2019, 2021** (Sanger/APHA) | Cattle–badger transmission, ~thousands of UK isolates, high-confidence SNP calls | ENA, whole datasets deposited |
| **New Zealand AgResearch datasets** | NZ TB eradication program; distinct clonal complex | SRA |
| **Biek et al. 2012** (*Science*) | Foundational WGS paper, SNP-level data | Available |

**Build recipe for the M. bovis truth VCF:**

1. Pull ~100–200 M. bovis isolates spanning ≥ 3 clonal complexes (UK1, UK2, European, African) from the above sources.
2. Call with GATK HaplotypeCaller + MAGMA's `-profile mbovis` pipeline (same pipeline we are building), using NC_002945v4 as the reference.
3. Hard-filter for high-confidence sites: `QUAL > 200`, `DP > 20`, `AF = 1.0` (population-level fixed calls).
4. Cross-validate with a second caller (FreeBayes or `bcftools mpileup`).
5. Intersect → `mbovis_truth_snp.vcf.gz` and `mbovis_truth_indel.vcf.gz`.
6. Add to `resources/truth/` and enable in `mbovis.config`.

**Reference GVCF (EXIT-RIF analogue) — arguably the higher priority:**

EXIT-RIF adds lineage diversity when the study cohort is small or clonal. For M. bovis (often < 20 samples per study, and highly clonal within a herd outbreak), this matters more than it does for Mtb. Build from the same APHA/Loiseau isolates:
- Call 15–20 diverse isolates spanning different clonal complexes as gVCFs with MAGMA `-profile mbovis`
- Store as `resources/ref_gvcfs/Mbovis_ClonalComplexes.g.vcf.gz`
- Re-enable `use_ref_gvcf = true` in `mbovis.config` once available

**Config stubs for when the above are built** (currently commented out in `mbovis.config`):

```groovy
// VQSR truth VCFs — uncomment and set skip_variant_recalibration = false
// coll_mbovis_vcf     = "${projectDir}/resources/truth/mbovis_truth_snp.vcf.gz"
// coll_mbovis_vcf_tbi = "${params.coll_mbovis_vcf}.tbi"

// Reference GVCF for small/clonal M. bovis cohorts (EXIT-RIF analogue)
// use_ref_gvcf = true
// ref_gvcf     = "${projectDir}/resources/ref_gvcfs/Mbovis_ClonalComplexes.g.vcf.gz"
// ref_gvcf_tbi = "${params.ref_gvcf}.tbi"
```

**Ti/Tv target for VQSR:** Published M. bovis WGS studies report Ti/Tv ≈ 1.7–1.9 (similar range to Mtb's 1.85). Confirm from the actual truth-set data before setting `--target-titv`. Not needed until `skip_variant_recalibration = false`.

---

### Other open questions

- **Test data SRA accession:** `scripts/create_test_data.sh` exists but needs an *M. bovis* SRA accession. Candidate: any well-characterised APHA isolate from Crispell et al. / Loiseau et al. supplementary tables.
- **TBprofiler *M. bovis* support:** Unknown. If the WHO TB database includes M. bovis drug resistance markers, partial DR profiling may be possible without changes to the pipeline. Check `tb-profiler version` and the bundled DB metadata.
- **ISMapper queries for M. bovis:** IS900 and IS1081 are M. bovis-specific IS elements. Query sequences need to be assembled and added to `resources/regions/` before `skip_ismapper = false` can be used.
