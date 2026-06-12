# Plan 2 — *Mycobacterium bovis* Integration

**Status:** code assets ready on `feat/plan2-mbovis`, awaiting integration.
**Base branch:** `nf-core-port` on `mycobactopia-org/MAGMA` (post-Plan-1 / `v0.3.0`).
**Working branch:** `feat/plan2-mbovis` (2 commits ahead of an early scaffold ancestor).
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

## 3. What's missing (the actual Plan-2 work)

The branch was started before Plan 1 merged, so most of the work boils down to **rebasing + finishing the wiring**:

1. **Rebase `feat/plan2-mbovis` onto `nf-core-port`.** The pre-Plan-1 commits will conflict in `subworkflows/local/prepare_cohort_vcf.nf` (was edited heavily during the nf-core port) and `workflows/magma.nf` (now down from 662 → 204 lines, completely restructured). Resolution strategy: discard the old `workflows/magma.nf` diff (organism switch belongs in a config profile now, not the workflow file), keep the `prepare_cohort_vcf.nf` SnpEff parameterization.
2. **Add an `mbovis` profile** in `conf/` selecting:
   - `ref_fasta_basename = "NC_002945v4"`
   - `ref_fasta_dir = "${projectDir}/resources/genome"`
   - `excluded_loci_list = "${projectDir}/resources/regions/mbovis/excluded_loci.list"`
   - `snpeff_genome_id = "mbovis_AF2122"`
   - `dbsnp_vcf` + `coll2018_vcf` — **TODO: build or import M. bovis equivalents.** Without these BQSR + VQSR can't run on *M. bovis* data. Stopgap: set `skip_base_recalibration = true` and skip VQSR / use hard filters.
   - `skip_tbprofiler_fastq = true`, `skip_tbprofiler_vcf = true`, `skip_ntmprofiler = true`, `skip_spotyping = true` (none of these have *M. bovis* support).
3. **Carry the `prepare_cohort_vcf.nf` SnpEff genome-ID parameterization forward** into Plan 1's refactored version. The 18-line tweak from `7e90a04` needs to be replayed onto the new file structure.
4. **Run the subsampled *M. bovis* test dataset end-to-end** (`-profile mbovis,abc_cluster` or `mbovis,test`). Confirm variants come out, SnpEff annotations look correct, region exclusions kick in.
5. **Document the gaps** in the README — no DR profiling, VQSR depends on user-supplied truth set, etc.

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

## 7. Open questions for the M. bovis chat

- **Truth set strategy for VQSR on *M. bovis*?** Build one from public *M. bovis* VCFs (Salvador et al.; SB number databases), or skip VQSR and fall back to GATK hard filters? Both are valid; hard filters are simpler but less defensible at low coverage.
- **Is there a target Ti/Tv for *M. bovis*?** *M. tuberculosis* uses 1.85; *M. bovis* may differ. Quick literature check needed.
- **What test data?** The `scripts/create_test_data.sh` subsampling script exists — needs an *M. bovis* SRA accession to pull from. Candidate datasets in the literature.
- **Does TBprofiler have any *M. bovis* support?** Unknown as of now. If yes, partial DR profiling may be possible. If not, leave that for a future plan.
