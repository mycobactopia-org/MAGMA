# nf-core Module Migration Tracker (Plan 1)

Branch: `feat/plan1-implement-magma`.
Reference: `analysis/pipelines/torch-magma/NFCORE_MODULE_MIGRATION_ASSESSMENT.md`.

**Current state:** code complete, SciVer green at `v0.2.4-sciver` (commit `64de4b3`).
Ready to merge into `main`.

## Pattern established (POC: SNPDISTS, commit d341654)

For each Tier A module the migration is:

1. `nf-core modules install <tool>` (requires `JAVA_HOME` set — use the SDKMAN java 21).
2. Inspect `modules/nf-core/<tool>/main.nf` for the script template + input shape.
3. If the local caller passes a 2-channel `(prefix_ch, tuple val(meta), …)` shape
   that the standard nf-core single-channel input can't accept directly, add a
   small **adapter** at the call site — `prefix_ch.combine(other_ch).map { … }`
   — to fold extra context into `meta` (e.g. `meta + [phylo_prefix: prefix]`).
   This avoids refactoring the whole subworkflow until all sibling modules are
   ready to migrate at once.
4. In `conf/modules.config`, set `ext.args` + `ext.prefix` so the emitted
   `.command.sh` matches the torch-magma reference (modulo cosmetic arg order).
5. Update emit consumers if the local module's emit shape differs from
   `nf-core <tool>.out.<emit>`.
6. `git rm` the local module + empty directory.
7. `nextflow lint` — must pass on the changed files.
8. Commit one swap per module.

## Verification

After each commit, the byte-identity check is:
- Run the pipeline through the migrated process (cheap with `-resume`).
- Diff the emitted `.command.sh` against the cached torch-magma reference
  pool at `/tmp/cmp/tm/` (see `/tmp/contentdiff.sh` from the prior session).
- Accept only "byte-identical" or "cosmetic arg-order divergence with
  identical output".

End-to-end validation is **SciVer** against `magma-subsampled-test-v8`
(torch-magma baseline). Tags `v0.2.0-sciver` → `v0.2.4-sciver` document
successive green checkpoints; current SciVer surface (v0.2.4):

- **Byte-identical**: 2 snp_distances TSVs, joint.merged_cohort_stats.tsv,
  6 TBprofiler itol summaries (major + minor).
- **Topology-identical** (tip-set): 2 phylogeny treefiles + 4 cluster
  picks (5SNP/12SNP × Ex/Inc).
- **Record-count match**: 5 of 6 cohort VCFs (SNP RawIndels, lofreq,
  SNP exc/inc-rRNA, raw annotated).
- **Known divergence**: delly cohort VCF (41 vs 39) + the 3 delly
  TBprofiler itol summaries that derive from it. Root cause: DELLY 1.7.3
  vs the bundled DELLY in `magma-container-1`. Tracked, accepted, not
  blocking.

## Status — modules migrated

All Tier A modules migrated to nf-core registry; aliases in parens.

| Module (nf-core path) | Aliases | Commit |
|---|---|---|
| `snpdists`                       | 1     | `9c4d846` |
| `samtools/stats`                 | 1     | `e5491db` |
| `samtools/index`                 | 3 (`_DELLY`, `_LOFREQ`) | `2eff953` |
| `samtools/merge`                 | 2 (`_DELLY`) | `7cc8ef3` |
| `lofreq/indelqual`               | 1     | `8726efd` |
| `lofreq/filter`                  | 1     | `6dacd7e` |
| `delly/call`                     | 1     | `e15ee0d` |
| `bcftools/view`                  | 1 (`_DELLY`) | `7220e3c` |
| `gatk4/index_feature_file`       | 2 (`_LOFREQ`) | `dc3af24` |
| `gatk4/merge_vcfs`               | 1 (`_INC`) | `c56fa18` |
| `gatk4/apply_vqsr`               | 1 (`_SNP`) | `295889c` |
| `gatk4/mark_duplicates`          | 2 (`_DELLY`) | `f8ba380` |
| `gatk4/genotype_gvcfs`           | 1     | `7958159` |
| `gatk4/base_recalibrator`        | 2 (`_DELLY`) | `a8a5143` |
| `gatk4/apply_bqsr`               | 2 (`_DELLY`) | `a8a5143` |
| `gatk4/haplotype_caller`         | 1     | `aeef5f4` |
| `gatk4/select_variants`          | 2 of 5 (`_SNP`, `_INDEL`) | `91a0c72` |
| `gatk4/combine_gvcfs` (patched)  | 1 — adds `tbi` emit | `809c3b1` (+ `def7455` adapter) |
| `gatk4/select_variants` (patched) | EXCLUSION_SNP/INDEL — `ext.intervals_mode='exclude'` | `84aa3e6` |
| `gatk4/variantstotable` + new local `UTILS_VARIANT_TABLE_TO_FASTA` | 1 | `e2ee6fd` |

In addition: **21 remaining local modules** were refactored to nf-core
style with named `emit:` blocks (commit `92b319f`) so they compose
identically with the migrated nf-core modules.

## Deferred (intentionally still local)

| Module | Reason for deferral |
|---|---|
| `gatk4/variantrecalibrator` (ANN2) | Local emits use `--output-model`, `--rscript-file`, and capture `.command.log` for downstream parsing — none of which the nf-core wrapper exposes. Migration is invasive. |
| `snpsites` | nf-core variant has no `meta` in its input tuple; would force a phylogeny-subworkflow restructure. Defer until phylogeny rewrite. |
| `gatk4/select_variants` PHYLOGENY alias | Multi-file `-XL` with paired tbi staging; nf-core wrapper expects a single intervals file. Would need either a patched wrapper or a small adapter that pre-concats the exclusion lists. |

All three are isolated; each can be picked up independently after Plan 1
lands. The non-migrated `modules/local/<tool>/` directories follow the
same nf-core-style `emit:` conventions as the migrated ones, so the
swap-when-ready cost stays small.

## Subworkflow layout (torch-magma alignment)

`workflows/magma.nf` (~204 lines) composes 8 subworkflows that mirror
torch-magma's `workflows/*_wf.nf` 1:1:

| Subworkflow file | Process group display |
|---|---|
| `subworkflows/local/validate_fastqs_wf.nf`              | `MAGMA:VALIDATE_FASTQS_WF:*` |
| `subworkflows/local/quality_check_wf.nf`                | `MAGMA:QUALITY_CHECK_WF:*` |
| `subworkflows/local/map_wf.nf`                          | `MAGMA:MAP_WF:BWA_MEM` |
| `subworkflows/local/call_wf.nf`                         | `MAGMA:CALL_WF:*` |
| `subworkflows/local/minor_variants_analysis_wf.nf`      | `MAGMA:MINOR_VARIANTS_ANALYSIS_WF:*` |
| `subworkflows/local/structural_variants_analysis_wf.nf` | `MAGMA:STRUCTURAL_VARIANTS_ANALYSIS_WF:*` |
| `subworkflows/local/merge_wf.nf`                        | `MAGMA:MERGE_WF:*` (includes existing nested `PHYLOGENY_ANALYSIS_*` / `CLUSTER_ANALYSIS_*` / `SNP_ANALYSIS:OPTIMIZE_VARIANT_RECALIBRATION:*`) |
| `subworkflows/local/reports_wf.nf`                      | `MAGMA:REPORTS_WF:MULTIQC` |

When `MAGMA` is later included as a subworkflow of `NF_CORE_TBANALYZER`,
processes will display as `NF_CORE_TBANALYZER:MAGMA:<WF>:<process>`.

## Container strategy after migration (preserved, not disrupted)

The port already uses the idiomatic nf-core pattern: `withName: '<pattern>' { container = '…' }`
blocks in `conf/docker.config`, `conf/singularity.config`, `conf/apptainer.config`,
`conf/podman.config`, and `conf/abc_cluster.config`. **These patterns match by
process alias name, not by module file path**, so they keep working unchanged
when a module moves from `modules/local/` to `modules/nf-core/`.

Proof in-repo: BWA_MEM is already an nf-core module and `conf/docker.config`
routes it into `magma-container-2` via `withName: 'BWA.*|…' { container = … }`.
This is exactly how the rest of Tier A behaves post-migration — emitted
commands stay byte-identical AND the user keeps the bundled-container option.

After Tier A the user has three container modes:
1. **Bundled MAGMA containers** (current default) — any of the existing
   `-profile docker / singularity / apptainer / podman / abc_cluster` profiles.
   No change.
2. **Per-process biocontainers** (new, free) — the nf-core module's baked-in
   `quay.io/biocontainers/<tool>:<version>` defaults kick in when no profile
   overrides container. Suits laptop / Tower / nf-core launch usage.
3. **Mix-and-match** — user writes a config that applies `withName` for some
   tool families and leaves others on biocontainers.

The two container overrides added on `abc_cluster.config` after Tier A
finished (DELLY 1.7.3, bcftools 1.21) are surgical responses to upstream
nf-core modules requiring `--threads` (which the bundled DELLY didn't
support). They're documented at the override site.

## Gotchas

- `nf-core modules install` auto-creates `conf/containers_*.config` files
  (one per OS/arch/runtime). They define `withName: 'FASTQC' { container = '…' }`
  per probe, so if loaded they would **shadow the bundled-container
  profiles** that route tools into `magma-container-1/2`. Treat them as
  scratch: do **not** register as profiles and do **not** commit. (The
  nf-core module's in-file container default is enough fallback for
  mode 2 above.) Same applies to the scratch `magma-results/` dir and any
  `modules/nf-core/<tool>/` dirs from probe installs of deferred modules.
- Lint shows pre-existing deprecation warnings in unrelated modules (stub-block
  variable scope on `tbprofiler/collate.nf` and `utils/eliminate_annotation.nf`).
  Not migration concerns; ignore.
- The POC's `ext.prefix` closure references `meta.phylo_prefix` — a custom key
  we added at the adapter. This pattern (custom meta keys + `ext.prefix` closure)
  is what lets us drive standard nf-core modules without changing the wider
  workflow shape.
- **Workflow renames invalidate the Nextflow task hash** (process path is in
  the hash) — each refactor of subworkflow layout costs one full re-execute.
  Observed twice (v0.2.3 flatten, v0.2.4 8-WF split): ~29 min each on
  abc-cluster with the 600k subsampled dataset.
- **Refactor that mass-replaces emit-name patterns must scope to the
  `output:` block only.** An early sweep mangled `input:` sections;
  always anchor the regex on the output block (8-space indent for the
  module body).

## Next milestones (after Plan 1 merges)

1. **Plan 1 PR → main** (this PR). After merge: tag `v0.3.0`.
2. **NF_CORE_TBANALYZER integration** — include `MAGMA` as a subworkflow
   in `NF_CORE_TBANALYZER`. Display becomes
   `NF_CORE_TBANALYZER:MAGMA:<WF>:<process>`. No changes to MAGMA itself.
3. **Plan 2 — M. bovis integration** (`feat/plan2-mbovis`, 8 stages:
   organism profile, `mbovis.config`, SnpEff param, resource files).
   Blocked on Plan 1 merge.
4. **Pick up deferred modules** (variantrecalibrator, snpsites,
   select_variants PHYLOGENY) opportunistically — each is independent.
