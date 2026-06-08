# nf-core Module Migration Tracker (feat/plan1-nfcore-modules)

Worktree branch off `feat/plan1-implement-magma` for the Tier A migration plan.
Reference: `analysis/pipelines/torch-magma/NFCORE_MODULE_MIGRATION_ASSESSMENT.md`.

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

## Status

| # | Module | Status | Commit | Notes |
|---|---|---|---|---|
| 1 | `snpdists` | ✅ MIGRATED | `d341654` | POC. Adapter folds `prefix_ch` into `meta.phylo_prefix`. Filename reproduced exactly via `ext.prefix`. Pending real-run verification once n9 completes. |
| 2 | `iqtree` | pending | — | Same 2-channel pattern as snpdists; same adapter approach. Plus params-driven `if/else` → push into `ext.args` closure (case-on `params.iqtree_*`). |
| 3 | `snpsites` | pending | — | Same 2-channel adapter as snpdists. |
| 4 | `samtools/index` | pending | — | Standalone, no aliases. |
| 5 | `samtools/stats` | pending | — | Standalone. |
| 6 | `delly/call` | pending | — | Standalone; DELLY constraint in `abc_cluster.config` is by alias so stays. |
| 7 | `lofreq/call` | pending | — | nf-core: `lofreq/callparallel`. ext.args carries the per-call flags (NTM has `-r region`). |
| 8 | `lofreq/filter` | pending | — | Direct match. |
| 9 | `lofreq/indelqual` | pending | — | Direct match (we already cleared bogus `-m 60`). |
| 10 | `gatk/combine_gvcfs` | pending | — | `--variant ${input_list}` shape matches. |
| 11 | `gatk/genotype_gvcfs` | pending | — | Direct match. |
| 12 | `gatk/variants_to_table` | pending | — | Direct match. |
| 13 | `gatk/index_feature_file` (×3 aliases) | pending | — | Direct; per-alias `-O` already plumbed via `ext.args`. |
| 14 | `gatk/merge_vcfs` | pending | — | Direct match. |
| 15 | `gatk/apply_vqsr` | pending | — | Direct match. |
| 16 | `gatk/apply_bqsr` | pending | — | Direct match. |
| 17 | `gatk/base_recalibrator` | pending | — | Direct match. |
| 18 | `gatk/variant_recalibrator` (ANN2 alias) | pending | — | Per-mode `ext.args` already in modules.config. |
| 19 | `gatk/select_variants` (×5 aliases) | pending | — | Each alias's args in `ext.args` already. |
| 20 | `gatk/mark_duplicates` (×2 aliases) | pending | — | Direct match. |
| 21 | `bcftools/view` (×2 aliases) | pending | — | Direct match. |
| 22 | `samtools/merge` (×2 aliases) | pending | — | Tier B note: arg shape differs (`bams/*` glob vs file-list) — may need adapter or accept text divergence. |
| 23 | `gatk/haplotype_caller` | LAST | — | gatk4's wrapper has many feature flags (DBSNP, intervals, DRAGSTR, bamout) plumbed via separate vars — careful `ext.args` mapping required. |

## Gotchas

- `nf-core modules install` auto-creates `conf/containers_*.config` files and a
  scratch `magma-results/` dir. They're scaffolding noise and **should not be
  committed** — our container strategy is per-process via `abc_cluster.config`.
  Leave them untracked.
- Lint shows pre-existing deprecation warnings in unrelated modules (stub-block
  variable scope on `tbprofiler/collate.nf` and `utils/eliminate_annotation.nf`).
  Not migration concerns; ignore.
- The POC's `ext.prefix` closure references `meta.phylo_prefix` — a custom key
  we added at the adapter. This pattern (custom meta keys + `ext.prefix` closure)
  is what lets us drive standard nf-core modules without changing the wider
  workflow shape.

## Resume point after n9 validation

Once the validation pipeline completes green on `feat/plan1-implement-magma`,
return here and:
1. Cherry-pick or rebase the migration commits onto the validated base.
2. Run the snpdists migration through the pipeline (one process re-runs).
3. Confirm byte-identity (modulo `-b` arg order) against torch-magma.
4. If green, batch the next 5 trivial migrations (`samtools/index`,
   `samtools/stats`, `delly/call`, `lofreq/filter`, `lofreq/indelqual`,
   `snpsites`) in one sweep — they have no adapter requirements.
