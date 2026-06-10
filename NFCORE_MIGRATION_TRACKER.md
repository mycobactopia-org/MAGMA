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
| 1 | `snpdists` | ✅ MIGRATED | `9c4d846` | POC. Adapter folds `prefix_ch` into `meta.phylo_prefix`. Filename reproduced exactly via `ext.prefix`. Real-run verification pending (next end-to-end test). |
| 5 | `samtools/stats` | ✅ MIGRATED | `e5491db` | Adapter: pad bam-only channel with `[]` for index slot; wrap ref as `[[id:'ref'], fasta, fai]` value tuple. Cosmetic diffs vs torch-magma (added `--threads`, long-form `--reference`, arg order) — same stats output. `ext.prefix = '${meta.id}.SamtoolStats'`. |
| 9 | `lofreq/indelqual` | ✅ MIGRATED | `8726efd` | Adapter: wrap ref as `[[id:'ref'], fasta]` value tuple. `ext.args = '--dindel'`, `ext.prefix = '${meta.id}.dindel'`. Emits same dindel BAM; only flag order changes vs torch-magma. |
| 2 | `iqtree` | pending | — | Same 2-channel pattern as snpdists; same adapter approach. Plus params-driven `if/else` → push into `ext.args` closure (case-on `params.iqtree_*`). |
| 3 | `snpsites` | pending | — | Same 2-channel adapter as snpdists. |
| 4 | `samtools/index` | ✅ MIGRATED | `2eff953` | 3 aliases (`SAMTOOLS_INDEX`, `_DELLY`, `_LOFREQ`). nf-core emits 2-tuple; downstream wants 3-tuple `[meta, bai, bam]`. Resolved with `SAMTOOLS_INDEX.out.index.join(bam_ch)` adapter at each of the 3 call sites. Cosmetic command diff: nf-core adds `-@ ${task.cpus}`. |
| 6 | `delly/call` | ✅ MIGRATED | `e15ee0d` | Most involved adapter so far — pad 3-tuple to 6-tuple with `[]` for unused optionals (vcf/vcf_index/exclude_bed), wrap fasta+fai as separate value tuples, add `'bcf'` suffix selector. Cosmetic command diff: `--genome` (long) vs `-g`, `--outfile` vs `-o`, `--threads` added. |
| 7 | `lofreq/call` | pending | — | nf-core: `lofreq/callparallel`. ext.args carries the per-call flags (NTM has `-r region`). |
| 8 | `lofreq/filter` | ✅ MIGRATED | `6dacd7e` | Direct input match. nf-core writes bgzipped `.vcf.gz`, port wrote uncompressed `.vcf` — change is publish-only (`LOFREQ_FILTER.out` not consumed downstream; UTILS_REFORMAT_LOFREQ takes raw LOFREQ_CALL output). |
| 13 | `gatk/index_feature_file` (×2 active aliases) | ✅ MIGRATED | `dc3af24` | nf-core emits `[meta, tbi]`; downstream wants `[meta, tbi, vcf]` (COHORT) or `[tbi, vcf]` (LOFREQ). Both rebuilt via `.index.join(bgzip_ch)` adapter. Per-alias `ext.args` for `-O <name>.tbi` still applies. Cosmetic command diff: `--input`/`--tmp-dir .`/`-Xmx${mem}M`/`-XX:-UsePerfData` added. |
| 14 | `gatk/merge_vcfs` | ✅ MIGRATED | `c56fa18` | 1 alias (`GATK_MERGE_VCFS_INC`). Adapter: drop unused tbi paths, wrap VCFs into a list `[meta, [snpVcf, indelVcf]]`. Empty dict tuple `[[id:'none'], []]`. Output rebuild: `.out.tbi.join(.out.vcf)` for the downstream `[meta, tbi, vcf]` shape. |
| 15 | `gatk/apply_vqsr` | ✅ MIGRATED | `295889c` | 1 alias (`GATK_APPLY_VQSR_SNP`). Tuple order swap: `[meta, tbi, vcf, recalTbi, recalVcf, tranches]` → `[meta, vcf, tbi, recalVcf, recalTbi, tranches]`. Move `-mode SNP` from positional val into `ext.args`. Output rebuild: `.out.tbi.join(.out.vcf)`. |
| 20 | `gatk/mark_duplicates` (×2 aliases) | ✅ MIGRATED | `f8ba380` | 2 aliases (`GATK_MARK_DUPLICATES`, `_DELLY`). Pad with `[]` for the optional fasta/fai (CRAM-conversion-only). Emit rename: `.bam_tuple` → `.bam`. `ext.prefix` must include the `.bam` extension (nf-core uses it verbatim as output filename). Metrics file renamed (publish-only, harmless). |
| 11 | `gatk/genotype_gvcfs` | ✅ MIGRATED | `7958159` | 1 alias. Reorder combined_ch `[meta, tbi, vcf]` → `[meta, vcf, tbi, [], []]`. 5 ref value tuples (fasta/fai/dict/dbsnp/dbsnp_tbi — last 2 empty). Emit consumer: `.out` → `.out.vcf`. |
| 16-17 | `gatk/base_recalibrator` + `gatk/apply_bqsr` (×4 aliases) | ✅ MIGRATED | `a8a5143` | Disabled by default; low-risk migration. Pad `[meta, bam]` with `[]` for bai + intervals (gatk auto-indexes). 5 ref value tuples for BaseRecal, 3 bare paths for ApplyBQSR. ApplyBQSR forced to BAM via `ext.suffix='bam'` (default is CRAM). ext.prefix preserves `\${meta.id}.baserecalibrated`. |
| 22 | `samtools/merge` (×2 aliases) | ✅ MIGRATED | `7cc8ef3` | 2 aliases. Pad with `[]` for index slot + empty fasta value tuple. Move `-f` from hardcoded script into `ext.args`. ext.prefix appends `.\${file_type}` automatically (gives `.sorted_reads.bam`). Emit consumer: `.out` → `.out.bam`. |
| 21 | `bcftools/view` (×1 active alias — _DELLY) | ✅ MIGRATED | `7220e3c` | Local ran 3 commands (view+bgzip+index); nf-core does it in ONE via `-Ob` + `--write-index=csi`. Adapter joins DELLY's `.bcf` + `.csi` emits to feed nf-core's `[meta, vcf, index]` input. Output rebuild: `.out.vcf.join(.out.index)` to match the local's 3-tuple flatten chain. |
| 23 | `gatk/haplotype_caller` | ✅ MIGRATED | `aeef5f4` | 1 alias. Reorder `[meta, bai, bam]` → `[meta, bam, bai, [], []]`. 5 ref value tuples (fasta/fai/dict/dbsnp/dbsnp_tbi — last 2 empty). `-ERC GVCF` and `--bam-output ${meta.id}.haplotype_caller.bam` moved into `ext.args`; `ext.prefix = '${meta.id}.g'` preserves the `.g.vcf.gz` filename. Output rebuild: `.out.tbi.join(.out.vcf)` for the downstream gvcf_ch chain. |
| 19 | `gatk/select_variants` (×2 of 5 aliases) | ✅ MIGRATED | `91a0c72` | SNP + INDEL aliases migrated. PHYLOGENY (multi-file `-XL`), EXCLUSION_SNP, EXCLUSION_INDEL (need `--exclude-intervals` instead of nf-core's auto `--intervals`) stay local. 8 positional inputs collapsed to 1 tuple; `--select-type-to-include` + per-alias filename via `ext.args`/`ext.prefix`. Output rebuild via `.tbi.join(.vcf)`. |
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

## Container strategy after migration (preserved, not disrupted)

The port already uses the idiomatic nf-core pattern: `withName: '<pattern>' { container = '…' }`
blocks in `conf/docker.config`, `conf/singularity.config`, `conf/apptainer.config`,
`conf/podman.config`, and `conf/abc_cluster.config`. **These patterns match by
process alias name, not by module file path**, so they keep working unchanged
when a module moves from `modules/local/` to `modules/nf-core/`.

Proof in-repo: BWA_MEM is already an nf-core module and `conf/docker.config`
routes it into `magma-container-2` via `withName: 'BWA.*|…' { container = … }`.
This is exactly how the rest of Tier A will behave post-migration — emitted
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

## Gotchas

- `nf-core modules install` auto-creates `conf/containers_*.config` files and a
  scratch `magma-results/` dir. They define `withName: 'FASTQC' { container = '…' }`
  per OS/arch/runtime — if loaded they would **shadow the bundled-container
  profiles** that route tools into `magma-container-1/2`. Leave them untracked
  and do not register them as profiles. The nf-core module's in-file container
  default is enough fallback for mode 2 above.
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
