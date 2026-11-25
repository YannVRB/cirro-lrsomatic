# Timepoint DMR Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Cirro process that compares every other tumor in selected lrsomatic-dev result datasets with a named baseline tumor (S4 vs S1, S6 vs S1) using `dev`'s DMR modules.

**Architecture:** A Python preprocess step (pure function, unit-tested) turns the selected datasets' files into `config/comparisons.csv`. A new Nextflow entry script in the lrsomatic fork reads that CSV and runs `dev`'s `HTSLIB_BGZIPTABIX`, `DMR_HAPLOTYPE_REGIONS`, `MODKIT_DMR` and `DMR_NEAREST_GENE` unchanged, preceded by one new 5hmC check. A second Cirro process definition (separate folder) wires it up. Nothing is pushed, registered or run on Cirro.

**Tech Stack:** Python 3 + pandas + pytest (config repo, env `nf-core`), Nextflow 26.04.6 DSL2 + nf-test (fork, env `nf-env` for Java 21 and Nextflow, `nf-test` from `/home/yannvrb/.conda/envs/nf-core/bin`), bash + awk.

**Spec:** `docs/superpowers/specs/2026-11-25-timepoint-dmr-design.md` (in the config repo `cirro-lrsomatic`)

## Global Constraints

- Repos: config repo `/net/bmc-lab3/data/bcc/yannvrb/GBM_TeamLab/Cirro_LRsomatic/cirro-lrsomatic` (branch `master`); fork `/net/bmc-lab3/data/bcc/yannvrb/GBM_TeamLab/Cirro_LRsomatic/lrsomatic-fork`, new branch `timepoint-dmr` created from local `cirro-dev` (tip `c816fd53fd24411a3d7e7214d224e34fe560b357`).
- Process id `process-intgenomicslab-timepoint-dmr-dev-test`; parent `process-intgenomicslab-lrsomatic-dev-test`; `allowMultipleSources` true; GRCh38 only; `executorVersion` `26.04.6`; config folder `community/timepoint-dmr/0.1`.
- Inputs are whole-sample TUMOR pileups only: `data/<sample>/methylation/tumor/modkit_pileup/<sample>.bed.gz` or `<sample>_combined.bed.gz`. Never `_hp1`/`_hp2`, never anything under `normal/`.
- Baseline = the full Cirro sample name typed in the form; every other tumor is compared with it. Treatment = the later (non-baseline) sample, control = baseline. Comparison label `<sample>_vs_<baseline>`.
- Output per comparison: `<label>/dmr/<label>.nearestGene.bed` and `<label>.log`. No extra filtering or columns.
- DMR modules and `HTSLIB_BGZIPTABIX` are reused UNCHANGED. The only edits to existing fork files are two additive lines in `nextflow.config` (Task 6) and one existing process definition in the config repo (Task 4).
- Nextflow and nf-test must not run on the login node: use `srun`/`sbatch` (default request <= 64G, <= 4h). pytest on the config repo is light and may run on the login node.
- Commits: author/committer from `git log -1` of the repo, message ends with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`. Show `git status` and `git diff` before each commit. Never `git push`. Never register in Cirro. Never touch `community/lrsomatic/1.0`.

## Review Focus

Inputs and failure modes the spec implies; each has a test in the task named.

1. Baseline `S1` when `S10` also exists: exact match only, never a prefix (Task 1).
2. The same sample name in two selected datasets: an error naming both files, no guess (Task 1).
3. `_hp1`, `_hp2` and the normal's pileup (named after the tumor) are never picked (Task 1).
4. A sample with both `<s>.bed.gz` and `<s>_combined.bed.gz`: an error (Task 1).
5. A HAC pileup that has no 5hmC rows but also no `N_other` rows must NOT trigger the check's failure (Task 5).

---

### Task 1: Preprocess logic (pure functions)

**Files:**
- Create: `community/timepoint-dmr/0.1/preprocess.py`
- Test: `tests/test_timepoint_preprocess.py`

**Interfaces:**
- Produces (in `preprocess.py`):
  - `COMPARISON_COLUMNS = ["comparison", "treatment_name", "treatment_bed", "control_name", "control_bed"]`
  - `tumor_pileups(files: pd.DataFrame) -> tuple[dict[str, list[str]], list[str]]`: from a frame with a `file` column (columns matched case-insensitively), returns `{sample: [paths of whole-sample tumor pileups]}` and the list of ignored `.bed.gz`/other paths. A path counts when it matches `(?:^|/)data/(?P<sample>[^/]+)/methylation/tumor/modkit_pileup/(?P<stem>[^/]+)\.bed\.gz$` and `stem` is `sample` or `sample + "_combined"`. S3 URIs and relative paths both work.
  - `build_comparisons(files: pd.DataFrame, baseline: str) -> pd.DataFrame`: one row per non-baseline sample, sorted by sample name, columns `COMPARISON_COLUMNS`; `comparison = f"{sample}_vs_{baseline}"`, treatment = sample, control = baseline. Raises one `ValueError("Invalid timepoint DMR input:\n- " + "\n- ".join(problems))` listing every problem.
- Problems that must be reported (one message each, naming the samples or files): baseline blank; baseline not found (list the samples that were found); no pileup files at all (say the datasets are probably not lrsomatic-dev outputs, or were run without "Phased methylation pileup"); no other tumor besides the baseline; a sample with more than one candidate file (covers the same name in two datasets, and plain plus `_combined`).

- [ ] **Step 1: Write the failing tests** in `tests/test_timepoint_preprocess.py` (load `preprocess.py` by path like `tests/test_preprocess.py` does). Helper `frame(*paths)` builds `pd.DataFrame({"file": paths})`; helper `pile(sample, suffix=".bed.gz", kind="tumor", prefix="s3://b/ds1/data")` returns `f"{prefix}/{sample}/methylation/{kind}/modkit_pileup/{sample}{suffix}"`. Tests, each with the stated assertion:
  - `test_three_samples_give_two_comparisons`: samples S1, S4, S6, baseline `S1`: `comparison == ["S4_vs_S1", "S6_vs_S1"]`, `treatment_name == ["S4","S6"]`, `control_name == ["S1","S1"]`, `control_bed` both equal S1's path.
  - `test_combined_file_is_accepted`: files `pile("S1", "_combined.bed.gz")` and `pile("S4", "_combined.bed.gz")` give one row.
  - `test_baseline_is_an_exact_match` (Review Focus 1): samples S1, S10, S4 with baseline `S1`: comparisons `["S10_vs_S1","S4_vs_S1"]`; with baseline `S` raises and the message lists `S1`, `S10`, `S4`.
  - `test_same_sample_in_two_datasets_is_an_error` (2): `pile("S4", prefix="s3://b/ds1/data")` and `pile("S4", prefix="s3://b/ds2/data")` plus S1: raises, message contains both paths.
  - `test_hp_and_normal_files_are_never_used` (3): add `pile("S4","_hp1.bed.gz")`, `pile("S4","_hp2.bed.gz")`, `pile("S4", kind="normal")` next to the real S4 and S1 files: the result uses only the tumor `S4.bed.gz`; `tumor_pileups(...)[1]` contains the hp and normal paths.
  - `test_plain_and_combined_for_one_sample_is_an_error` (4): S4 has both names: raises naming both.
  - `test_missing_baseline_blank_baseline_no_files_no_other_sample`: each raises, with the stated message fragments ("not found", "blank", "no ... pileup", "no other tumor").
  - `test_all_problems_are_listed_together`: a blank baseline with no files reports both problems in one exception.
  - `test_relative_paths_work`: paths starting `data/S1/methylation/tumor/modkit_pileup/S1.bed.gz` (no URI) give the same result.

- [ ] **Step 2: Run to verify failure**: `/home/yannvrb/.conda/envs/nf-core/bin/python -m pytest tests/test_timepoint_preprocess.py -q -p no:cacheprovider` — Expected: FAIL (module missing).
- [ ] **Step 3: Implement** `tumor_pileups` and `build_comparisons` in `preprocess.py` (module docstring in the style of the 1.0 file; `cirro` imported only inside `main()` in Task 2; reuse the 1.0 helpers' approach for stripped text and lower-cased columns, copied, not imported).
- [ ] **Step 4: Run to verify pass**: same command — Expected: all pass.
- [ ] **Step 5: Commit** (`Add timepoint DMR preprocess logic`).

---

### Task 2: Preprocess main and extra parameters

**Files:**
- Modify: `community/timepoint-dmr/0.1/preprocess.py`
- Test: `tests/test_timepoint_preprocess.py`

**Interfaces:**
- Consumes: `build_comparisons`, `COMPARISON_COLUMNS` (Task 1).
- Produces:
  - `_PROTECTED_PARAMS: frozenset[str]` = `{"comparisons", "outdir", "igenomes_base", "monochrome_logs", "custom_config_base", "custom_config_version", "validate_params", "igenomes_ignore", "genome"}`.
  - `extra_params(raw: object) -> dict[str, object]`: the same contract and messages as `community/lrsomatic/1.0/preprocess.py::extra_params` (blank gives `{}`, invalid JSON, non-object, over 10,000 characters, or a protected key raise `ValueError`), with this module's protected set.
  - `main() -> None`: `PreprocessDataset.from_running()`; reads form param `baseline_sample`; writes `build_comparisons(ds.files, baseline)` as CSV to `ds.params["comparisons"]` (an S3 URI mapped in `process-input.json`); logs the comparison table and any ignored files; applies `extra_params(ds.params.get("extra_params_json"))` with `ds.add_param(key, value, overwrite=True)`; removes `baseline_sample` and `extra_params_json` with `ds.remove_param(..., force=True)`.

- [ ] **Step 1: Write failing tests**: `test_extra_params_blank_is_empty`, `test_extra_params_rejects_protected_key` (key `comparisons`), `test_extra_params_rejects_non_object_and_bad_json`, `test_extra_params_length_limit` (10,001 characters), `test_protected_set_matches_the_form_description` is deferred to Task 3. Same style as the 1.0 `tests/test_preprocess.py` extra-params tests.
- [ ] **Step 2: Run**: Expected FAIL. **Step 3: Implement** `extra_params` and `main()` (not unit-tested: needs Cirro). **Step 4: Run**: Expected PASS.
- [ ] **Step 5: Commit** (`Add timepoint DMR preprocess entry point`).

---

### Task 3: Cirro process folder and consistency tests

**Files:**
- Create: `community/timepoint-dmr/0.1/process-definition.json`, `process-form.json`, `process-input.json`, `process-output.json`, `process-compute.config`, `README.md`
- Test: `tests/test_timepoint_process.py`

**Interfaces:**
- Consumes: `_PROTECTED_PARAMS` (Task 2).
- Produces: a registrable process folder (`register_lrsomatic.py --folder community/timepoint-dmr/0.1` dry-runs).

Contents (decisions):
- `process-definition.json`: copy the structure of `community/lrsomatic-dev/0.1/process-definition.json`; `id`, `name` "Timepoint DMR (lrsomatic dev, test)", `parentProcessIds` `["process-intgenomicslab-lrsomatic-dev-test"]`, `code.uri` `yannvrb/lrsomatic`, `code.script` `timepoint_dmr.nf`, `code.version` a placeholder of 40 zeros until Task 7, `executorVersion` `26.04.6`, `usesSampleSheet` false, `fileMappingRules` with one rule for `(?P<sampleName>[^/]+)/dmr/[^/]+\.nearestGene\.bed$` and one for `.../dmr/[^/]+\.log$`; keep the `<RESOURCES_BUCKET>` placeholders and the `documentationUrl` placeholder as in the other definitions.
- `process-form.json`: a required text field `baseline_sample` ("Baseline sample", description: the full sample name of the earliest timepoint; every other tumor is compared with it) and the `extra_params_json` textarea (description naming every protected param, as in 1.0).
- `process-input.json`: `baseline_sample` and `extra_params_json` from `$.dataset.params...`; `comparisons` = `$.dataset.s3|/config/comparisons.csv`; `outdir` = `$.dataset.dataPath`; `genome` = `GRCh38`; `monochrome_logs` = true.
- `process-compute.config`: the same content as `community/lrsomatic-dev/0.1/process-compute.config` with the header comment adapted (small job: minutes; label sizes from `dev`).
- `README.md`: status "prepared, not registered, not run"; how it works; what must be true of the source runs (launched with "Phased methylation pileup" ticked); the unverified list from the spec.

- [ ] **Step 1: Write failing tests** in `tests/test_timepoint_process.py` (same helpers and skip-if-SHA-missing style as `tests/test_dev_process.py`):
  - `test_definition`: id and parent as above; `code["script"] == "timepoint_dmr.nf"`; `executorVersion` >= the pinned manifest's `nextflowVersion` when the SHA exists in the fork (skip while it is the zero placeholder); `usesSampleSheet is False`; `allowMultipleSources is True`.
  - `test_every_dataset_param_path_resolves_to_a_form_property` and `test_baseline_is_required_text` (`"baseline_sample" in form["required"]`, type string).
  - `test_form_description_names_every_protected_param` (whole-word, as in the 1.0 test).
  - `test_input_mapping`: exactly the keys listed above, `genome == "GRCh38"`.
  - `test_output_rules_match_examples`: each rule's example name matches its own pattern (convert `(?P<` as in `tests/test_dev_process.py`), and `s1_vs_s0/dmr/s1_vs_s0.log` is not matched by the nearestGene rule.
  - `test_compute_config_matches_dev_test_process`: after removing comments, equal to the dev test config's code.
  - `test_preprocess_loads_and_exposes_the_contract`: `COMPARISON_COLUMNS` is the five columns the pipeline script reads (Task 6 asserts the same list from the other side).
- [ ] **Step 2: Run** `pytest tests -q -p no:cacheprovider`: new tests FAIL, old ones pass. **Step 3: Create the six files.** **Step 4: Run all**: Expected all pass (104 + new). **Step 5: Dry run** `/home/yannvrb/.conda/envs/cirro-sdk/bin/python docs/registration/register_lrsomatic.py --project d29bb6ed-297c-47f6-acd2-64cd6ed46c9b --folder community/timepoint-dmr/0.1` — Expected: prints the request and "DRY RUN"; nothing sent. **Step 6: Commit** (`Add the timepoint DMR process definition (not registered)`).

---

### Task 4: Expose the tumor pileups from the lrsomatic-dev test process

**Files:**
- Modify: `community/lrsomatic-dev/0.1/process-definition.json` (append one `fileMappingRules` entry)
- Modify: `tests/test_dev_process.py` (extend `test_file_mapping_rules_cover_the_new_outputs`)

**Interfaces:** the rule `"Tumor methylation pileup (bedMethyl)"`, `min` 0, one pattern `(?P<sampleName>[^/]+)/methylation/tumor/modkit_pileup/[^/]+\.bed\.gz$`, example `sampleName/methylation/tumor/modkit_pileup/sampleName.bed.gz`.

- [ ] **Step 1: Failing test**: the rule exists; its example matches; `s1/methylation/normal/modkit_pileup/s1.bed.gz` does not match; `s1/methylation/tumor/modkit_pileup/s1_combined.bed.gz` and `..._hp1.bed.gz` do match the Cirro rule (the stricter selection happens in preprocess).
- [ ] **Step 2: Run**: FAIL. **Step 3: Add the rule.** **Step 4: Run all tests**: PASS. **Step 5: Commit** (`Expose the tumor pileup from the lrsomatic-dev test process`).

---

### Task 5: Fork — toolchain and the 5hmC check

**Files (fork, new branch `timepoint-dmr` from `cirro-dev`):**
- Create: `modules/local/dmr/check_modifications/main.nf`, `modules/local/dmr/check_modifications/templates/check_modifications.sh`, `modules/local/dmr/check_modifications/tests/check_modifications.test.sh`

**Interfaces:**
- Produces: process `DMR_CHECK_MODIFICATIONS`: input `tuple val(meta), path(bedmethyl)`; output `tuple val(meta), path(bedmethyl), emit: checked` (the same file passed through) and `path "versions.yml", emit: versions`; label `process_single`; container: the one `DMR_HAPLOTYPE_REGIONS` uses (it has awk, zcat).
- `check_modifications.sh BEDGZ` exit codes: 0 when the pileup is usable, 1 with the message `"<file>: rows with other modifications (N_other > 0) and no 5hmC rows; re-run lrsomatic with 'Phased methylation pileup' ticked"` on stderr.

Decision for the script: read the first 200,000 lines; fail only when some row has bedMethyl column 14 (`Nother_mod`) > 0 AND no row has column 4 equal to `h`. Verify the column numbers against modkit's bedMethyl format (column 10 is `Nvalid_cov`, which `DMR_HAPLOTYPE_REGIONS` already uses) before relying on 14.

- [ ] **Step 1: Toolchain check** (once, on a compute node): `srun -p <partition used by earlier sbatch scripts> -c 2 --mem 8G -t 00:20:00 bash -lc 'export NXF_VER=26.04.6; PATH=/home/yannvrb/.conda/envs/nf-env/bin:/home/yannvrb/.conda/envs/nf-core/bin:$PATH; nextflow -v; nf-test version'` — Expected: both print versions. If the 26.04.6 download is blocked on compute nodes, download it on the login node first (`NXF_VER=26.04.6 nextflow -v`, a download, not a run) and re-run. Record the working command in the README of the fork branch commit message.
- [ ] **Step 2: Write the failing shell test** `check_modifications.test.sh`: builds three gzipped fixtures with `printf ... | gzip` and asserts exit codes: (a) rows with code `h` present and `N_other` > 0 → 0; (b) only code `m`, `N_other` > 0 → 1 and stderr contains "5hmC"; (c) only code `m`, all `N_other` 0 (HAC, Review Focus 5) → 0. Run `bash modules/local/dmr/check_modifications/tests/check_modifications.test.sh` — Expected FAIL (script missing).
- [ ] **Step 3: Implement** the script and `main.nf` (Nextflow `template 'check_modifications.sh'`; `stub:` runs `touch` and emits the input). **Step 4: Run the shell test** — Expected: three PASS lines.
- [ ] **Step 5: Commit in the fork** (`Add DMR_CHECK_MODIFICATIONS`).

---

### Task 6: Fork — the workflow, config and entry script

**Files (fork):**
- Create: `subworkflows/local/timepoint_dmr/main.nf`, `conf/timepoint_dmr.config`, `timepoint_dmr.nf`, `subworkflows/local/timepoint_dmr/tests/main.nf.test`
- Modify: `nextflow.config` (exactly two added lines: `comparisons = null` inside the params block next to the DMR options, and `includeConfig 'conf/timepoint_dmr.config'` after the igenomes include). This is the one deviation from "new files only" in the spec: the `ext.args` and `publishDir` settings `dev` keeps in `conf/modules.config` select `.*:DMR:...` and do not match a differently named workflow. Record this in the commit message.

**Interfaces:**
- Consumes: `DMR_CHECK_MODIFICATIONS` (Task 5); `HTSLIB_BGZIPTABIX(tuple val(meta), path(infile), path(tbi), path(regions); val action; val make_index; val out_ext)`; `DMR_HAPLOTYPE_REGIONS(tuple val(meta), path(a), path(a_tbi), path(b), path(b_tbi); path cpg_islands_bed; tuple val(meta2), path(fai))`; `MODKIT_DMR(tuple val(meta), path(a), path(a_tbi); tuple val(meta2), path(b), path(b_tbi); tuple val(meta3), path(regions); tuple val(meta4), path(fasta), path(fai))`; `DMR_NEAREST_GENE(tuple val(meta), path(bed); path gencode_gene_bed)`; `SAMTOOLS_FAIDX(tuple val(meta), path(fasta), path(fai); val get_sizes)`.
- Produces: `workflow TIMEPOINT_DMR { take: comparisons /* [meta, treatment_bed, control_bed], meta = [id: comparison] */, fasta /* [[:], fasta] */, fai /* [[:], fai] */, cpg_islands_bed /* path */, gencode_gene_bed /* path */; emit: nearest_gene_bed /* [meta, bed] */, log /* [meta, log] */, versions }`.

Decisions:
- Per comparison: run `DMR_CHECK_MODIFICATIONS` on both files (`meta + [role: 'treatment'|'control']`), index both with `HTSLIB_BGZIPTABIX` ('compress', true, 'bed') following the flatMap/join pattern of `subworkflows/local/dmr/main.nf`, run `DMR_HAPLOTYPE_REGIONS` with the CONTROL (baseline) in the first slot and the treatment in the second, drop empty region files with a `log.warn` naming the comparison (as `dev` does), then `MODKIT_DMR` with the treatment as `-a` and the control as `-b`, then `DMR_NEAREST_GENE`.
- `conf/timepoint_dmr.config`: `withName: '.*:TIMEPOINT_DMR:MODKIT_DMR' { ext.args = '--base C' }`; publish `<outdir>/${meta.id}/dmr` for `MODKIT_DMR` (log files only), `DMR_NEAREST_GENE` (not `versions.yml`); everything else `publishDir enabled: false`.
- `timepoint_dmr.nf`: errors with a clear message if `params.comparisons` is null; reads the CSV with `splitCsv(header: true)`; builds `[[id: row.comparison], file(row.treatment_bed), file(row.control_bed)]`; fasta `file(params.genomes[params.genome].fasta)` with `SAMTOOLS_FAIDX([[:], fasta, []], false)` for the fai; BEDs `file(params.genomes[params.genome].dmr_cpg_islands_bed)` and `...dmr_gencode_gene_bed`; entry workflow calls `TIMEPOINT_DMR`.
- The CSV columns are exactly `COMPARISON_COLUMNS` from Task 1 (`comparison,treatment_name,treatment_bed,control_name,control_bed`).

- [ ] **Step 1: Write the failing nf-test** `subworkflows/local/timepoint_dmr/tests/main.nf.test`: `options "-stub"`, inputs two comparisons (`S4_vs_S1`, `S6_vs_S1`) from tiny fixture files, `fasta`/`fai` placeholders; assertions: `workflow.success`; `workflow.out.nearest_gene_bed` has 2 entries whose meta ids are `S4_vs_S1` and `S6_vs_S1`; each `*.nearestGene.bed` file name starts with its meta id. A second test with one comparison whose regions file is empty (use the non-stub region module only if stubbing cannot produce an empty file; otherwise mark the empty-regions path as covered by `dev`'s own `DMR` subworkflow test and say so in the test file) asserts the run succeeds and emits no result for it.
- [ ] **Step 2: Run** (on a compute node, env from Task 5): `nf-test test subworkflows/local/timepoint_dmr/tests/main.nf.test` — Expected FAIL (workflow missing).
- [ ] **Step 3: Implement** the subworkflow, config, entry script and the two `nextflow.config` lines.
- [ ] **Step 4: Run** the nf-test — Expected PASS. Also run the entry script end to end in stub mode on a compute node: `nextflow run . -main-script timepoint_dmr.nf -stub-run --comparisons <fixture csv> --outdir <scratch>` — Expected: completes and `<scratch>/S4_vs_S1/dmr/S4_vs_S1.nearestGene.bed` exists. If `params.comparisons` is rejected as unknown, keep the declared `comparisons = null` line (that is what it is for).
- [ ] **Step 5: Confirm `dev` still parses**: `nextflow config .` on the branch — Expected: prints the config without error.
- [ ] **Step 6: Commit in the fork** (`Add the timepoint DMR workflow`).

---

### Task 7: Pin, verify, document

**Files:** `community/timepoint-dmr/0.1/process-definition.json` (config repo), the memory note `cirro-lrsomatic-wrapper.md`, `docs/superpowers/specs/2026-11-25-timepoint-dmr-design.md` (status line only)

- [ ] **Step 1:** Replace the zero SHA in `process-definition.json` with the tip of fork branch `timepoint-dmr` (`git -C <fork> rev-parse timepoint-dmr`).
- [ ] **Step 2:** Run the whole config-repo suite: `pytest tests -q -p no:cacheprovider` — Expected: all pass, none skipped for the SHA. Re-run the Task 3 dry run — Expected: request shows the real SHA and `executorVersion` `26.04.6`.
- [ ] **Step 3:** Add to the README of the process folder the exact push and registration commands (fork: `git push yannvrb timepoint-dmr`; config repo: `git push origin master:main`; registration: the dry-run command with `--create`), marked as needing the user's approval, and the unverified list from the spec in order.
- [ ] **Step 4:** Update the spec's status line to "built locally, not pushed, registered or run". Update the memory note and index with the branch name, tip SHA, commits, and the open items. Commit the config repo (`Pin the timepoint DMR process to the fork commit`).
- [ ] **Step 5:** Report to the user: what was built, test results, and the three approval gates (two pushes, registration, first chr4 test).

## Self-review (done when writing)

- Spec coverage: baseline selection, whole-sample file names, validation, 5hmC check, module reuse, outputs, parent rules, reference sources, form, order of work, README with unverified items — each maps to Tasks 1-7. Early 5hmC check lives in the pipeline (Task 5). Fallback repo for a non-`main.nf` entry point is NOT built; it is only triggered if the platform check fails later.
- Types: `COMPARISON_COLUMNS` is defined in Task 1 and consumed by Tasks 3 and 6; module signatures in Task 6 were read from `cirro-dev`.
- Proportion: the plan is longer than the spec because it carries file paths, signatures and test assertions; it contains no function bodies.
