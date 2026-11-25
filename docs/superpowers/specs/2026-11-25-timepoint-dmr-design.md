# Timepoint DMR on Cirro — design

Date: 2026-11-25. Status: built locally (config repo and fork branch `timepoint-dmr`); not pushed, registered or run on Cirro.

## Intent and why

Yann's lab wants to see how tumor methylation changes between timepoints of the same patient: for a
patient with tumors S1, S4 and S6, differentially methylated regions (DMRs) of **S4 vs S1** and
**S6 vs S1**. The lrsomatic `dev` branch already calls DMRs, but only between the two haplotypes of
one sample. Timepoint comparison is missing there. The same comparison existed in the earlier NCP
pipeline (`BTC_ONT_test_pipeline`, `DMR_ONLY`, S6 vs S1) but not on Cirro.

Lab members launch pipelines from the Cirro web form, so this must run as a Cirro process.

## Goals

- A Cirro process that takes finished lrsomatic-dev datasets, and for one named **baseline** tumor
  sample compares every other tumor sample in the selected datasets against it.
- One DMR result per comparison: DMR calls restricted to CpG islands covered in both samples,
  annotated with the nearest gene.
- Cheap and quick: minutes, no re-run of lrsomatic.
- A clear, early error when the inputs cannot work.

## Non-goals

- Timepoint order or dates: the user names the baseline; all other samples are treated as later.
- Comparisons other than baseline vs each other sample (no all-pairs, no adjacent pairs).
- The NCP HTML report, promoter-restricted plots and IGV snapshots (possible follow-up).
- Haplotype DMR (already in lrsomatic `dev`), tumor-vs-normal, PacBio-specific handling, other genomes
  than GRCh38.
- Verifying that all samples are from the same patient (see Constraints).

## Anti-goals

- Do not fork the lrsomatic logic or make `cirro-dev` diverge from upstream `dev` more than by new
  files. Do not change the registered 1.0 process or the lrsomatic-dev test process's behaviour.
- Do not produce a comparison that silently contains no real data (empty region sets, mismatched
  samples, a modification mix `modkit dmr pair` rejects).

## Constraints

- Inputs are the tumor bedMethyl files that lrsomatic `dev` publishes in
  `data/<sample>/methylation/tumor/modkit_pileup/` (layout seen in run 6774e393 on the 1.0 pin, names from
  `docs/output.md` of `dev`). The whole-sample file is `<sample>.bed.gz` without `modkit_phased`, and
  `<sample>_combined.bed.gz` with it (then there is no plain `<sample>.bed.gz`, only `_hp1`, `_hp2` and
  `_combined`). `_hp1` and `_hp2` must never be used. The normal's files under `.../methylation/normal/` are
  named after the tumor and must be ignored.
- `modkit dmr pair` rejects rows with `N_other > 0`, and SUP basecalls also call 5hmC. lrsomatic `dev`
  adds 5hmC to the pileup only when `modkit_phased` is on, `skip_dmr` is off and the DMR BEDs are
  known. The source runs must therefore have been launched with "Phased methylation pileup" ticked
  (or `--modkit_args` including 5hmC).
- The registered sample metadata (`patient`, `status`) is linked to the **input** BAM datasets, not
  to lrsomatic's output datasets (checked on run 6774e393). The child process cannot read the patient,
  so it cannot verify that the samples belong to one patient.
- Nextflow 26.04.x (the code reuses `dev`'s modules, which need it). Cirro stores the Nextflow version
  per process (`executorVersion`); 26.04.0 and 26.04.6 are in use by other processes on the instance.
  Not yet run by us.
- Repository rules: nothing is pushed, registered or run without the user's approval.

## Approaches considered

- **A. Inside the lrsomatic run.** Add a `timepoint` samplesheet column and a subworkflow in the fork
  that pairs same-patient tumors after the pileup. Rejected: all timepoints must be launched
  together (about 20 h and $75 per sample), the fork diverges further from upstream, and the
  samplesheet schema needs a column upstream lacks.
- **B. A separate Cirro process that reads finished lrsomatic datasets.** Chosen by the user.
  Timepoints are usually sequenced at different times; the comparison costs minutes; the fork only gains
  new files. Cost: one more process to register and keep in step.
- Re-running the NCP `DMR_ONLY` on the cluster: rejected, it does not run on Cirro.

## The design

### Flow

1. The user selects the lrsomatic-dev result datasets (one or several) and fills in the form:
   **Baseline sample** (name, for example `S1`).
2. **Preprocess** (Python, like `community/lrsomatic/1.0/preprocess.py`, with a pure function that is
   unit-tested without Cirro) lists the selected files, keeps the whole-sample tumor pileup
   (`data/<sample>/methylation/tumor/modkit_pileup/<sample>.bed.gz` or `<sample>_combined.bed.gz`), takes the
   sample name from the path, and writes `config/comparisons.csv` with columns
   `comparison,treatment_name,treatment_bed,control_name,control_bed`, one row per non-baseline sample:
   `S4_vs_S1` (treatment S4, control S1) and `S6_vs_S1`.
3. **Pipeline** (Nextflow, new entry workflow) reads that CSV and for each row runs the existing
   `dev` modules unchanged:
   - `HTSLIB_BGZIPTABIX` to index both bedMethyl files,
   - `DMR_HAPLOTYPE_REGIONS` (the module is symmetric: it keeps CpG islands covered in both files; baseline
     goes in the "hp1" slot, the later sample in "hp2"),
   - `MODKIT_DMR` with the later timepoint as `-a` (treatment) and the baseline as `-b` (control),
   - `DMR_NEAREST_GENE`.
4. Output dataset: `data/<comparison>/dmr/<comparison>.nearestGene.bed` and `<comparison>.log`.

### Validation (preprocess, one error listing every problem)

- The baseline name matches exactly one tumor pileup; it is a required field.
- At least one other tumor pileup exists.
- Sample names are unique across the selected datasets.
- A sample with both a plain and a `_combined` file, or with more than one candidate, is an error naming the files.
- Files that are not whole-sample tumor pileups (including `_hp1`, `_hp2` and everything under `normal/`) are ignored,
  with a log line naming them.
- No pileup files found at all: say the datasets are probably not lrsomatic-dev outputs (or were run
  without the methylation step).

### Early check in the pipeline

Before running `modkit dmr pair`, read the first part of each bedMethyl and stop with a clear message
("re-run lrsomatic with Phased methylation pileup ticked") when there are rows with `N_other > 0` and no
5hmC rows. The exact test is left to the builder; the requirement is a readable failure in minutes
rather than an opaque modkit error.

### The Cirro process

- Id `process-intgenomicslab-timepoint-dmr-dev-test`; parent
  `process-intgenomicslab-lrsomatic-dev-test`; allow multiple sources; Nextflow 26.04.6; GRCh38 only.
- Form: Baseline sample (text), plus the usual extra-parameters JSON (same protection rules as 1.0).
- Output rules for `*.nearestGene.bed` and `.log`.
- The lrsomatic-dev test process gains output rules for the whole-sample tumor pileup (plain and `_combined` names) so Cirro
  exposes the files to the child process. This is a change to the unregistered `community/lrsomatic-dev/0.1`.

## Decisions

Decided by Yann: approach B; compare tumors of the same patient across timepoints; baseline vs every other
sample (S4 vs S1, S6 vs S1); the baseline is named in the form; more than two samples per run are possible;
the inputs are NEW lrsomatic-dev runs launched with "Phased methylation pileup" ticked (each full-size sample about
$75 and 20 h). Existing bedMethyl files, such as the NCP S1/S6 outputs, are not an input and a way to upload them is out of scope.

My calls (proposals, say so if wrong):

- **Code location:** new files in the fork on a branch from `cirro-dev` (an entry script `timepoint_dmr.nf` and a small
  workflow), reusing `dev`'s modules at their existing paths with no copying. Fallback if Cirro cannot use a
  non-`main.nf` entry point: a new repository with the three modules copied in.
- **Direction:** later sample = treatment, baseline = control, so a positive difference means more
  methylation at the later timepoint (to confirm against modkit's `-a` vs `-b` convention at build time).
- **Treated as the same patient** by the user's choice of datasets; the output README and logs say that this is not checked.
- **Only the lrsomatic-dev test process as parent** (not 1.0 runs: their pileup requests 5mC only, so `modkit dmr pair` would reject SUP data).
- **Baseline is the full Cirro sample name** as it appears in the pileup path (for example `P4_tumor_full`), not a suffix.
  Output labels are `<sample>_vs_<baseline>`. Names are unique across the selected datasets; a collision is an
  error, not resolved by guessing. Several patients in one selection is the user's responsibility.
- **Reuse is firm for the three DMR modules and `HTSLIB_BGZIPTABIX`** (read and confirmed symmetric: they take two
  bedMethyl+index pairs, a CpG-island BED and a fai). New workflow and subworkflow files are allowed.
- **References:** GRCh38 only, fixed (not a form field). CpG-island and gene BEDs come from `dev`'s genome defaults
  (the `dmr_*` entries in `igenomes.config`, raw GitHub URLs); fasta and fai from `dev`'s reference preparation
  on the same iGenomes source lrsomatic uses. The builder reuses that preparation rather than adding a staging task, unless it cannot.
- **Output content is whatever `dev`'s modules produce**: no extra significance filter, no minimum CpG count, no new columns.
- **Order of work:** build and test everything that needs no platform first (preprocess, stub test); then the unverified items below are
  checked with one minimal registration, each step needing the user's approval at that time. If item 2 fails, switch to the fallback.
- **The early 5hmC check** lives in the pipeline (it reads the bedMethyl), not in preprocess, which sees no file contents.
- Output names carry the comparison label (`S6_vs_S1`).

## Unverified (to check during the build, in this order)

1. That a Cirro child process sees the pileup files and their paths in its preprocess step, and that the
   parent's output rules are what exposes them.
2. That Cirro accepts a non-`main.nf` entry point for a GitHub pipeline; otherwise use the fallback above.
3. That Nextflow 26.04.6 runs on Cirro at all.
4. The modkit direction convention, and that `modkit dmr pair` works on the unphased pileup once 5hmC is included.
5. How the fasta and fai (needed by `MODKIT_DMR` and `DMR_HAPLOTYPE_REGIONS`) are staged on Cirro, as `dev` prepares
   them inside its main workflow; the builder reuses that preparation or stages them like our other references.

## Left to the builder

- File and module layout, process and parameter names, the exact 5hmC check, the CSV details, and the form wording.
- Tests: preprocess unit tests (many samples, missing baseline, duplicate names, no pileups, a pair of datasets),
  consistency tests like `tests/test_dev_process.py`, and a Nextflow stub test for the new workflow.
- Docs: a README in the process folder stating the untested items, and a how-to paragraph.

## Out of scope for the first build

Registering the process, pushing anything, running on Cirro. Each needs the user's approval.
