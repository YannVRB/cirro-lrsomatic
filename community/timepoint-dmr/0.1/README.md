# Timepoint DMR (TEST process)

Differentially methylated regions (DMRs) of each tumor against a baseline tumor of the same patient, for example
S4 vs S1 and S6 vs S1. Process id `process-intgenomicslab-timepoint-dmr-dev-test`, a child of
`process-intgenomicslab-lrsomatic-dev-test` (`community/lrsomatic-dev/0.1`).

Status: **prepared, not registered, not run.** Nothing here has been tested on Cirro. Design:
`docs/superpowers/specs/2026-11-25-timepoint-dmr-design.md`; plan: `docs/superpowers/plans/2026-11-25-timepoint-dmr.md`.

## How it works

1. You select the finished lrsomatic-dev result datasets (one or several) and type the **baseline sample**, the full
   sample name of the earliest timepoint (for example `P4_tumor_full`).
2. `preprocess.py` keeps each sample's whole-sample tumor bedMethyl
   (`data/<sample>/methylation/tumor/modkit_pileup/<sample>.bed.gz` or `<sample>_combined.bed.gz`; never `_hp1`, `_hp2`
   or anything under `normal/`) and writes `config/comparisons.csv`: one row per sample other than the baseline,
   `comparison,treatment_name,treatment_bed,control_name,control_bed`, with the later sample as treatment.
3. `timepoint_dmr.nf` (in the lrsomatic fork) runs the dev branch's DMR modules on each row: it checks the pileup for
   5hmC, indexes both files, keeps the CpG islands covered in both, runs `modkit dmr pair` (treatment `-a`, baseline
   `-b`) and adds the nearest gene.
4. Output per comparison: `data/<sample>_vs_<baseline>/dmr/<sample>_vs_<baseline>.nearestGene.bed` and a `.log`.

## What must be true of the source runs

- Launched with **Phased methylation pileup** ticked (or `--modkit_args` including 5hmC). `modkit dmr pair` rejects
  rows with other modifications, and SUP basecalls also call 5hmC. The pipeline stops with a clear message otherwise.
- Launched with the lrsomatic-dev test process. 1.0 runs only request 5mC.
- The samples belong to the same patient. **This is not checked**: the registered `patient` metadata is linked to the
  input BAM datasets, not to lrsomatic's output datasets, so this process cannot see it.

## Not verified (check in this order)

1. That a Cirro child process sees the pileup files and their paths in its preprocess step, and that the parent's output
   rules (`community/lrsomatic-dev/0.1/process-definition.json`) are what exposes them.
2. That Cirro accepts the entry point `timepoint_dmr.nf` (not `main.nf`). Fallback: a new repository with the modules copied in.
3. That Nextflow 26.04.6 runs on Cirro at all.
4. The modkit direction convention (positive difference = more methylated in the later sample), and that `modkit dmr pair`
   works on the whole-sample pileup once 5hmC is included.
5. How the fasta and fai are staged on Cirro (the workflow reads the fasta from the iGenomes S3 bucket and indexes it).

## Registering (needs the user's approval; not done)

```
python -u docs/registration/register_lrsomatic.py --project PROJECT_ID --folder community/timepoint-dmr/0.1            # dry run
python -u docs/registration/register_lrsomatic.py --project PROJECT_ID --folder community/timepoint-dmr/0.1 --create   # register and sync
```

The pinned fork commit must be pushed to `YannVRB/lrsomatic` and these files to the `main` branch of
`YannVRB/cirro-lrsomatic` first; the sync reads them from GitHub.
