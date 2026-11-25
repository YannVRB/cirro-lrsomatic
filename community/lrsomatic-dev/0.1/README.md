# LRsomatic dev (TEST process)

A separate Cirro process for the `dev` branch of IntGenomicsLab/lrsomatic, next to the working
`community/lrsomatic/1.0` process, which is left alone. Process id
`process-intgenomicslab-lrsomatic-dev-test`.

Status: **prepared, not registered, not run.** Nothing here has been tested on Cirro.

## What differs from 1.0

| | 1.0 | dev test |
|---|---|---|
| Pipeline | fork commit 42458112 (Nextflow 25.10.4) | fork branch `cirro-dev`, pinned by SHA (Nextflow 26.04.6) |
| Why 26.04 | | `dev` declares `nextflowVersion '!>=26.04.0'`; on 25.10 the LRSOMATIC workflow body is over Nextflow's 65,535-character limit |
| New analyses | | CoRAL (ecDNA), AmpliconClassifier (ecDNA/BFB call; BFBArchitect runs inside it), haplotype DMR, Padfoot, ReConPlot, ShatterSeek |
| Form | | adds "ecDNA and BFB", "Haplotype methylation" and three more skip boxes |
| `download_sigprofiler_genome` | not set (pipeline default false) | set false (`dev` defaults to true) |
| Compute config | label sizes, retryScale, SAVANA_RUN and SEVERUS overrides | only retry rule, resource cap and on-demand queue: `dev`'s `base.config` sizes by process name |
| `preprocess.py`, `process-output.json` | | identical copies (a test checks this) |

## Things to know before the first run

- **Haplotype DMR needs "Phased methylation pileup" ticked.** It is off by default, so DMR does
  nothing unless the box is ticked. The pileup also needs small-variant calling and BAMs that carry
  modification tags (MM/ML).
- **Haplotype DMR compares haplotype 1 with haplotype 2 inside one sample.** It is not a comparison
  between timepoints.
- **Downloads during the run (not staged):** the AmpliconArchitect GRCh38 data repository for
  AmpliconClassifier (md5-checked), plus ClinVar and AlphaMissense for VEP plugins. The DMR CpG-island
  and gene BEDs come from `raw.githubusercontent.com`. Turn off the first with "Skip AmpliconClassifier only".
- **Unmeasured resources:** CORAL_RECONSTRUCT, CORAL_CYCLE, CORAL_PLOT, AMPLICONCLASSIFIER, the DMR
  steps and the new plotting steps run at their label sizes. Check `trace.txt` after the first
  full-size run.
- **Optional steps ignore their own failures** (`conf/modules.config` in the fork): a CoRAL or
  AmpliconClassifier failure does not stop the run, so look at the trace and not only the run status.
- **Not verified:** that Cirro runs this on Nextflow 26.04.6 (18 other processes on the instance pin
  26.04.x; none was run by us), the new form fields in the web form, and every new step on Cirro.

## Registering (needs approval, not done)

```
python -u docs/registration/register_lrsomatic.py --project PROJECT_ID --folder community/lrsomatic-dev/0.1            # dry run
python -u docs/registration/register_lrsomatic.py --project PROJECT_ID --folder community/lrsomatic-dev/0.1 --create   # register and sync
```

The pinned fork commit must be pushed to `YannVRB/lrsomatic` and these files to the `main` branch of
`YannVRB/cirro-lrsomatic` first; the sync reads them from GitHub.
