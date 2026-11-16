# Somatic Variant Calling for Long Reads (LRsomatic)

DRAFT. This folder is a Cirro wrapper for the long-read somatic variant calling pipeline
[IntGenomicsLab/lrsomatic](https://github.com/IntGenomicsLab/lrsomatic) (tumor/normal, ONT or PacBio, AWS Batch).
It has never been registered in Cirro or run on AWS. Read "Not verified" before relying on it.

## What this is

A set of Cirro configuration files plus a `preprocess.py` that turns a Cirro dataset (BAM files and per-sample
metadata) into the pipeline samplesheet `sample,bam_tumor,bam_normal,platform,sex,fiber`, one row per tumor sample,
written to `config/pipeline_samplesheet.csv` of the output dataset. The pipeline code is not in this folder: Cirro runs
a pinned fork of it.

| File | Purpose |
|---|---|
| `process-definition.json` | The process record: pinned fork, Nextflow version, parents, output file rules |
| `process-input.json` | Maps form values and Cirro paths to pipeline parameters |
| `process-form.json` | The launch form |
| `process-compute.config` | Process/executor settings for AWS Batch (see "Compute") |
| `process-output.json` | Web-optimized outputs: none (`{"commands": []}`) |
| `preprocess.py` | Builds the samplesheet and applies `extra_params_json` |

### Pinned fork

`code.version` is the full commit `42458112ee927a37c1879a3f23e0fb82ab7b8fdb` of the fork's `cirro` branch. The fork is
upstream IntGenomicsLab/lrsomatic commit `7644f58f3dcbd577caef6e90f40fe8bd34a88b02` (in upstream `dev` history, 158
commits ahead of the 1.1.0 release; it is not 1.1.0) plus three commits:

1. `Map HAC v4.3.0 basecall model to Clair3/ClairS models`: adds `dna_r10.4.1_e8.2_400bps_hac@v4.3.0` to the Clair3 and
   ClairS model maps in `workflows/lrsomatic.nf` (upstream maps only SUP models).
2. `Fix METAEXTRACT awk SIGPIPE`: fixes an awk SIGPIPE failure in `modules/local/metaextract/main.nf`.
3. `Use SigProfiler container 1.3.6-chm13-7894689`: switches the three SigProfiler modules to that image tag, which
   accepts the GRCh38 files upstream SigProfilerMatrixGenerator serves today (the older `28a9ce8` image does not).

**Before registering:** `code.uri` is `yannvrb/lrsomatic` and the repository is declared `GITHUBPUBLIC`, so the fork
must be a public GitHub repository and that commit must be pushed there. The other
placeholders (`<RESOURCES_BUCKET>`, `<PROCESS_DIRECTORY>`, `<DOCS_SITE>`) are the ones the Cirro-pipelines catalogue
uses and are filled in at registration. Nextflow `25.10.4` or newer is required (`code.executorVersion`).

### Process definition notes

- Standard (AWS Batch) execution only; `executionModes` is omitted.
- `allowMultipleSources` is true because a patient's tumors and normal are probably separate Cirro datasets.
- `parentProcessIds` is `ingest_aligned_bam` and `ingest_unaligned_bam`. The lab's BAMs are unaligned
  (`*_merged_unaligned.bam`), so the unaligned-BAM intake type is listed as a parent. That fit was read from the
  Cirro-pipelines catalogue (its file patterns are `<sample>.unaligned.bam` and `<sample>.bam`) and was not run on
  Cirro. `custom_dataset` is deliberately not a parent: a custom dataset carries no sample association, so the pairing
  by sample name cannot work on it.
- `fileMappingRules` describe three optional (`min` 0) outputs per tumor sample, modelled on the layout of two finished
  runs (`<sample>/variants/phased/somatic_smallvariants.vcf.gz`, `<sample>/vep/somatic/<sample>_SOMATIC_VEP.vcf.gz`,
  each with its `.tbi`, and `<sample>/report/<sample>_report.html`). The VEP files exist only when VEP is not skipped.
  Everything else the pipeline writes (ASCAT, Wakhan, Savana, Severus, QC, methylation, SV and germline VEP files) is
  still in the output dataset but has no rule.

## Metadata contract

The samplesheet is built from the Cirro per-sample metadata. Column names are matched case-insensitively; values are
matched case-insensitively after trimming. Any other value is an error naming the sample, never a default. The accepted
spellings are one dict (`SPELLINGS`) at the top of `preprocess.py`; if your real Cirro metadata uses other column names
or spellings, only that dict (and the column list next to it) needs to change.

| Column | Meaning | Accepted values |
|---|---|---|
| `patient` | Groups a patient's samples; any non-blank text | free text |
| `status` | Tumor or normal sample | `tumor`, `tumour` (tumor); `normal`, `pbmc` (normal) |
| `sex` | Patient sex | `male`, `m` (male); `female`, `f` (female) |
| `platform` | Sequencing platform | `ont` (ont); `pb`, `pacbio` (pb) |
| `fiber` | Fiber-seq sample | `y`, `yes` (y); `n`, `no` (n) |

Rules, each enforced with an error that names the sample or patient:

- Only dataset files whose name ends in `.bam` (case-sensitive) are treated as BAMs. Any other file, such as a
  `<sample>.bam.bai` index, is ignored and listed in a warning in the preprocess log. An input with no `.bam` file at
  all is an error.
- One BAM per sample. A sample with more than one BAM is an error (no merge logic in v1).
- Every BAM has a metadata row and every metadata row has a BAM (a sample whose only file is not a `.bam` has no BAM).
  The sample name must not contain whitespace.
- Each patient needs exactly one normal. A patient with none, or with two or more, is an error. The normal is shared by
  all of that patient's tumors and never across patients. Tumor-only runs are not supported in v1.
- A patient's tumor and normal must have the same `platform` and the same `sex`. A mismatch is an error naming the
  patient, both sample names and both values (a differing `sex` is the one free sign of a sample filed under the wrong
  patient). `fiber` may differ, because the pipeline has `skip_normalfiber`; the samplesheet row takes the tumor
  row's `fiber` value, and the normal's `fiber` is not used.
- A patient with only normals produces no rows and a warning in the preprocess log (one per such patient); if the
  whole run has no tumor, it is an error.
- `preprocess.py` logs the full resulting samplesheet, so the pairing can be checked in the Cirro log.

## Form parameters

| Form field | Pipeline parameter | Default |
|---|---|---|
| Reference Genome | `genome` | `GRCh38` (the only option in v1; CHM13 is not offered until tested) |
| Skip ASCAT / Savana / Wakhan / VEP / VEP plugins | `skip_ascat`, `skip_savana`, `skip_wakhan`, `skip_vep`, `skip_vep_plugins` | unticked |
| Skip Fiber-seq / Fiber-seq on normals / m6A / modification calling / QC | `skip_fiber`, `skip_normalfiber`, `skip_m6a`, `skip_modcall`, `skip_qc` | unticked |
| Skip mutational signatures | `skip_signatures` | unticked (signatures run) |
| VEP cache version | `vep_cache_version` | `113` |
| Extra Parameters (JSON) | `extra_params_json` | empty |

Not shown in the form and set by `process-input.json`: `outdir` (the dataset data path), `input` (the samplesheet),
`sigprofiler_genome_dir`, `vep_revel`, `vep_revel_tbi` (all under the references bucket) and
`monochrome_logs`. `igenomes_base` and `vep_cache` are not set: the fork's defaults point at the public buckets
`s3://ngi-igenomes/igenomes/` and `s3://annotation-cache/vep_cache/`. Every other pipeline parameter keeps its
pipeline default.

### `extra_params_json`

A JSON object, applied in `preprocess.py` after the form values, so its keys override them.

- Rejected (error, not a warning): the protected keys `input`, `outdir`, `igenomes_base`, `vep_cache`,
  `sigprofiler_genome_dir`, `vep_revel`, `vep_revel_tbi`, `monochrome_logs`, `custom_config_base`,
  `custom_config_version`, `validate_params`, `igenomes_ignore`, `genome`, `download_vep_cache` and
  `download_sigprofiler_genome`. The first group is what `process-input.json` sets, plus `igenomes_base` and
  `vep_cache`, which stay at the fork's public defaults; `custom_config_base` and
  `custom_config_version` would pull a remote config; `validate_params` and `igenomes_ignore` switch off checks or the
  references tree; `genome` is set only by the form's Reference Genome field (GRCh38 in v1); the two `download_*` keys
  would fetch references that are staged instead. There is no allow-list: any other key is passed on.
- Rejected: more than 10000 characters (the form's `maxLength` is client-side only, so it is checked again in
  `preprocess.py`).
- Rejected: text that is not valid JSON, or valid JSON that is not an object.
- Blank means no extra parameters.

Because the reference paths are protected, a user cannot work around a wrong reference path below
through this box; the paths have to be fixed in `process-input.json`.

## Prerequisites (references the launch needs)

The default launch reads two references from public S3 buckets (fork defaults, listed anonymously on 2026-11-16; not
yet read from a Cirro Batch job) and two from the Cirro references bucket (`$.referencesBucketUri`). Nothing has been
staged in the references bucket by this work. A default launch fails until the REVEL and SigProfiler items below are in
place at the mapped paths.

| What | Param | Expected folder layout | Status |
|---|---|---|---|
| iGenomes GRCh38 tree | `igenomes_base`, not mapped: fork default `s3://ngi-igenomes/igenomes/` | `Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta` and `Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/{GC,RT}_G1000_hg38.zip`, plus the `Homo_sapiens/ClairSTO/GRCh38/Annotation/` files named in the fork's `conf/igenomes.config` (read from that file; which of them a paired run uses was not checked) | public bucket; the FASTA, `GC_G1000_hg38.zip`, `RT_G1000_hg38.zip` and four ClairSTO annotation files were listed anonymously on 2026-11-16 |
| VEP cache 113 | `vep_cache`, not mapped: fork default `s3://annotation-cache/vep_cache/`; `vep_cache_version` = 113 (form) | `113_GRCh38/homo_sapiens/113_GRCh38/` in that bucket (`subworkflows/local/prepare_annotation.nf` adds the `<version>_<genome>/` prefix only when `vep_cache` equals exactly this URL, so do not map it to another bucket without the unprefixed `homo_sapiens/113_GRCh38/` layout; `homo_sapiens_merged` or `homo_sapiens_refseq` if `--merged` or `--refseq` is in `vep_args`) | public bucket; that folder was listed anonymously on 2026-11-16 and matches what the fork code looks for; the code was read, not run |
| REVEL | `vep_revel` = `$.referencesBucketUri\|/vep_revel/revel_grch38.tsv.gz` and `vep_revel_tbi` = `$.referencesBucketUri\|/vep_revel/revel_grch38.tsv.gz.tbi` | two files in the `vep_revel/` folder: the prepared, bgzipped GRCh38 score file and its tabix index, named as mapped. The lab's P4 run already produced exactly this pair (made by the pipeline's prepare step from the REVEL release zip): `/net/bmc-lab3/data/bcc/yannvrb/GBM_TeamLab/P4.LRsomatic/vep_plugins/revel_grch38.tsv.gz` (about 606 MB) and `.../vep_plugins/revel_grch38.tsv.gz.tbi`. They exist there; they were not integrity-checked | upload requested from the references-bucket admin on 2026-11-16 at these paths; not confirmed staged |
| SigProfiler volume | `sigprofiler_genome_dir` = `$.referencesBucketUri\|/sigprofiler/volume/` | the volume directory itself, containing `tsb/GRCh38/` with 25 files | upload requested from the references-bucket admin on 2026-11-16 at this path; not confirmed staged |

Reference-bucket paths copied exactly from `process-input.json`. They are the paths requested from the bucket admin; if
the files end up elsewhere, change the matching param in `process-input.json`:

- `$.referencesBucketUri|/sigprofiler/volume/`
- `$.referencesBucketUri|/vep_revel/revel_grch38.tsv.gz`
- `$.referencesBucketUri|/vep_revel/revel_grch38.tsv.gz.tbi`

SigProfiler volume: a verified copy is at
`/net/bmc-lab3/data/bcc/yannvrb/GBM_TeamLab/Cirro_LRsomatic/reference/sigprofiler_volume` (about 2.9 GB;
`tsb/GRCh38/`, 25 files). The 25 files were checked 25/25 against the checksums in upstream SigProfilerMatrixGenerator's
current GRCh38 registry, and SigProfilerMatrixGenerator then ran on them with exit 0 (SBS96 total 15528 of 15528 input
SNVs, with the Singularity image). `--sigprofiler_genome_dir` skips only the pipeline's own assertion; the tool
re-checks against its registry, so the volume works only with the corrected container (the third fork commit). The
location in the bucket is the one requested, not a confirmed one.

REVEL: the fork sends any `vep_revel` value containing `://` (so an `s3://` URI too) to `wget` (`subworkflows/local/prepare_vep_plugins.nf`), which cannot read S3, so the REVEL release zip cannot be used from the references bucket. The prepared `.tsv.gz` plus its `.tbi` avoids that: the fork treats a `.gz` as already prepared and requires the index, and Nextflow stages S3 paths natively. Not run on Batch. REVEL is free for non-commercial use only.

## Compute

| Setting | Value | Why |
|---|---|---|
| `errorStrategy` / `maxRetries` | pipeline predicate plus unknown exit status (retry on exit 130-145, 104, 175-177 or `Integer.MAX_VALUE`, else `finish`) / 2 | The pipeline's predicate from the fork's `conf/base.config` is kept exactly and one condition is added, `task.exitStatus == Integer.MAX_VALUE` (Nextflow's value when no exit code was recorded). A spot reclaim on AWS Batch MAY surface that way instead of as a container exit code; that is unverified (nothing ran on Batch). Real non-zero exit codes outside the list are still not rerun. `maxRetries` is raised from the pipeline's 1 to 2 (three attempts), which also bounds retries of an unknown exit status. `aws.batch.maxSpotAttempts` is not set because its behaviour on Cirro cannot be verified. Overridden memory scales with `task.attempt`: `SAVANA_RUN` 220 / 440 / 488 GB (third attempt capped by `resourceLimits`), `SEVERUS` 96 / 192 / 288 GB. |
| `resourceLimits` | `cpus: 64, memory: 488.GB, time: 240.h` | **Assumed.** These are the sarek values; the largest AWS Batch instance of the user's compute environment is not known. Replace them with the real limits. |
| `SAVANA_RUN` memory | `220.GB * task.attempt` | Peak RSS 175.4 GB observed vs 72 GB from `process_high`. |
| `SEVERUS` memory | `96.GB * task.attempt` | Peak RSS 75.4 GB observed vs 36 GB from `process_medium`. |
| `MINIMAP2_ALIGN`, `LONGPHASE_HAPLOTAG` queue | `env('PW_ONDEMAND_JOB_QUEUE')` | Tasks above 4 h in the traces; run on-demand so a spot interruption does not discard hours of work. The variable is Cirro-specific, copied from the sarek config, and **not verified for this pipeline**. |

There is no `params {}` block (the pipeline requires parameters on the command line, never via
`-c`) and no Slurm or Singularity settings (Cirro supplies the AWS Batch ones).

The sizes come from the Slurm execution traces of two real runs (five tumor/normal pairs, about 200 GB input BAMs per
sample), a small sample; the headroom factor is 1.25 on the observed maximum. Full evidence and the list of processes
with no evidence are in `docs/compute-sizing.md` of the repository this folder lives in.

## Not verified

None of the following was ever exercised:

- Cirro itself: registering this process, the form, the `preprocess.py` hook and the output rules were not tried.
- AWS Batch: no run on Batch, so the assumed `resourceLimits` and `PW_ONDEMAND_JOB_QUEUE` (see "Compute") are untested.
- Spot-reclaim behaviour on Cirro/AWS Batch: whether a reclaimed task reaches Nextflow as an unknown exit status
  (`Integer.MAX_VALUE`, which the retry predicate now covers), as one of the listed exit codes, or in another way, was
  not observed. `aws.batch.maxSpotAttempts` is not set.
- The Docker image pull of `ghcr.io/ljwharbers/sigprofiler:1.3.6-chm13-7894689`: only the Singularity image was run.
- PacBio: there is no PacBio data; `platform=pb` and mixed ONT/PacBio runs are accepted but untested.
- `--plot` and SigProfilerAssignment in the signatures step.
- The real shape of `ds.files` for BAM inputs (the `sample` and `file` columns `preprocess.py` expects), and the real
  names of the Cirro metadata columns.
- `main()` of `preprocess.py`: the pairing logic is unit-tested, `main()` was never run against Cirro (its logging
  was only exercised with a stand-in for the Cirro dataset object).
- Whether the pipeline is offered on, and works with, datasets uploaded as "Unaligned BAM" (`ingest_unaligned_bam`):
  read from the catalogue, not run on Cirro (see "Process definition notes").
- The `fileMappingRules` regexes were only checked against file names from the finished Slurm runs, and only with
  Python's `re`, not Cirro's regex engine.
- The prepared-REVEL route (the `.tsv.gz` and `.tbi` staged in the references bucket and read by Nextflow over S3) was
  not run on Batch, and the integrity of the P4 `revel_grch38.tsv.gz` / `.tbi` pair was not checked.
- The Cirro references bucket was never listed and nothing was staged there. The two public buckets (iGenomes, VEP
  cache) were listed anonymously from the lab cluster on 2026-11-16, not from a Cirro Batch job.

## Verification

Run on 2026-11-13 against fork commit `42458112ee927a37c1879a3f23e0fb82ab7b8fdb`. Full suite after the final review
fixes: `python -I -m pytest tests -q -W error` gives `93 passed`, none skipped. Passing here
means the checks below passed; it does not mean the process works in Cirro (see "Not verified").

| # | Check | Result |
|---|---|---|
| 1 | Cirro `validate.py` on a throwaway copy of Cirro-pipelines with this folder added at `community/lrsomatic/1.0/`: `python -I validate.py` | PASS. Exit 0, no ERROR lines, "All processes pass validation". Proves the form, definition and input JSON match Cirro's schemas. Does not prove Cirro accepts the process at registration or that the paths and IDs exist. |
| 2 | `python -I docs/verification/verify_real_samplesheets.py` rebuilds the real P3 and P4 samplesheets from a Cirro-style `files` table (tumor and normal in different datasets) and a metadata table, both patients in one `build_pipeline_samplesheet` call | PASS, 5/5 (P3: PT3.S1.L2, PT3.S6.L1; P4: PT4.S1.L1, PT4.S4.L2, PT4.S6.L6): same samples, same platform/sex/fiber, each `bam_normal` is the URI of that tumor's real normal BAM. Proves the pairing logic on the real P3/P4 layout. Does not prove the real shape of Cirro's `ds.files` or the real metadata column names: the table was built by the script, not by Cirro. |
| 3 | Failure messages: `pytest tests/test_preprocess.py -v -k "patient_without_normal or patient_with_two_normals or bad_sex or bad_platform"` | PASS, 4 passed (`test_patient_without_normal_raises`, `test_patient_with_two_normals_raises`, `test_bad_sex_names_sample`, `test_bad_platform_raises`). The messages are shown below. Proves each case is rejected with a message naming the patient or sample; does not prove how Cirro displays an exception raised in the hook. |
| 4 | Mapped-param consistency: `pytest tests/test_config_consistency.py -v` | PASS, 10 passed. Every param mapped in `process-input.json` exists in `nextflow_schema.json` of the fork working copy (`lrsomatic-fork/nextflow_schema.json`, branch `cirro`, clean tree at the pinned commit), and the form field types agree with it. Does not prove the pipeline accepts the values at runtime. |
| 6 | Behaviour added by the final review fixes: `pytest tests -v -k "bai or non_bam or no_bam or ignored or case_sensitive or sex_mismatch or fiber_difference or normal_only or protected or unknown_exit or parents"` | PASS, 31 passed. `.bam.bai` rows are ignored and listed, an input with no `.bam` raises, a tumor/normal `sex` mismatch raises (a `fiber` difference does not, the tumor's value is used), `normal_only_patients` finds patients without a tumor, the seven added protected keys are rejected and named in the form text, the retry predicate contains `Integer.MAX_VALUE`, the parents are exactly `ingest_aligned_bam` and `ingest_unaligned_bam`. Proves the code and files behave as stated on small in-memory tables and on the JSON; does not prove how Cirro supplies `ds.files` (whether a `.bai` row ever appears), that a spot reclaim surfaces as an unknown exit status, or that the unaligned-BAM intake offers this pipeline. |

Re-run after the final review fixes, 2026-11-13:

- Full suite: `python -I -m pytest tests -q -W error`: 93 passed, none skipped.
- Check 2 (`python -B -I docs/verification/verify_real_samplesheets.py`): `TOTAL 5/5 PASS, 5 rows built, 0 failure(s)`.
  The script needed no change; its fixtures hold only `.bam` rows and the P3/P4 tumor and normal have the same sex.
- Check 1 (`validate.py` on a fresh throwaway copy of Cirro-pipelines with this folder at `community/lrsomatic/1.0/`):
  exit 0, no ERROR lines, "All processes pass validation" (the form, definition and input JSON schemas).
- Nextflow parse of `process-compute.config` in the fork through `srun` (`nextflow -c process-compute.config config
  -profile singularity`): `CONFIG_OK`. The resolved `process.errorStrategy` is
  `{ (task.exitStatus in ((130..145) + 104 + (175..177)) || task.exitStatus == Integer.MAX_VALUE) ? 'retry' : 'finish' }`
  with `maxRetries = 2`. The cherry-pick part of check 5 was not repeated (the fork did not change).
- Still not proven by any of these: the "Not verified" list above is unchanged except for the additions made with these
  fixes (spot-reclaim behaviour, the unaligned-BAM intake).
| 5 | Fork: in a fresh `git clone --shared`, detach at `7644f58f3dcbd577caef6e90f40fe8bd34a88b02` and cherry-pick the three `cirro` commits in order; then `nextflow -c process-compute.config config -profile singularity` in the fork through `srun` (Nextflow 26.04.6) | PASS. All three applied cleanly and `git diff cirro` was empty (identical tree `491f3873f899d19bef8eda607e2ac8ba6afce486`; the commit SHAs differ because cherry-picking sets new committer dates, so the pin is the original SHA). The parse printed `CONFIG_OK`. Proves the three changes are a clean patch on that upstream commit and that the fork plus the compute config load in Nextflow. Does not prove a run succeeds (no process was executed). |

What a user sees for each failure case (the text of the `ValueError`, from a small in-memory P3 example):

```
patient with no normal:
Invalid samplesheet input:
- patient 'P3' has 0 normal sample(s) [], expected exactly 1

patient with two normals:
Invalid samplesheet input:
- patient 'P3' has 2 normal sample(s) ['P3_PBMC', 'P3_PBMC2'], expected exactly 1

bad sex:
Invalid samplesheet input:
- sample 'PT3.S1.L2': invalid sex 'unknown' (accepted: f, female, m, male)

bad platform:
Invalid samplesheet input:
- sample 'PT3.S1.L2': invalid platform 'illumina' (accepted: ont, pacbio, pb)

sex differs between tumor and normal:
Invalid samplesheet input:
- patient 'P3': sex mismatch, tumor 'PT3.S1.L2' is female but normal 'P3_PBMC' is male
```

## Before registering this process

Already done: the fork (branch `cirro`, commit `42458112ee927a37c1879a3f23e0fb82ab7b8fdb`) is pushed to the public
GitHub repository `yannvrb/lrsomatic` (the `code.uri` in `process-definition.json`), and this config repo is pushed to
`yannvrb/cirro-lrsomatic`. In a local clone of the fork, `origin` may still point at upstream
IntGenomicsLab/lrsomatic; do not push there.

What remains; nothing below has been done:

1. Replace the `documentationUrl` anchor in `process-definition.json`; it is an invented guess.
2. Confirm the references are staged (see "Prerequisites") at the paths in `process-input.json`, or change the paths to
   wherever they were put:
   - `$.referencesBucketUri|/sigprofiler/volume/`
   - `$.referencesBucketUri|/vep_revel/revel_grch38.tsv.gz`
   - `$.referencesBucketUri|/vep_revel/revel_grch38.tsv.gz.tbi`
3. Confirm in Cirro that a dataset uploaded as "Unaligned BAM" is offered this pipeline (`parentProcessIds` lists
   `ingest_unaligned_bam`, read from the catalogue and not run; see "Process definition notes").
4. Set `resourceLimits` in `process-compute.config` to the real AWS Batch limits.
5. Decide whether to keep `PW_ONDEMAND_JOB_QUEUE` for `MINIMAP2_ALIGN` and `LONGPHASE_HAPLOTAG`.
