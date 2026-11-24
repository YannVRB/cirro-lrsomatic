# Launching LRsomatic from the Cirro web form

For lab members. This describes how to run the tumor/normal long-read somatic pipeline
(IntGenomicsLab/lrsomatic, our pinned fork) from Cirro.

Status: the process has been run end to end on Cirro, but only through the SDK. Nobody has used the web form yet.
The field names below come from the form definition (`community/lrsomatic/1.0/process-form.json`); the menu wording
and button names are Cirro's and may differ. If something does not match, tell Yann.

## Where

- Cirro instance: `breakthroughcancer.cirro.bio`
- Project: **BTC-GBM-Development**. The process is registered there only, not in BTC-GBM-Production.
- Process name: **Somatic Variant Calling for Long Reads (LRsomatic)**

## What you need first

1. **Two datasets of type "Unaligned BAM"** in the project: one for the tumor, one for the normal. The tumor and
   normal can be separate datasets. The BAM file name must end in `.bam`.
2. **Sample metadata on each sample**, filled in before the dataset reaches COMPLETED (a sample sheet added later does
   not fill the sample registry; ask Yann if that happened):

| Column | Values |
|---|---|
| `patient` | any text, the same for the tumor and its normal |
| `status` | `tumor` or `normal` |
| `sex` | `male` or `female` (same for the tumor and its normal) |
| `platform` | `ont` (same for the tumor and its normal) |
| `fiber` | `y` or `n` (the tumor's value is used) |

3. **Rules the pipeline checks at the start.** A mistake stops the run within a minute and the error names the sample.
   - One BAM per sample (no merging).
   - Every patient has exactly one normal.
   - Sample names contain no spaces.
   - Tumor-only runs are not supported.

## Launching

1. Open the project and select the tumor dataset and the normal dataset together as the inputs.
2. Choose **Run analysis** and pick the LRsomatic process.
3. Fill in the form (below), name the run, and start it.

### Form fields

| Field | What to do |
|---|---|
| Reference Genome | GRCh38 is the only option. Leave it. |
| Skip options (ASCAT, Savana, Wakhan, VEP, VEP plugins, Fiber-seq, Fiber-seq on normals, m6A, modification calling, QC) | Leave everything unticked for a full analysis. Tick a box to skip that step. |
| Skip mutational signatures | Leave unticked to run SigProfiler. |
| VEP cache version | Leave at `113`. |
| Extra Parameters (JSON) | Leave empty. A JSON object such as `{"skip_dmr": true}` overrides other pipeline parameters. Reference paths, the input and the output location cannot be overridden and are rejected. |

Testing on a small slice? Tick **Skip ASCAT**, **Skip Savana** and **Skip Wakhan**. They fail on a few megabases of data
because there are too few reads per window.

## How long it takes and what it costs

- Full-size ONT P4 (82 GB tumor and 183 GB normal BAMs, everything on): about 20 hours.
- chr4 subset (a few hundred MB): 12 to 19 minutes.
- Cost: the full no-cache P4 run cost $75.84 of EC2 plus $2.22 of EBS on the AWS account (November 2026), so plan on
  about $75 to $80 for a full-size run like it. Cirro's own estimate for that run was $25.12, about a third of the real
  cost, so do not trust it for full-size runs. The AWS bill goes to the project, and the project has a monthly budget.
  Check with Yann before launching several full-size runs.
- Cirro reuses finished tasks from earlier runs of the same process with the same inputs and parameters. A rerun after a
  failure is therefore much faster, and a test run on the same data may complete in minutes.

## Where the results are

The run creates a new dataset. In it, under `data/<sample>/`:

- `variants/phased/somatic_smallvariants.vcf.gz` (+ `.tbi`): phased somatic small variants
- `vep/somatic/<sample>_SOMATIC_VEP.vcf.gz` (+ `.tbi`): the same, VEP-annotated (absent if VEP was skipped)
- `report/<sample>_report.html`: the summary report

Everything else the pipeline writes (ASCAT, Wakhan, Savana, Severus, QC, methylation, SV and germline VEP files,
signatures) is in the same dataset but not highlighted. `config/pipeline_samplesheet.csv` shows the tumor and normal
pairing the pipeline was given. Check it if the pairing looks wrong.

## If a run fails

- **Fails within a minute:** a metadata or file problem. Open the run log: the preprocess step prints the samplesheet
  and a message naming the sample or patient.
- **A task died because of an AWS host loss or ran out of memory:** it is retried automatically, at most twice. A host
  loss is retried at the same size, an out-of-memory kill at a bigger size. You do not need to do anything.
- **A task failed with an ordinary error:** it is not retried and the run stops. Send Yann the run name and the run log.
- If the run is hung or you need to stop it, tell Yann. Stopping and restarting a full-size run wastes money.

## Not supported or not tested

- Tumor-only runs: not supported.
- Reference genomes other than GRCh38: not offered.
- PacBio data (`platform` `pb`): the form accepts it but it has never been run.
- Mixed ONT and PacBio samples: never run.
- The pipeline's `--plot` option: not exposed.
