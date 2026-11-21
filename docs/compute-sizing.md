## Compute

`community/lrsomatic/1.0/process-compute.config` keeps the pipeline's own per-label
resources (`conf/base.config` in the fork) and changes only what the execution traces of two
real Slurm runs justify.

### What the config sets

| Setting | Value | Why |
|---|---|---|
| `errorStrategy` / `maxRetries` | pipeline predicate plus unknown exit status (retry on exit 130-145, 104, 175-177 or `Integer.MAX_VALUE`, else `finish`) / 2 | The pipeline's predicate from the fork's `conf/base.config` is kept exactly and one condition is added, `task.exitStatus == Integer.MAX_VALUE` (Nextflow's value when no exit code was recorded). A lost host on AWS Batch does surface that way instead of as a container exit code (exit `-` in the trace of run b7730d18, 2026-11-17). Real non-zero exit codes outside the list are still not rerun. `maxRetries` is raised from the pipeline's 1 to 2 (three attempts), which also bounds retries of an unknown exit status. `aws.batch.maxSpotAttempts` is not set because its behaviour on Cirro cannot be verified. Overridden memory scales with `task.ext.retryScale`, so only after a kill (see the retry scaling row): `SAVANA_RUN` 220 / 440 / 488 GB (third attempt capped by `resourceLimits`), `SEVERUS` 96 / 192 / 288 GB; after a lost host the next attempt keeps the same size. |
| `resourceLimits` | `cpus: 64, memory: 488.GB, time: 240.h` | Sarek's values, checked against the on-demand compute environment on 2026-11-21 (screenshot `notes/Cirro_Compute_Environment_Overview.png`). It lists 9 families: c5ad, c5d, c6id, m5ad, m5d, m6id, r5ad, r5d, r6id. The largest sizes of the memory-heavy ones (r5d and r5ad 24xlarge: 96 vCPU, 768 GiB; r6id 32xlarge: 128 vCPU, 1024 GiB) can hold 64 cpus and 488 GiB. Not verified: the instance sizes are from the AWS catalogue, not from AWS for this account, and the environment's Max vCPUs (cut off in the screenshot) is unknown. If it is below about 96, a 64-cpu retry could wait in the queue. The largest measured peak is 80 GiB (SAVANA_RUN, run 6774e393), so the cap only matters on a retry after a memory kill. |
| retry scaling | `ext.retryScale`: `task.attempt` if the previous attempt exited 130-145, 104 or 175-177, else 1; multiplies cpus, memory and time of the 7 labels the workflow uses | `conf/base.config` multiplies by `task.attempt` on every retry, so a lost host (exit unknown) also asked for a bigger instance: in run c6012ed5 CLAIR3 went 24 -> 48 -> 64 cpus and 144 -> 288 -> 432 GB, and its peak memory grew with the threads (63.9 GB on Slurm, 170.9 GB at 64 cpus). Tested on Nextflow 25.10.4 with probe processes (attempt 1 base size, after exit 137 x2, after exit 1 back to x1); a real lost host (exit `-`) is not simulated. Known gap: after a memory kill followed by a lost host the next attempt is at the base size. |
| `SAVANA_RUN` memory | `220.GB * task.ext.retryScale` | Peak RSS 175.4 GB observed vs 72 GB from `process_high`. |
| `SEVERUS` memory | `96.GB * task.ext.retryScale` | Peak RSS 75.4 GB observed vs 36 GB from `process_medium`. |
| every task: queue | `env('PW_ONDEMAND_JOB_QUEUE')` | The first full-size P4 run (2026-11-17) lost 18 hosts on the default queue ("Host EC2 ... terminated"), restarting multi-hour tasks that each stage a 180 GB BAM; only `MINIMAP2_ALIGN` was on-demand and was not interrupted. The variable is Cirro-specific, copied from the sarek config. A plain process-level `queue = ...` was ignored in run c6012ed5 (all 36 re-run tasks ran on `Cirro-Spot-d29bb6ed`; the `queue` column of `artifacts/trace.txt` shows it), so the queue is set in a `withName: '.*'` block. That this block wins is not yet confirmed by a run. |

There is no `params {}` block (the pipeline requires parameters on the command line, never via
`-c`) and no Slurm or Singularity settings (Cirro supplies the AWS Batch ones).

### Evidence

Execution traces of the P4 run (three trace files, 15 + 35 + 159 rows) and the P3 run (one
non-empty trace, 109 rows), Slurm, about 200 GB input BAMs per sample. Both runs used the
pipeline's default resources: `P4.LRsomatic/lrsomatic_patched/conf/{base,modules}.config` are
identical to the fork's, and the `cluster.config` of both runs only sets executor, queue and
Singularity. Process names are stripped to the bare name used by `withName`. Nextflow trace
sizes are binary (1 GB = 1024^3 bytes). The traces cover 5 tumor/normal pairs (SAVANA_RUN and
SEVERUS show 7 and 5 tasks); this is a small sample.

Allocation is the first attempt from the process label in `conf/base.config` (no
`conf/modules.config` entry sets cpus, memory or time except `WGET_EVE`, never run). A process
is overridden only if its peak RSS exceeded the allocation, it exited 137 or 140, or its
realtime reached 75 % of the allocated time. Statuses seen: 265 COMPLETED, 46 CACHED, 4 FAILED.
No task exited 137 or 140. The FAILED tasks are: `SAVANA_RUN` twice (exit 143, about 20 min in,
on the largest sample; retried by the pipeline and then completed), `WGET` (exit 8) and
`SIGPROFILER_INSTALL` (exit 1), the last two not resource related.

| Process | Tasks | Max peak RSS (GB) | Max realtime (h) | Not COMPLETED/CACHED | Exit 137/140 | Label | Alloc mem (GB) | Peak / alloc | Realtime / alloc |
|---|---|---|---|---|---|---|---|---|---|
| SAVANA_RUN | 7 | 175.4 | 0.49 | 2 (exit 143) | 0 | process_high | 72 | 244 % | 3 % |
| SEVERUS | 5 | 75.4 | 2.69 | 0 | 0 | process_medium | 36 | 209 % | 34 % |
| LONGPHASE_PHASE_SOMATIC | 5 | 33.7 | 0.46 | 0 | 0 | process_medium | 36 | 94 % | 6 % |
| LONGPHASE_PHASE_GERMLINE | 5 | 27.7 | 0.17 | 0 | 0 | process_medium | 36 | 77 % | 2 % |
| NANOPLOT_POST | 11 | 8.2 | 0.02 | 0 | 0 | process_low | 12 | 68 % | 0 % |
| MINIMAP2_ALIGN | 11 | 64.5 | 4.24 | 0 | 0 | process_very_high | 144 | 45 % | 27 % |
| CLAIR3 | 5 | 63.9 | 1.83 | 0 | 0 | process_very_high | 144 | 44 % | 11 % |
| LRSOMATICREPORT | 5 | 15.0 | 0.07 | 0 | 0 | process_medium | 36 | 42 % | 1 % |
| SAVANA_CNA | 5 | 30.5 | 3.69 | 0 | 0 | process_high | 72 | 42 % | 23 % |
| LONGPHASE_MODCALL_SOMATIC | 5 | 29.1 | 1.41 | 0 | 0 | process_high | 72 | 40 % | 9 % |
| NANOPLOT_PRE | 16 | 4.7 | 0.01 | 0 | 0 | process_low | 12 | 39 % | 0 % |
| LONGPHASE_MODCALL_GERMLINE | 5 | 23.3 | 1.40 | 0 | 0 | process_high | 72 | 32 % | 9 % |
| MULTIQC | 2 | 1.9 | 0.01 | 0 | 0 | process_single | 6 | 32 % | 0 % |
| WHATSHAP_STATS | 5 | 1.3 | 0.04 | 0 | 0 | process_single | 6 | 22 % | 1 % |
| ASCAT | 5 | 7.2 | 0.40 | 0 | 0 | process_medium | 36 | 20 % | 5 % |
| WAKHAN | 5 | 7.2 | 0.93 | 0 | 0 | process_medium | 36 | 20 % | 12 % |
| MODKIT_PILEUP | 11 | 14.2 | 0.49 | 0 | 0 | process_high | 72 | 20 % | 3 % |
| LONGPHASE_HAPLOTAG | 10 | 7.3 | 5.35 | 0 | 0 | process_medium | 36 | 20 % | 67 % |
| CRAMINO_POST | 11 | 2.0 | 0.37 | 0 | 0 | process_low | 12 | 17 % | 9 % |
| CRAMINO_PRE | 16 | 1.3 | 0.41 | 0 | 0 | process_low | 12 | 11 % | 10 % |
| SAMTOOLS_STATS | 10 | 4.0 | 1.15 | 0 | 0 | process_medium | 36 | 11 % | 14 % |
| CLAIRS | 5 | 13.9 | 3.44 | 0 | 0 | process_very_high | 144 | 10 % | 21 % |
| SAVANA_CLASSIFY | 5 | 3.4 | 0.02 | 0 | 0 | process_medium | 36 | 9 % | 0 % |
| MOSDEPTH | 11 | 2.7 | 0.22 | 0 | 0 | process_medium | 36 | 8 % | 3 % |
| SOMATIC_VEP / SV_VEP / GERMLINE_VEP / VEP_SAVANA | 5 each | 3.0 / 2.4 / 3.4 / 1.2 | 0.07 / 0.04 / 0.84 / 0.02 | 0 | 0 | process_very_high | 144 | at most 2 % | at most 5 % |
| SAMTOOLS_INDEX, SAMTOOLS_IDXSTATS, SAMTOOLS_FLAGSTAT, BCFTOOLS_SORT, BCFTOOLS_CONCAT, BCFTOOLS_VIEW, SAMTOOLS_FAIDX, METAEXTRACT, UNZIP_* , UNTAR, VEPPLUGIN_REVEL | 3 to 22 | at most 0.5 | at most 0.43 | 0 | 0 | single / low / medium | 6 to 36 | at most 2 % | at most 11 % |
| WGET | 4 | 0.01 | 0.00 | 1 (exit 8) | 0 | process_single | 6 | n/a | n/a |
| SIGPROFILER_INSTALL | 1 | n/a | 0.02 | 1 (exit 1) | 0 | process_single + process_long | 6 | n/a | 0 % |

Headroom factor: **1.25** on the observed maximum, rounded up to the next multiple of 4 GB
(175.4 GB -> 220 GB, 75.4 GB -> 96 GB). Reason: the maxima come from 5 to 7 tasks on inputs of
about 200 GB, so a larger sample can exceed them, and AWS Batch enforces the container memory
limit as a hard OOM kill (on the Slurm cluster the same tasks ran above their requested memory
without being killed). A larger factor was not applied because the retry multiplier
(`* task.ext.retryScale`, which equals `task.attempt` after a kill and 1 otherwise, capped by `resourceLimits`, at most two retries) covers a miss. CPU and time are left at the label
defaults: no process reached 75 % of its allocated time.

### Left at defaults, worth knowing

- `LONGPHASE_PHASE_SOMATIC` used 94 % of its 36 GB (33.7 GB) and `LONGPHASE_PHASE_GERMLINE` 77 %.
  Not overridden because the peak did not exceed the allocation; a larger tumour may need the
  first retry.
- `LONGPHASE_HAPLOTAG` reached 67 % of its 8 h (5.35 h), below the 75 % threshold.
- `MINIMAP2_ALIGN` ran up to 4.24 h on Slurm. All tasks now use the on-demand queue (see the table above).

### No evidence (left at pipeline defaults)

These processes never ran in the P3/P4 traces: `SIGPROFILER_MATRIXGENERATOR`,
`SIGPROFILER_ASSIGNMENT` (signatures were skipped), `VEPPLUGIN_EVE`, `ENSEMBLVEP_DOWNLOAD`,
`CLAIRSTO`, `CLAIRSTO_CNA_RESOURCES`, `CLAIRSTO_VERDICT_TAG` and `SAVANA_TO` (tumour-only paths),
`DEEPVARIANT_*` and `DEEPSOMATIC_*`, `MODKIT_DMR`, `DMR_HAPLOTYPE_REGIONS`, `DMR_NEAREST_GENE`,
`HTSLIB_BGZIPTABIX`, `FIBERTOOLSRS_*`, `BCFTOOLS_NORM`, `BCFTOOLS_ISEC`, `BCFTOOLS_ANNOTATE`,
`BCFTOOLS_QUERY`, `BCFTOOLS_MERGE`, `SAMTOOLS_CAT`, `SAMTOOLS_MERGE`, `VCFSPLIT`,
`MINIMAP2_INDEX`, `PIGZ_UNCOMPRESS`.
`SIGPROFILER_INSTALL` ran once and failed (exit 1, not a resource problem).

### How the table was produced

A short pandas script over `P4.LRsomatic/pipeline_info/execution_trace_*.txt` and
`P3.LRsomatic/pipeline_info/execution_trace_*.txt` (columns `name, status, exit, realtime,
peak_rss`), grouped by the bare process name, compared with the label allocation in the fork's
`conf/base.config`. The script is not kept in the repository.
