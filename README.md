# HydEn-seq analysis pipeline

A single Snakemake pipeline that goes from raw paired-end HydEn-seq/Ribo-seq
reads all the way to strand-specific bedGraphs, per-base ribonucleotide
counts, and origin-centered heatmaps/metaplots.

## Environment setup

The pipeline is split across **two conda environments** because the raw-read
processing tools (Snakemake, `bowtie`, `cutadapt`) need Python ≥3.11, while
`pybedtools`/`pysam` are only available as prebuilt Python 3.10 packages.
Don't try to merge them into one env — the dependency solver can't satisfy
both at once.

### `bio_env` — analysis and plotting (Python 3.10)

Used for `count_bases.py` and `origin_metaplot.py`.

```bash
conda create -n bio_env python=3.10
conda install -n bio_env -c bioconda -c conda-forge pybedtools pysam bedtools samtools
conda install -n bio_env -c conda-forge matplotlib
```

### `hyden_pipeline` — raw read processing (Python 3.11)

Used to run the `Snakefile`.

```bash
conda create -n hyden_pipeline -c bioconda -c conda-forge python=3.11 snakemake bowtie samtools bedtools perl-dbi perl-dbd-mysql
conda run -n hyden_pipeline pip install cutadapt
```

`cutadapt` is installed via `pip`, not `conda`, because the bioconda/conda-forge
build currently pulls a broken `xopen` package that fails to solve on macOS
(it spuriously requires a Windows-only virtual package). `pip install cutadapt`
inside the env sidesteps this.

### Verifying an environment

```bash
conda run -n bio_env python -c "import pybedtools, pysam, matplotlib; print('ok')"
conda run -n hyden_pipeline snakemake --version
conda run -n hyden_pipeline cutadapt --version
conda run -n hyden_pipeline bowtie --version
```

### Troubleshooting: `python` resolves to the wrong interpreter

If a conda env is active (prompt shows `(bio_env)`) but
`python -c "import sys; print(sys.executable)"` doesn't point inside that
env's `bin/`, check for a shell alias or function shadowing `python`
(`type python` in zsh/bash) — e.g. a stray `alias python=/usr/bin/python3`
in `~/.zshrc`. Remove it and open a new terminal.

## Data layout

Paths default to this machine's local layout but are all overridable via
`HYDEN_*` environment variables (see [AWS deployment](#aws-deployment) below
for why) — set one before running `snakemake`/the scripts, or edit the
`CONFIGURATION`/variable block defaults directly for a permanent local change.

| Path | Env var | Contents |
|---|---|---|
| `/Users/xranea/raw/` | `HYDEN_RAW_DIR` | Full-depth raw paired-end fastq files, named `{sample}_end1.fastq` / `{sample}_end2.fastq` |
| `/Users/xranea/raw_test/` | — | A small subsampled fastq pair for fast sanity-checking pipeline changes before running the full dataset |
| `/Users/xranea/genome/sacCer3*` | `HYDEN_GENOME` (prefix, no extension) | Bowtie1 index (`.ebwt`), reference fasta (`.fa`), and fasta index (`.fa.fai`) |
| `bin/oligo/oligo*` | `HYDEN_OLIGOS` (prefix, no extension) | Bowtie1 index used to filter out oligo-matching reads. Defaults to this repo's own `bin/oligo/oligo`, resolved relative to the Snakefile's location |
| `or200.txt` | `HYDEN_ORIGINS_FILE` | Replication origin coordinates (OriDB), used by `origin_metaplot.py`. Defaults to this repo's own copy |
| `/Users/xranea/bedgraphs/` | `HYDEN_OUT_DIR` | All pipeline intermediates and final bedGraphs |
| `/Users/xranea/bedgraphs/processed_results/` | (derived from `HYDEN_OUT_DIR`) | Output of `count_bases.py`; per-sample subfolders hold `origin_metaplot.py` output |

To switch which dataset the pipeline runs on locally, change `RAW_DIR`'s
default at the top of the `Snakefile` (e.g. back to `/Users/xranea/raw_test`
for a quick test run), or export `HYDEN_RAW_DIR` before running `snakemake`.

## Running the pipeline

The whole thing — raw reads through final plots — runs from one command:

```bash
conda activate hyden_pipeline
cd /Users/xranea/git/HydEn-seq_analysis
snakemake -n              # dry run — sanity-check the DAG first
snakemake --cores 8       # actually run it
```

`--cores 8` matches this machine's core count and what `align`'s `threads: 8`
declares (`bowtie -p{threads}`) — the only rule that's actually parallel;
everything else is single-threaded. With one sample the DAG is fully linear
anyway, so `--cores` mostly only matters for the align step's own thread count
and for scheduling multiple samples concurrently if you add more later.

### Pipeline stages

1. `cutadapt` trimming
2. oligo filtering (`bowtie`)
3. mate extraction (`seqkit pair`) — keeps R2 in sync with R1 after the
   oligo filter, which only runs on R1
4. alignment (`bowtie -S`, SAM output)
5. `samtools sort`
6. mate1 extraction (`samtools view -f 64 -F 4`)
7. `bedtools genomecov -bg -5 -strand +/-` → bedGraph, shifted 1bp to report
   the incorporated ribonucleotide's position (one base 5' of each read's
   mapped start) rather than the read start itself
8. `count_bases.py` — per-base counts (see below)
9. `origin_metaplot.py` — origin-centered heatmaps/metaplots (see below)

Steps 8 and 9 run in `bio_env`, not `hyden_pipeline`, even though the
`Snakefile` itself runs in `hyden_pipeline` — their shell commands go through
`conda run -n bio_env python ...` (`BIO_ENV_PYTHON` in the Snakefile), not a
direct path to `bio_env`'s python binary. That distinction matters:
`count_bases.py` uses `pybedtools`, which shells out to the `bedtools` CLI,
and a direct binary path skips `bio_env`'s own `PATH` setup, so `bedtools`
silently isn't found. `conda run` activates the environment properly.

### Rerunning after an edit

Snakemake only reruns rules whose *input file* mtimes changed — editing the
Snakefile's logic (e.g. swapping which fastq maps to read1) without any raw
file changing won't be detected automatically. If you edit the Snakefile,
force a rerun of the affected rules:

```bash
snakemake --cores 8 -R <rule_name>   # rerun one rule and everything downstream
snakemake --cores 8 -F               # force a full rerun from scratch
```

Careful with `-R` on a rule several steps downstream: most intermediates
through `bedgraph` are `temp()` and get deleted once a run completes
successfully. If you `-R` a late rule (e.g. `bedgraph`, `count_bases`,
`origin_metaplot`) after temp cleanup has already happened, Snakemake
doesn't just reuse the last non-temp file that still exists
(`{sample}_mate1.bam`, or the bedGraphs themselves) — it regenerates the
*entire* upstream chain from raw fastq, even though only the final step
actually needs to change. On a full-size dataset that's a ~20 minute rebuild
for what might be a one-line fix. If you only need to fix something at or
after the last non-temp output, it's faster to run that step's shell command
directly (find the exact command via `snakemake -n -p <rule_name>`) rather
than go through Snakemake's `-R`.

`bowtie2bedgraph_t1.pl` is the old Perl-based bedGraph generator this
Snakefile used before switching to `bedtools genomecov`. It's kept in the
repo for reference but is no longer called by the pipeline — it has a known
2bp coordinate bug on the reverse strand, so don't reuse it for new data.

`rectify_trimmed_pairs.pl` (re-syncing R1/R2 by read ID, plus an unconditional
1bp 3'-end trim) is likewise no longer called. Its pairing logic was already
redundant — `cutadapt`'s paired mode keeps R1/R2 in lockstep on its own
(verified by diffing read IDs) — and the 1bp trim has no effect on the
ribonucleotide's mapped position, which is read from mate1's 5' end.

### `count_bases.py` — bedGraph → per-base counts

Runs automatically as part of the pipeline (`rule count_bases`), or standalone:

```bash
conda activate bio_env
python count_bases.py
```

Reads every `*.bedGraph`/`*.bedgraph` in `INPUT_DIR`, groups
`{sample}__forward.bedgraph` / `{sample}__reverse.bedgraph` pairs by sample
name, and for each base position looks up the reference sequence
(`REFERENCE_FASTA`) via `pybedtools`. Writes, per sample, a per-position
detail file and a combined `base_count_totals.txt` summary (both strands
combined) to `OUTPUT_DIR`. No CLI args — edit the `CONFIGURATION` block at
the top of the script to point it elsewhere.

### `origin_metaplot.py` — bedGraph + origins → heatmap/metaplot

Runs automatically per sample as part of the pipeline (`rule origin_metaplot`,
writing to `processed_results/{sample}/`), or standalone:

```bash
conda activate bio_env
python origin_metaplot.py \
  --forward-bedgraph /path/to/{sample}__forward.bedgraph \
  --reverse-bedgraph /path/to/{sample}__reverse.bedgraph \
  --output-dir /path/to/output/dir \
  --origins-file or200.txt   # optional, defaults to or200.txt in this repo
```

All flags are optional and fall back to the `CONFIGURATION` block's defaults
at the top of the script if omitted.

Centers every origin in `--origins-file` (arabic `chr1..chr16` naming is
mapped to the roman-numeral `chrI..chrXVI` used by the bedGraphs) on its
midpoint, bins a ±`WINDOW` bp window (default 2000bp) into `BIN_SIZE` bp bins
(default 50bp), and produces two sets of outputs — one for all origins, one
restricted to origins with an assigned ARS name (suffixed `_named`, i.e.
excluding OriDB's unnamed/`null` entries):

- `origin_heatmap[_named].png` / `origin_metaplot[_named].png` — Watson (+)
  and Crick (-) strand signal, origins stacked as rows. The heatmap uses one
  blue (low) → red (high) color scale (one shared colorbar), but the
  gradient direction flips exactly at the origin midpoint — upstream reads
  blue→red, downstream reads red→blue — so the same count value renders as
  opposite colors on either side, making the origin boundary visually
  obvious. The Watson and Crick panels are colored as inverted mirrors of
  each other (where Watson is blue→red, Crick is red→blue, and vice versa),
  so a symmetric strand-switch pattern also shows up as a color inversion
  between the two panels
- `origin_ratio_heatmap[_named].png` / `origin_ratio_metaplot[_named].png` —
  log2(Watson/Crick) strand-bias ratio; the metaplot version is the clearest
  signal of an active origin — it shows a sharp sign flip (negative
  upstream, positive downstream) right at the midpoint, from the
  leading/lagging-strand switch as replication forks diverge
- `origin_matrix_watson[_named].tsv` / `origin_matrix_crick[_named].tsv` —
  the raw origin × bin matrices, for re-analysis without rerunning the script

### `braid_plot.py` — bedGraph → genome-wide strand-ratio "braid plot"

Runs automatically per sample as part of the pipeline (`rule braid_plot`,
writing to `processed_results/{sample}/`), or standalone:

```bash
conda activate bio_env
python braid_plot.py \
  --forward-bedgraph /path/to/{sample}__forward.bedgraph \
  --reverse-bedgraph /path/to/{sample}__reverse.bedgraph \
  --genome-fai /path/to/genome.fa.fai   # optional, defaults to $HYDEN_GENOME + .fa.fai
  --origins-file or200.txt              # optional, defaults to or200.txt in this repo
  --output-dir /path/to/output/dir \
  --bin-size 500                        # optional, defaults to 500bp
```

Unlike `origin_metaplot.py`, this isn't centered on origins — it bins the
*entire* genome into fixed-size windows (default 500bp) and plots
log2(Watson/Crick) per bin as one panel per chromosome, with each origin
from `--origins-file` marked as a vertical line - matching
`origin_metaplot.py`'s own Watson/Crick convention. The ratio's sign flips
across replication origins and termini as the leading/lagging strand
switches, so plotted as a continuous line it alternates above and below
zero in a repeating pattern — the "braid" this plot is named for.

- `braid_plot.png` — one horizontal panel per chromosome, each windowed to
  that chromosome's own length (and drawn at a proportional on-page width,
  so relative chromosome lengths stay visible at a glance rather than every
  panel sharing one genome-wide axis); red fill where Watson > Crick, blue
  where Crick > Watson; dashed vertical lines mark each origin's midpoint
- `braid_matrix.tsv` — the raw per-bin `chrom, bin_start, watson, crick,
  log2_watson_crick` values, for re-analysis without rerunning the script

## AWS deployment

The pipeline also runs on AWS Batch (Fargate), reading input and writing
output through S3, in AWS account `069509443906` (`us-east-1`, profile
`ec2-pipeline` locally). Everything below is already provisioned and live —
this section documents what exists and how to use/rebuild it.

### Architecture

```
S3 landing/{sample}/*.fastq[.gz]
        │  (S3 event on upload)
        ▼
Lambda: hyden-seq-landing-watcher (source: lambda/landing_watcher_lambda.py)
  - waits until both _end1 and _end2 exist for {sample}
  - DynamoDB conditional update (table hyden-seq-sample-status, key
    sample_id) so the two near-simultaneous upload events don't both
    submit a job - only the invocation that flips status to SUBMITTED
    proceeds
        │  batch:SubmitJob
        ▼
AWS Batch (queue: hyden-seq-job-queue, job def: hyden-seq-pipeline-job)
  Fargate task running the hyden-seq-pipeline:latest image (ECR)
  docker/batch_entrypoint.sh:
    1. sync s3://.../reference/        → /data/reference  (genome, oligo index, or200.txt)
    2. sync s3://.../landing/{sample}/ → /data/raw
    3. HYDEN_* env vars point the unmodified Snakefile at /data/*
    4. snakemake --cores 8
    5. sync /data/output → s3://.../output/{sample}/
    6. move s3://.../landing/{sample}/ → s3://.../processed_raw/{sample}/
```

### Provisioned resources

**Storage**

| Resource | Notes |
|---|---|
| S3 bucket `kunkel-ribo-data-2026` | All data — see prefixes below |
| ECR repo `hyden-seq-pipeline` | One combined image (both conda envs) |
| DynamoDB table `hyden-seq-sample-status` | Landing-watcher idempotency + job status. Partition key `sample_id`; attributes `status` (`PENDING`/`SUBMITTED`/`COMPLETED`/`FAILED`), `job_id`, `mate1_key`, `mate2_key`, `created_at`, `updated_at`. On-demand billing (matches the pipeline's scale-to-zero cost design) |

S3 prefixes:

| Prefix | Contents |
|---|---|
| `reference/` | Genome index, oligo index, `or200.txt` — synced to `/data/reference` at container startup. Update the reference by re-uploading here; no image rebuild needed |
| `landing/{sample}/` | Drop zone — upload a sample's paired fastqs here to trigger automatic processing |
| `processed_raw/{sample}/` | Raw fastqs moved here after a landing-triggered job completes |
| `output/{sample}/` | Pipeline output, same layout as local `OUT_DIR` |

**Compute**

| Resource | Notes |
|---|---|
| Batch compute environment `hyden-seq-fargate-ce` | Fargate, serverless, scales to zero — no cost when idle, max 16 vCPU |
| Batch job queue `hyden-seq-job-queue` | |
| Batch job definition `hyden-seq-pipeline-job` | 8 vCPU / 16GB, command = `[input_uri, output_uri, reference_uri?, processed_raw_uri?]` |
| Lambda `hyden-seq-landing-watcher` | Triggered by S3 `ObjectCreated` under `landing/`; auto-submits Batch jobs. Source: `lambda/landing_watcher_lambda.py` |

**IAM** (each role scoped to only what it needs)

| Role | Grants |
|---|---|
| `hyden-seq-batch-execution-role` | Pulls the ECR image, writes CloudWatch logs |
| `hyden-seq-batch-job-role` | Container's own S3 read/write, scoped to the bucket |
| `hyden-seq-batch-service-role` | Batch's own AWS-managed service role |
| `hyden-seq-landing-watcher-role` | Lambda's S3 `ListBucket`/`GetObject` on the bucket, `dynamodb:GetItem`/`PutItem`/`UpdateItem` scoped to `hyden-seq-sample-status`, and `batch:SubmitJob` |
| `hyden-seq-github-actions-ecr-push` | GitHub Actions' OIDC role for CI image builds — ECR push only, this repo's `main` branch only |
| OIDC provider `token.actions.githubusercontent.com` | Lets GitHub Actions assume the role above without any stored AWS credentials |

> **Rollout status:** live and cleaned up. The DynamoDB table is created,
> `hyden-seq-landing-watcher-role` has the IAM grants above (the old
> `s3:PutObject` lock-file permission has been removed), and
> `lambda/landing_watcher_lambda.py` is deployed as the function's code.
> Verified end-to-end against a real sample: the fixed code correctly
> claimed the sample via DynamoDB, submitted a Batch job, and the job
> succeeded, with output landing in `output/{sample}/` and raw input
> archived to `processed_raw/{sample}/`. Old `locks/{sample}.lock` objects
> have been deleted.
>
> Two issues only surfaced under this real test, not the moto-based unit
> tests, since they're AWS API validation rules rather than logic bugs:
> Batch job names reject characters like `.` (sample names may contain
> them), and `submit_job` rejects an empty string as a command array
> element. Both are fixed in `landing_watcher_lambda.py` and covered by
> regression tests in `lambda/test_landing_watcher_lambda.py`.

### What's automatic vs manual

- **Fully automatic**: merging code to `main` → image rebuild (CI); uploading a
  complete fastq pair to `landing/{sample}/` → job submission → output in S3
  → raw data archived to `processed_raw/`
- **Manual, by choice**: submitting a job directly via `aws batch submit-job`
  (for reprocessing existing data, or input that's already elsewhere in S3)

### Cost shape

Nothing runs, and nothing is billed, while idle. Ongoing cost is just S3 and
ECR storage (both cheap at this data scale) plus whatever Fargate compute a
job actually consumes while running. Lambda invocations and GitHub Actions
minutes stay comfortably inside the free tier at any normal usage frequency.

### Using the landing zone (fully automatic)

```bash
aws s3 cp my_sample_end1.fastq s3://kunkel-ribo-data-2026/landing/my_sample/my_sample_end1.fastq --profile ec2-pipeline
aws s3 cp my_sample_end2.fastq s3://kunkel-ribo-data-2026/landing/my_sample/my_sample_end2.fastq --profile ec2-pipeline
```

That's it — the Lambda picks up the second upload, submits a Batch job once
both mates are present, and the raw files end up in `processed_raw/my_sample/`
once it succeeds. Watch progress via:

```bash
aws logs tail /aws/lambda/hyden-seq-landing-watcher --since 5m --profile ec2-pipeline
aws batch list-jobs --job-queue hyden-seq-job-queue --job-status RUNNING --profile ec2-pipeline
```

### Submitting a job manually

Bypasses the landing zone — useful for reprocessing, or input that's already
in S3 somewhere else:

```bash
aws batch submit-job --profile ec2-pipeline --region us-east-1 \
  --job-name hyden-seq-my-sample \
  --job-queue hyden-seq-job-queue \
  --job-definition hyden-seq-pipeline-job \
  --container-overrides '{"command":["s3://kunkel-ribo-data-2026/some/input/prefix/","s3://kunkel-ribo-data-2026/output/my_sample/"]}'
```

The command list is `[input_uri, output_uri]`, with two optional trailing
args: a reference-data URI (defaults to `s3://kunkel-ribo-data-2026/reference/`)
and a processed-raw destination URI (if omitted, raw input is left in place).

### Rebuilding and redeploying the image

`docker/Dockerfile` builds by `git clone`-ing this repo at build time (arg
`GIT_REF`) rather than copying a local checkout — the image always runs the
exact code at a real, traceable commit. There is no "repo on AWS" to keep in
sync: edit code locally, open a PR, merge to `main` like any other change.

**This now happens automatically.** `.github/workflows/build-and-push-image.yml`
rebuilds and pushes the image to ECR (tag `:latest`) on every push to `main`
that touches pipeline code (`Snakefile`, `*.py`, `bin/`, `or200.txt`,
`docker/`) — pinned to the exact triggering commit SHA. It authenticates to
AWS via OIDC federation (role `hyden-seq-github-actions-ecr-push`, scoped to
this repo's `main` branch and ECR push on this one repository — no AWS keys
stored in GitHub), and uses GitHub Actions' cache backend for Docker layers,
so only the final `git clone` layer actually rebuilds on most runs rather
than reinstalling both conda environments from scratch. This is decoupled
from job execution: Batch jobs always run whatever image currently has the
`:latest` tag, regardless of how often data lands in the S3 landing zone or
how many jobs run — only a code change on `main` triggers a rebuild.

To trigger a rebuild without a code change (e.g. after editing the workflow
itself), use **Actions → Build and push pipeline image → Run workflow** on
GitHub, or `gh workflow run build-and-push-image.yml`.

To build manually instead (e.g. for local testing before merging):

```bash
cd docker
docker build --platform linux/amd64 -t hyden-seq-pipeline:latest -f Dockerfile .
# pin a specific commit instead of latest main:
#   docker build --platform linux/amd64 --build-arg GIT_REF=<commit-sha> -t hyden-seq-pipeline:latest -f Dockerfile .

aws ecr get-login-password --region us-east-1 --profile ec2-pipeline \
  | docker login --username AWS --password-stdin 069509443906.dkr.ecr.us-east-1.amazonaws.com
docker tag hyden-seq-pipeline:latest 069509443906.dkr.ecr.us-east-1.amazonaws.com/hyden-seq-pipeline:latest
docker push 069509443906.dkr.ecr.us-east-1.amazonaws.com/hyden-seq-pipeline:latest
```

New Batch jobs pick up the `:latest` tag automatically — no job definition
change needed unless resource requirements (vCPU/memory) change.

## Snowflake warehouse

Pipeline bedGraph output also lands in Snowflake, read directly from the
same S3 bucket the AWS pipeline already writes to — no data duplication,
no second copy step. This is a practice/portfolio addition on top of the
pipeline itself (the pipeline's own outputs don't depend on it), demonstrating
a cross-cloud storage integration pattern: Snowflake (account region
`AWS_US_WEST_2`) assumes an IAM role in the pipeline's AWS account to read
S3 directly, rather than data being pushed into Snowflake from the AWS side.

### Architecture

```
S3 kunkel-ribo-data-2026/output/{sample}/{sample}__forward.bedgraph
S3 kunkel-ribo-data-2026/output/{sample}/{sample}__reverse.bedgraph
        │
        │  read via cross-account IAM role assumption
        │  (external ID-scoped trust, Snowflake-initiated)
        ▼
Snowflake storage integration: HYDEN_S3_INT
        │
        ▼
External stage: HYDEN_SEQ.PIPELINE.HYDEN_OUTPUT_STAGE
        │  COPY INTO, pattern-matched on *__(forward|reverse).bedgraph
        │  (not hardcoded to specific sample names - covers any past or
        │  future pipeline output under output/, no re-run needed)
        ▼
HYDEN_SEQ.PIPELINE.BEDGRAPH_SIGNAL_RAW  (landing table, keeps METADATA$FILENAME
        │                                as src_file for provenance/re-parsing)
        ▼
HYDEN_SEQ.PIPELINE.BEDGRAPH_SIGNAL  (view: sample_name/strand parsed out of
                                      src_file via regex - the queryable table)
```

### Provisioned resources

| Resource | Notes |
|---|---|
| Warehouse `HYDEN_WH` | XSMALL, `AUTO_SUSPEND = 60` / `AUTO_RESUME = TRUE` - zero compute cost while idle, matching the AWS side's scale-to-zero design |
| Database/schema `HYDEN_SEQ.PIPELINE` | Holds everything below |
| Storage integration `HYDEN_S3_INT` | `STORAGE_ALLOWED_LOCATIONS` scoped to `s3://kunkel-ribo-data-2026/output/` only - can't reach `landing/`, `reference/`, or `processed_raw/` |
| File format `BEDGRAPH_TSV` | Tab-delimited, no header row (matches the pipeline's raw `bedtools genomecov` bedGraph output exactly) |
| External stage `HYDEN_OUTPUT_STAGE` | Points at the storage integration + `output/` prefix |
| Table `BEDGRAPH_SIGNAL_RAW` | Landing table: `chrom, start_pos, end_pos, score, src_file, loaded_at` |
| View `BEDGRAPH_SIGNAL` | `sample_name, strand, chrom, start_pos, end_pos, score, loaded_at` - `sample_name`/`strand` parsed from `src_file` via `REGEXP_SUBSTR` |

AWS side, in the pipeline's own account (`069509443906`):

| Resource | Notes |
|---|---|
| IAM role `snowflake-hyden-seq-s3-access` | Trust policy allows only Snowflake's generated IAM user (`STORAGE_AWS_IAM_USER_ARN` from `DESC INTEGRATION HYDEN_S3_INT`), conditioned on `sts:ExternalId` matching Snowflake's generated external ID - the standard Snowflake storage-integration security pattern, not a broad cross-account trust |
| Inline policy `s3-read-output` | `s3:GetObject`/`s3:GetObjectVersion` + prefix-scoped `s3:ListBucket`, both limited to `output/*` - read-only, can't write or touch other prefixes |

### What's automatic vs manual

- **Automatic**: the warehouse suspending/resuming itself; the view staying in sync with whatever's in the landing table
- **Manual**: re-running `COPY INTO BEDGRAPH_SIGNAL_RAW ...` after a new pipeline run lands new bedGraphs in S3. `COPY INTO` only loads files it hasn't already loaded (tracked via Snowflake's own load history), so it's safe to re-run repeatedly - it just does nothing on files already landed. There's no S3-event-triggered auto-load wired up (the AWS side's Lambda watcher triggers *pipeline runs*, not this load) - a natural next step would be an SNS/SQS-triggered Snowpipe for that, not yet built here

### Cost shape

Nothing runs, and nothing is billed for warehouse compute, while idle -
`HYDEN_WH` auto-suspends after 60 seconds. Ongoing cost is Snowflake account
credits only consumed during an actual `COPY INTO` or query (and even those
are usually a few seconds against XSMALL for this data volume), plus
whatever Snowflake's own storage cost is for the landed rows - not tracked
separately from the AWS side's own S3 storage cost, which is unaffected
(this reads S3, it doesn't copy data out of it).

### Querying

```sql
USE WAREHOUSE HYDEN_WH;
USE DATABASE HYDEN_SEQ;
USE SCHEMA PIPELINE;

SELECT sample_name, strand, COUNT(*) AS bins, SUM(score) AS total_score
FROM BEDGRAPH_SIGNAL
GROUP BY 1, 2
ORDER BY 1, 2;
```

### Loading new pipeline output

After a new sample's bedGraphs land in `s3://kunkel-ribo-data-2026/output/`:

```sql
COPY INTO BEDGRAPH_SIGNAL_RAW (chrom, start_pos, end_pos, score, src_file)
FROM (
    SELECT $1, $2, $3, $4, METADATA$FILENAME
    FROM @HYDEN_OUTPUT_STAGE
)
PATTERN = '.*__(forward|reverse)\\.bedgraph'
FILE_FORMAT = BEDGRAPH_TSV;
```

### Transforming the landed data: `dbt/hyden_seq/`

`BEDGRAPH_SIGNAL_RAW` above is the raw landing table; a small dbt project in
[`dbt/hyden_seq/`](dbt/hyden_seq/) turns it into a proper staging → marts
model (source declaration, tests, a SQL-native version of `braid_plot.py`'s
log2(Watson/Crick) ratio) instead of one hand-written view. See that
directory's own README for how to run it.
