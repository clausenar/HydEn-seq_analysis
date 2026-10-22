# Glue job: bedGraphs → Iceberg

`bedgraph_to_iceberg.py` is an AWS Glue 5.0 Spark job (Spark 3.5, Iceberg 1.7.1).
It turns the bedGraph text files the pipeline publishes into two Iceberg tables
in the stage's Glue database, which Athena (and any Iceberg-aware engine) can query.

```
tables/bedgraph/{sample}__{strand}.bedgraph      text, one file per sample and strand (Batch job)
        │
        ▼  Glue job, 2 × G.1X
hyden_seq_<stage>.bedgraph_signal   silver: sample_name, strand, chrom, start_pos, end_pos, score,
                                    source_file, processed_at
hyden_seq_<stage>.strand_balance    gold: per sample and chromosome, forward/reverse score and
                                    forward_share
        │
        ▼
Athena (saved query strand_balance_by_sample), files under warehouse/ in the stage's bucket
```

| | |
|---|---|
| Partitioning | both tables by `sample_name` |
| Reruns | `overwritePartitions()`: a rerun replaces that sample's rows, other samples are untouched |
| Bad input | stops before writing: unparseable lines (FAILFAST), files not named `{sample}__{forward\|reverse}.bedgraph`, no rows |
| Data Quality | `SILVER_RULESET` (Glue Data Quality, DQDL) runs on the silver rows before any write: completeness, strand values and both strands present, the 17 sacCer3 chromosomes, start ≥ 0, score > 0, end > start. A failed rule fails the job and nothing is written. Results: the job's *Data quality* tab and CloudWatch metrics (namespace `Glue Data Quality`). |
| History | snapshots older than 7 days are expired at the end of each run, keeping at least the last 5 |
| Table schema | created by the job (`CREATE TABLE IF NOT EXISTS`); Terraform owns the database, job, role and script upload |
| Logs | `/aws-glue/jobs/hyden-seq-<stage>-role/hyden-seq-<stage>-glue-job/output` (driver stdout, including the run's summary line) and `/error` (Spark logs), one stream per run and executor, kept 30 days. Glue names them after the security configuration and role. |

## Running it

Terraform uploads the script and creates the job (`terraform/modules/pipeline/glue_job.tf`).
It runs by itself for every sample the pipeline finishes: a successful Batch job
starts the Step Functions state machine `hyden-seq-<stage>-to-iceberg`
(`orchestration.tf`), which runs the job for that sample and sends an SNS alert
naming the sample if it fails. Start it by hand to reload everything or retry a sample:

```bash
# all samples in tables/bedgraph/
aws glue start-job-run --job-name hyden-seq-dev-bedgraph-to-iceberg

# one sample
aws glue start-job-run --job-name hyden-seq-dev-bedgraph-to-iceberg --arguments '{"--sample_name":"SRR1609190sub5"}'
```

A run costs about 2 DPU × a few minutes, roughly $0.05 (Glue bills $0.44 per
DPU-hour, per second with a 1-minute minimum).

## Tests

`tests/` runs the real job code on a local Spark with a file-based Iceberg
catalog instead of the Glue Data Catalog. Glue Data Quality only exists inside
Glue, so the tests check that the gate runs before any write and that a failed
gate leaves the tables unchanged; the ruleset itself is checked in a dev run. Needs Java 17, and network access the
first time (Spark downloads the Iceberg runtime jar). CI runs them in
`.github/workflows/glue-ci.yml`.

```bash
python3 -m venv .venv && .venv/bin/pip install -r glue/requirements-dev.txt
SPARK_LOCAL_IP=127.0.0.1 .venv/bin/pytest glue/tests -q
```
