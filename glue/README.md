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
| History | snapshots older than 7 days are expired at the end of each run, keeping at least the last 5 |
| Table schema | created by the job (`CREATE TABLE IF NOT EXISTS`); Terraform owns the database, job, role and script upload |

## Running it

Terraform uploads the script and creates the job (`terraform/modules/pipeline/glue_job.tf`).
In phase 2 it is started by hand; phase 3 will start it when a Batch job succeeds.

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
catalog instead of the Glue Data Catalog. Needs Java 17, and network access the
first time (Spark downloads the Iceberg runtime jar). CI runs them in
`.github/workflows/glue-ci.yml`.

```bash
python3 -m venv .venv && .venv/bin/pip install -r glue/requirements-dev.txt
SPARK_LOCAL_IP=127.0.0.1 .venv/bin/pytest glue/tests -q
```
