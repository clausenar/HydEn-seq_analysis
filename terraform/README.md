# HydEn-seq infrastructure (Terraform)

The same event-driven pipeline as the CDK stack in [`infra/`](../infra/README.md),
written in Terraform. Both describe the same design; use one or the other for a
given set of resources, never both.

**This is what runs.** Since 2026-10-15 both stages are deployed from here in
account `069509443906`, `us-east-1`: `prod` is the live pipeline and `dev` is a
permanent test copy next to it. `bootstrap` has been applied, and the CDK version
is not deployed. The hand-built pipeline this replaced was retired on 2026-10-16
(see "AWS deployment" in the [top-level README](../README.md#aws-deployment)).

| | dev | prod |
|---|---|---|
| Data bucket | `hyden-seq-dev-data-90322281e63f81697d01cae60c` | `hyden-seq-prod-data-49d313619a5b9de291931f8e93` |
| Other resources | `hyden-seq-dev-…` | `hyden-seq-prod-…` |
| Terraform state | `s3://hyden-seq-tfstate-069509443906/pipeline/dev.tfstate` | `…/pipeline/prod.tfstate` |

```
S3 landing/{sample}/*.fastq ──(ObjectCreated)──► Lambda landing watcher
                                                  │  DynamoDB conditional claim (exactly once)
                                                  ▼
                                 AWS Batch on Fargate (8 vCPU / 16 GB, image from ECR)
                                                  │  Snakemake
                                                  ▼
                     S3 output/{sample}/   +   raw input archived to processed_raw/{sample}/
                                    └─► tables/bedgraph/{sample}__{strand}.bedgraph ◄── Glue table ◄── Athena
                                              │
                                              ▼  Glue Spark job, started per sample by Step Functions
                                              │  (Data Quality rules gate the write; SNS alert on failure)
                                    Iceberg: bedgraph_signal, strand_balance ◄── Athena

Failed Batch job ──► EventBridge ──► SNS alert      Lambda errors ──► CloudWatch alarm ──► SNS alert
Finished Batch job ──► EventBridge ──► landing watcher ──► DynamoDB status SUCCEEDED / FAILED
```

## Layout

```
terraform/
  modules/pipeline/      the whole pipeline as one module
    storage.tf             S3 data bucket (private, encrypted, HTTPS only, lifecycle) + DynamoDB status table
    network.tf             VPC, 2 public subnets, no NAT gateway, free S3 gateway endpoint
    batch.tf               Fargate compute environment, job queue, job definition
    lambda.tf              landing watcher (zipped from ../lambda), S3 trigger on landing/ only, job-finished rule
    iam.tf                 least-privilege roles, written out in full
    alerting.tf            SNS topic, EventBridge rule for failed jobs, Lambda error alarm
    glue.tf                Glue database + bedGraph table, Athena workgroup and saved queries
    glue_job.tf            Glue Spark job (script from ../glue), its role, security configuration, log groups
    orchestration.tf       Batch job succeeded ─► EventBridge ─► Step Functions ─► Glue job; SNS alert on failure
    snowflake.tf           read-only IAM role for Snowflake's Glue catalog integration and external volume
    athena/                the saved queries' SQL
    tests/                 `terraform test` suite (mocked AWS provider)
  envs/dev, envs/prod    one root per stage, each with its own state file
  bootstrap/             one-time: state bucket + read-only CI plan role
```

## Stages

Each stage folder calls the module with its own `stage`. Every resource name
starts with `hyden-seq-<stage>-`, so dev runs next to prod without clashing.

| | dev | prod |
|---|---|---|
| Pipeline image | follows `:latest`, so every merge that rebuilds the image reaches dev at once | pinned with `image_digest` in `envs/prod/main.tf`; changes only through a release (below) |
| `terraform destroy` | deletes everything, bucket contents included | refuses while the bucket holds data; table has deletion protection |
| DynamoDB point-in-time recovery | off | on |

Keep dev deployed: the CI `plan-dev` job compares every pull request against it.
It costs well under $1 a month while idle (no NAT gateway, Batch scales to zero).

## Use

Set both of these first. The profile matters because the default AWS profile on
the development Mac is a different user. The alert email isn't committed: a
`plan` or `apply` without it removes the SNS email subscription.

```bash
export AWS_PROFILE=ec2-pipeline
export TF_VAR_alert_email=<the address that gets pipeline alerts>
```

```bash
# tests: no AWS access needed
cd terraform/modules/pipeline && terraform init -backend=false && terraform test

# change a stage: always plan to a file, read it, then apply that file
cd terraform/envs/dev        # or envs/prod
terraform init
terraform plan -out=dev.tfplan
terraform apply dev.tfplan

# already done once per account: state bucket + CI plan role (local state)
terraform -chdir=terraform/bootstrap init
terraform -chdir=terraform/bootstrap apply
```

A new stage needs `reference/` copied into its bucket before its first run (dev
and prod already have it).

### Releasing an image to prod

1. Merge the change. The image workflow pushes a new `:latest`, which dev runs.
2. Test it in dev (below).
3. Copy the digest from the workflow run's summary into `image_digest` in
   `envs/prod/main.tf`, open a pull request, merge it, then plan and apply
   `envs/prod`. The plan shows 3 changes: the job definition gets a new
   revision, and the Lambda's `BATCH_JOB_DEFINITION` and its `batch:SubmitJob`
   permission move to that revision's ARN.

### Testing a stage end to end

A small paired test sample is kept in the prod bucket at
`test-data/SRR1609190sub/`: the first 500,000 read pairs of SRR1609190, about
21 MB and 23 MB. Nothing under `test-data/` triggers the pipeline. Each mate is
over 8 MB on purpose, so `aws s3 mv` copies it in parts and the archive step's
tagging permissions are exercised.

```bash
P=s3://hyden-seq-prod-data-49d313619a5b9de291931f8e93
L=s3://<stage bucket>/landing/<name>   # <name>: new each time, e.g. SRR1609190sub5
aws s3 cp $P/test-data/SRR1609190sub/SRR1609190sub_end1.fastq.gz $L/<name>_end1.fastq.gz
aws s3 cp $P/test-data/SRR1609190sub/SRR1609190sub_end2.fastq.gz $L/<name>_end2.fastq.gz
```

Use a name that hasn't been used before: one whose status record is still
`SUBMITTED` is not claimed again. A pass takes about 2.5 minutes: the job
`SUCCEEDED`, the status record says `SUCCEEDED`, `output/<name>/` holds 21
files, `tables/bedgraph/` holds `<name>__forward.bedgraph` and
`<name>__reverse.bedgraph`, both FASTQs are in `processed_raw/<name>/` and
`landing/<name>/` is empty. In prod, delete the test's `output/`,
`processed_raw/`, its two `tables/bedgraph/` files and its status record
afterwards so the test data isn't loaded into Snowflake or counted by Athena.
A successful run also starts the to-Iceberg state machine; in prod, delete the
test sample's partition from the two Iceberg tables too
(`DELETE FROM bedgraph_signal WHERE sample_name = '<name>'` in Athena, and the
same for `strand_balance`).

### Querying with Glue and Athena

Each stage has a Glue database (`hyden_seq_dev`, `hyden_seq_prod`) with one
table, `bedgraph_signal_raw`: chrom, start_pos, end_pos, score, the same
columns as `BEDGRAPH_SIGNAL_RAW` in Snowflake. It reads `tables/bedgraph/`,
not `output/`, because Athena reads every file under a table's location and
`output/` also holds BAMs and plots. The Batch job copies each sample's two
bedGraphs there (`HYDEN_S3_BEDGRAPH_TABLE_PREFIX` on the job definition). The
sample and strand are in the file name, so queries read them from `"$path"`,
as `stg_bedgraph_signal` does in dbt. There are no partitions: the data is a
few tens of MB per sample.

Query in the stage's Athena workgroup (`hyden-seq-dev`, `hyden-seq-prod`).
It enforces its settings: results go to `athena-results/` (SSE-S3, deleted
after 7 days) and a query stops after scanning 1 GB. The console lists two
saved queries, `bedgraph_signal_summary` (the same columns as the dbt mart, to
compare with Snowflake) and `bedgraph_strand_share`.

```bash
aws athena start-query-execution --work-group hyden-seq-prod \
  --query-execution-context Database=hyden_seq_prod \
  --query-string 'SELECT count(*) FROM bedgraph_signal_raw'
```

Samples processed before the copy step existed are not in `tables/bedgraph/`;
copy their two bedGraphs from `output/<sample>/` once.

The Glue job `hyden-seq-<stage>-bedgraph-to-iceberg` turns those files into two
Iceberg tables in the same database, `bedgraph_signal` (silver) and
`strand_balance` (gold), both partitioned by sample, with their files under
`warehouse/`. It runs for each sample the pipeline finishes (state machine
`hyden-seq-<stage>-to-iceberg`), and by hand for all samples or a retry; see
[`glue/README.md`](../glue/README.md). The saved query `strand_balance_by_sample`
reads the gold table.

## Tests and CI

`modules/pipeline/tests/pipeline.tftest.hcl` mirrors the CDK template tests
(private bucket, lifecycle, table contract, prod vs dev, Fargate sizing, network,
trigger prefix, least-privilege IAM, alerting, Glue table prefix and columns,
Athena workgroup limits, Glue job script, version, inputs and role, Glue log
groups, the to-Iceberg trigger, state machine and roles, variable validation). The Glue job's own tests are in `glue/tests` (see `glue/README.md`). It runs
against `mock_provider "aws"`, so it creates nothing and needs no credentials.

`.github/workflows/terraform-ci.yml` runs on every pull request that touches
`terraform/` or `lambda/`: `fmt -check`, `validate`, `terraform test`, `tflint`
and `checkov`. Checks deliberately skipped are listed with reasons in
[`.checkov.yaml`](.checkov.yaml). Set the repository variable `TF_PLAN_ROLE_ARN`
(output of `bootstrap`) to also get a read-only `terraform plan` for dev on each
pull request, and the repository secret `ALERT_EMAIL` so that plan does not show
the alert email subscription being destroyed. Applying stays manual.

## Differences from the CDK version

- Resources get readable names (`hyden-seq-dev-…`) instead of CloudFormation-generated ones.
- Batch tasks may only make outbound HTTPS (port 443); CDK allows all outbound traffic.
- The VPC's default security group is emptied.
- Prod protection uses `force_destroy = false` and DynamoDB deletion protection
  (Terraform's `prevent_destroy` cannot depend on a variable).
- Prod is pinned to one image digest (`image_digest`); the CDK version always runs `:latest`.
