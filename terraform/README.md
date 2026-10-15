# HydEn-seq infrastructure (Terraform)

The same event-driven pipeline as the CDK stack in [`infra/`](../infra/README.md),
written in Terraform. Both describe the same design; use one or the other for a
given set of resources, never both.

```
S3 landing/{sample}/*.fastq ──(ObjectCreated)──► Lambda landing watcher
                                                  │  DynamoDB conditional claim (exactly once)
                                                  ▼
                                 AWS Batch on Fargate (8 vCPU / 16 GB, image from ECR)
                                                  │  Snakemake
                                                  ▼
                     S3 output/{sample}/   +   raw input archived to processed_raw/{sample}/

Failed Batch job ──► EventBridge ──► SNS alert      Lambda errors ──► CloudWatch alarm ──► SNS alert
```

## Layout

```
terraform/
  modules/pipeline/      the whole pipeline as one module
    storage.tf             S3 data bucket (private, encrypted, HTTPS only, lifecycle) + DynamoDB status table
    network.tf             VPC, 2 public subnets, no NAT gateway, free S3 gateway endpoint
    batch.tf               Fargate compute environment, job queue, job definition
    lambda.tf              landing watcher (zipped from ../lambda), S3 trigger on landing/ only
    iam.tf                 least-privilege roles, written out in full
    alerting.tf            SNS topic, EventBridge rule for failed jobs, Lambda error alarm
    tests/                 `terraform test` suite (mocked AWS provider)
  envs/dev, envs/prod    one root per stage, each with its own state file
  bootstrap/             one-time: state bucket + read-only CI plan role
```

## Stages

Each stage folder calls the module with its own `stage`. Every resource name
starts with `hyden-seq-<stage>-`, so a dev copy can run next to the live pipeline.

| | dev | prod |
|---|---|---|
| `terraform destroy` | deletes everything, bucket contents included | refuses while the bucket holds data; table has deletion protection |
| DynamoDB point-in-time recovery | off | on |

## Use

Credentials come from your AWS profile, e.g. `export AWS_PROFILE=ec2-pipeline`.

```bash
# once per account: state bucket + CI plan role (local state)
terraform -chdir=terraform/bootstrap init
terraform -chdir=terraform/bootstrap apply

# tests: no AWS access needed
cd terraform/modules/pipeline && terraform init -backend=false && terraform test

# deploy a dev copy
cd terraform/envs/dev
terraform init
terraform plan  -var alert_email=you@example.com
terraform apply -var alert_email=you@example.com
terraform destroy
```

Before a dev test run, copy `reference/` from the live bucket into the dev bucket
(`terraform output bucket_name`), then upload a sample to `landing/{sample}/`.

## Tests and CI

`modules/pipeline/tests/pipeline.tftest.hcl` mirrors the CDK template tests
(private bucket, lifecycle, table contract, prod vs dev, Fargate sizing, network,
trigger prefix, least-privilege IAM, alerting, variable validation). It runs
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
