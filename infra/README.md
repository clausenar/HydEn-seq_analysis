# HydEn-seq infrastructure (AWS CDK, Python)

The whole AWS pipeline as code: the same event-driven design that was first
set up by hand (see "AWS deployment" in the top-level README), now reviewable,
testable and reproducible with one command.

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

## What the stack creates

| Resource | Notes |
|---|---|
| S3 bucket | Private, encrypted, HTTPS only; `processed_raw/` moves to Infrequent Access after 30 days |
| DynamoDB table | `sample_id` key, on-demand billing; point-in-time recovery in prod |
| VPC | 2 public subnets, **no NAT gateway** (saves ~$30/month), free S3 gateway endpoint |
| Batch | Fargate compute environment (max 16 vCPU, scales to zero), job queue, job definition |
| Lambda | `lambda/landing_watcher_lambda.py`, triggered only by uploads under `landing/` |
| IAM | Least-privilege roles generated from code (tests enforce them) |
| Alerting | SNS topic, EventBridge rule for failed jobs, CloudWatch alarm on Lambda errors |

The container image is **not** built here. `.github/workflows/build-and-push-image.yml`
builds it and pushes it to the existing ECR repository `hyden-seq-pipeline`, which the stack imports.

## Stages

Every stage is a full, separate copy: `-c stage=dev` (default) or `-c stage=prod`.
Non-prod stages are disposable: `cdk destroy` deletes the bucket contents and the table.
Prod keeps both if the stack is ever deleted.

## Use

```bash
cd infra
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt -r requirements-dev.txt
npm install -g aws-cdk

.venv/bin/python -m pytest -q          # template tests, no AWS access needed
cdk synth                              # CloudFormation template, no AWS access needed
cdk diff  -c stage=dev --profile ec2-pipeline
cdk deploy -c stage=dev -c alert_email=you@example.com --profile ec2-pipeline
cdk destroy -c stage=dev --profile ec2-pipeline
```

The first deploy to the account needs a one-time `cdk bootstrap aws://069509443906/us-east-1 --profile ec2-pipeline`.
Before a dev test run, copy `reference/` from the live bucket into the dev bucket.

## CI

`.github/workflows/infra-ci.yml` runs the Lambda unit tests, the CDK template
tests and `cdk synth` on every pull request that touches `infra/` or `lambda/`.
Deploying stays a manual step.
