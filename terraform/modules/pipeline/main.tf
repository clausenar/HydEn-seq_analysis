# HydEn-seq pipeline infrastructure as code (Terraform).
#
# The same event-driven pipeline as the CDK stack in infra/ (and the one first
# provisioned by hand, see the "AWS deployment" section of the top-level README):
#
#   S3 landing/{sample}/*.fastq  --(ObjectCreated)-->  Lambda landing watcher
#     --(DynamoDB conditional claim)-->  batch:SubmitJob
#     -->  AWS Batch on Fargate running the pipeline image from ECR
#     -->  output/{sample}/ in S3, raw input archived to processed_raw/{sample}/
#
# The container image is not built here: GitHub Actions builds it and pushes it
# to the existing ECR repository, which this module only looks up.

locals {
  name = "hyden-seq-${var.stage}"

  # Production keeps its data: the bucket refuses to be deleted while it holds
  # objects and the table has deletion protection. Any other stage is
  # disposable, so `terraform destroy` leaves nothing behind (and nothing billed).
  is_prod = var.stage == "prod"

  landing_prefix       = "landing/"
  processed_raw_prefix = "processed_raw/"

  # Queryable copies of the bedGraphs (see glue.tf) and Athena's result files.
  bedgraph_table_prefix = "tables/bedgraph/"
  athena_results_prefix = "athena-results/"
}

data "aws_region" "current" {}

data "aws_caller_identity" "current" {}

data "aws_ecr_repository" "pipeline" {
  name = var.ecr_repository_name
}
