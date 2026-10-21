terraform {
  required_version = ">= 1.10"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 6.0"
    }
    archive = {
      source  = "hashicorp/archive"
      version = "~> 2.7"
    }
  }

  # State bucket created once by terraform/bootstrap. use_lockfile gives
  # S3-native state locking, so no DynamoDB lock table is needed.
  backend "s3" {
    bucket       = "hyden-seq-tfstate-069509443906"
    key          = "pipeline/dev.tfstate"
    region       = "us-east-1"
    encrypt      = true
    use_lockfile = true
  }
}

provider "aws" {
  region = "us-east-1"

  default_tags {
    tags = {
      Project   = "hyden-seq"
      Stage     = "dev"
      ManagedBy = "terraform"
    }
  }
}

module "pipeline" {
  source = "../../modules/pipeline"

  stage       = "dev"
  alert_email = var.alert_email
}

# The first dev runs of the Glue job (2026-10-21) created its two log groups
# before Terraform managed them; this adopts them. A no-op once applied, and
# prod never had them, so it needs no import.
import {
  for_each = toset(["output", "error"])
  to       = module.pipeline.aws_cloudwatch_log_group.glue_job[each.key]
  id       = "/aws-glue/jobs/hyden-seq-dev-role/hyden-seq-dev-glue-job/${each.key}"
}
