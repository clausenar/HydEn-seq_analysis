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
    key          = "pipeline/prod.tfstate"
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
      Stage     = "prod"
      ManagedBy = "terraform"
    }
  }
}

module "pipeline" {
  source = "../../modules/pipeline"

  stage       = "prod"
  alert_email = var.alert_email
}
