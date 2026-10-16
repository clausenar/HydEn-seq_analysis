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

  # Prod runs one exact image; dev keeps following :latest. To release a new
  # image, test it in dev, then copy its digest here (the Build and push
  # pipeline image workflow prints it in its run summary) and apply prod.
  # This one was built from a55b9ed and passed the dev run on 2026-10-16.
  image_digest = "sha256:5edb89829dbfe61eb7226f5c14971a9289898910d17d062c21deebf0d82afe05"
}
