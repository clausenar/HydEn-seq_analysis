# One-time setup that the pipeline stages depend on, applied by hand:
#
#   1. The S3 bucket that holds every stage's Terraform state.
#   2. A read-only IAM role GitHub Actions assumes (OIDC) to run
#      `terraform plan` on pull requests.
#
# This root keeps its own (tiny) state locally; it changes almost never.

terraform {
  required_version = ">= 1.10"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 6.0"
    }
  }
}

provider "aws" {
  region = "us-east-1"

  default_tags {
    tags = {
      Project   = "hyden-seq"
      ManagedBy = "terraform"
    }
  }
}

variable "github_repository" {
  description = "owner/repo allowed to assume the plan role from pull requests."
  type        = string
  default     = "clausenar/HydEn-seq_analysis"
}

data "aws_caller_identity" "current" {}

# ------------------------------------------------------------- state bucket

resource "aws_s3_bucket" "tfstate" {
  bucket = "hyden-seq-tfstate-${data.aws_caller_identity.current.account_id}"

  lifecycle {
    prevent_destroy = true
  }
}

resource "aws_s3_bucket_versioning" "tfstate" {
  bucket = aws_s3_bucket.tfstate.id

  # Lets a broken or deleted state file be recovered.
  versioning_configuration {
    status = "Enabled"
  }
}

resource "aws_s3_bucket_lifecycle_configuration" "tfstate" {
  bucket = aws_s3_bucket.tfstate.id

  # Every apply writes a new state version; keep old ones for 90 days.
  rule {
    id     = "ExpireOldStateVersions"
    status = "Enabled"

    filter {}

    noncurrent_version_expiration {
      noncurrent_days = 90
    }

    abort_incomplete_multipart_upload {
      days_after_initiation = 7
    }
  }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "tfstate" {
  bucket = aws_s3_bucket.tfstate.id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

resource "aws_s3_bucket_public_access_block" "tfstate" {
  bucket = aws_s3_bucket.tfstate.id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

# ----------------------------------------------------------- CI plan role

# The GitHub OIDC provider already exists in the account (the image build
# workflow uses it), so it is looked up rather than created.
data "aws_iam_openid_connect_provider" "github" {
  url = "https://token.actions.githubusercontent.com"
}

resource "aws_iam_role" "github_terraform_plan" {
  name        = "hyden-seq-github-actions-terraform-plan"
  description = "GitHub Actions: read-only terraform plan on pull requests"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Federated = data.aws_iam_openid_connect_provider.github.arn }
      Action    = "sts:AssumeRoleWithWebIdentity"
      Condition = {
        StringEquals = {
          "token.actions.githubusercontent.com:aud" = "sts.amazonaws.com"
          "token.actions.githubusercontent.com:sub" = "repo:${var.github_repository}:pull_request"
        }
      }
    }]
  })
}

# Plan has to read every resource it manages.
resource "aws_iam_role_policy_attachment" "plan_read_only" {
  role       = aws_iam_role.github_terraform_plan.name
  policy_arn = "arn:aws:iam::aws:policy/ReadOnlyAccess"
}

# Plan takes the S3 state lock, so it may write and delete only the .tflock files.
resource "aws_iam_role_policy" "plan_state_lock" {
  name = "terraform-state-lock"
  role = aws_iam_role.github_terraform_plan.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect   = "Allow"
      Action   = ["s3:PutObject", "s3:DeleteObject"]
      Resource = ["${aws_s3_bucket.tfstate.arn}/pipeline/*.tflock"]
    }]
  })
}

output "state_bucket" {
  value = aws_s3_bucket.tfstate.bucket
}

output "plan_role_arn" {
  description = "Set as the repository variable TF_PLAN_ROLE_ARN to turn on plans in CI."
  value       = aws_iam_role.github_terraform_plan.arn
}
