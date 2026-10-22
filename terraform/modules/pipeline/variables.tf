variable "stage" {
  description = "Deployment stage (dev, prod, ...). Every stage is a full, separate copy of the pipeline."
  type        = string

  validation {
    condition     = can(regex("^[a-z][a-z0-9-]{1,15}$", var.stage))
    error_message = "stage must be 2-16 characters: lowercase letters, digits and hyphens, starting with a letter."
  }
}

variable "ecr_repository_name" {
  description = "Existing ECR repository holding the pipeline image (built and pushed by GitHub Actions, not here)."
  type        = string
  default     = "hyden-seq-pipeline"
}

variable "image_tag" {
  description = "Image tag the Batch job definition runs. Ignored when image_digest is set."
  type        = string
  default     = "latest"
}

# CI pushes every merge to main as :latest and the repository's tags are
# mutable, so a tag can't hold a stage on a tested image. A digest can.
variable "image_digest" {
  description = "Optional image digest (sha256:...) to run instead of image_tag; pins the stage to one exact image."
  type        = string
  default     = null

  validation {
    condition     = var.image_digest == null || can(regex("^sha256:[0-9a-f]{64}$", var.image_digest))
    error_message = "image_digest must be sha256: followed by 64 lowercase hex characters."
  }
}

variable "alert_email" {
  description = "Optional email address subscribed to pipeline alerts. AWS sends a confirmation email first."
  type        = string
  default     = null
}

variable "availability_zones" {
  description = "Two AZs for the public subnets. Fixed here so plans and tests need no AWS lookup."
  type        = list(string)
  default     = ["us-east-1a", "us-east-1b"]

  validation {
    condition     = length(var.availability_zones) == 2
    error_message = "Exactly two availability zones are expected."
  }
}

variable "vpc_cidr" {
  description = "CIDR block for the pipeline VPC."
  type        = string
  default     = "10.0.0.0/16"
}

variable "max_vcpus" {
  description = "Upper limit for the Fargate compute environment. It scales to zero when idle."
  type        = number
  default     = 16
}

# Filled in after the Snowflake catalog integration and external volume exist
# (see snowflake.tf): DESC CATALOG INTEGRATION gives GLUE_AWS_IAM_USER_ARN and
# GLUE_AWS_EXTERNAL_ID, DESC EXTERNAL VOLUME gives STORAGE_AWS_IAM_USER_ARN and
# STORAGE_AWS_EXTERNAL_ID. Usually both IAM user ARNs are the same.
variable "snowflake_iceberg_trust" {
  description = "Snowflake IAM user ARN(s) and external IDs allowed to assume the Iceberg read role; null until linked."
  type = object({
    iam_user_arns = list(string)
    external_ids  = list(string)
  })
  default = null

  validation {
    condition = var.snowflake_iceberg_trust == null || (
      length(var.snowflake_iceberg_trust.iam_user_arns) > 0 &&
      length(var.snowflake_iceberg_trust.external_ids) > 0 &&
      alltrue([for a in var.snowflake_iceberg_trust.iam_user_arns : can(regex("^arn:aws:iam::[0-9]{12}:user/", a))])
    )
    error_message = "snowflake_iceberg_trust needs at least one IAM user ARN (arn:aws:iam::<12 digits>:user/...) and one external ID."
  }
}
