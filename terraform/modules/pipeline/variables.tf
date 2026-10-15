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
  description = "Image tag the Batch job definition runs."
  type        = string
  default     = "latest"
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
