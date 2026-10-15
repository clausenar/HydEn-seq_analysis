# Fargate compute: no instances to manage, scales to zero when the queue is
# empty, capped at var.max_vcpus. Batch uses its service-linked role.
resource "aws_batch_compute_environment" "fargate" {
  name = "${local.name}-fargate"
  type = "MANAGED"

  compute_resources {
    type               = "FARGATE"
    max_vcpus          = var.max_vcpus
    subnets            = aws_subnet.public[*].id
    security_group_ids = [aws_security_group.batch.id]
  }
}

resource "aws_batch_job_queue" "this" {
  name     = "${local.name}-queue"
  state    = "ENABLED"
  priority = 1

  compute_environment_order {
    order               = 1
    compute_environment = aws_batch_compute_environment.fargate.arn
  }
}

resource "aws_cloudwatch_log_group" "batch_jobs" {
  name              = "/aws/batch/${local.name}"
  retention_in_days = 30
}

resource "aws_batch_job_definition" "pipeline" {
  name                  = "${local.name}-pipeline"
  type                  = "container"
  platform_capabilities = ["FARGATE"]

  # The Lambda passes [input, output, reference, processed_raw] S3 URIs as the
  # command; the image's entrypoint is docker/batch_entrypoint.sh.
  container_properties = jsonencode({
    image = (var.image_digest == null
      ? "${data.aws_ecr_repository.pipeline.repository_url}:${var.image_tag}"
      : "${data.aws_ecr_repository.pipeline.repository_url}@${var.image_digest}"
    )
    jobRoleArn       = aws_iam_role.batch_job.arn
    executionRoleArn = aws_iam_role.batch_execution.arn
    resourceRequirements = [
      { type = "VCPU", value = "8" },
      { type = "MEMORY", value = "16384" },
    ]
    networkConfiguration         = { assignPublicIp = "ENABLED" }
    fargatePlatformConfiguration = { platformVersion = "LATEST" }
    logConfiguration = {
      logDriver = "awslogs"
      options = {
        "awslogs-group"         = aws_cloudwatch_log_group.batch_jobs.name
        "awslogs-region"        = data.aws_region.current.region
        "awslogs-stream-prefix" = "hyden-seq"
      }
    }
  })

  retry_strategy {
    attempts = 1
  }

  timeout {
    attempt_duration_seconds = 12 * 60 * 60
  }
}
