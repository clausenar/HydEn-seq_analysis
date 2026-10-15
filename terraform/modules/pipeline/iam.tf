# Least-privilege roles, written out in full. Policies are built with
# jsonencode() rather than aws_iam_policy_document so the tests in tests/ can
# decode and check them without AWS access.

locals {
  ecs_tasks_trust = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "ecs-tasks.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

# ------------------------------------------------------------ Batch execution
# Used by Fargate itself (not the container) to pull the image and write logs.
resource "aws_iam_role" "batch_execution" {
  name               = "${local.name}-batch-execution"
  description        = "Fargate: pull the HydEn-seq image and write its logs"
  assume_role_policy = local.ecs_tasks_trust
}

resource "aws_iam_role_policy" "batch_execution" {
  name = "pull-image-write-logs"
  role = aws_iam_role.batch_execution.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid      = "EcrAuth"
        Effect   = "Allow"
        Action   = ["ecr:GetAuthorizationToken"]
        Resource = ["*"] # this action has no resource-level permissions
      },
      {
        Sid      = "PullPipelineImage"
        Effect   = "Allow"
        Action   = ["ecr:BatchCheckLayerAvailability", "ecr:BatchGetImage", "ecr:GetDownloadUrlForLayer"]
        Resource = [data.aws_ecr_repository.pipeline.arn]
      },
      {
        Sid      = "WriteJobLogs"
        Effect   = "Allow"
        Action   = ["logs:CreateLogStream", "logs:PutLogEvents"]
        Resource = ["${aws_cloudwatch_log_group.batch_jobs.arn}:*"]
      },
    ]
  })
}

# ------------------------------------------------------------------ Batch job
# The container's own role: read and write the data bucket only
# (write + delete are needed for `aws s3 mv` to processed_raw/).
resource "aws_iam_role" "batch_job" {
  name               = "${local.name}-batch-job"
  description        = "HydEn-seq pipeline container: read/write its own data bucket only"
  assume_role_policy = local.ecs_tasks_trust
}

resource "aws_iam_role_policy" "batch_job" {
  name = "data-bucket-read-write"
  role = aws_iam_role.batch_job.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid      = "ListBucket"
        Effect   = "Allow"
        Action   = ["s3:ListBucket", "s3:GetBucketLocation"]
        Resource = [aws_s3_bucket.data.arn]
      },
      {
        # The tagging actions are for `aws s3 mv`: a file over the CLI's 8 MB
        # multipart threshold is copied part by part, and the CLI reads the
        # source's tags and writes them on the copy.
        Sid    = "ReadWriteObjects"
        Effect = "Allow"
        Action = [
          "s3:GetObject", "s3:PutObject", "s3:DeleteObject", "s3:AbortMultipartUpload",
          "s3:GetObjectTagging", "s3:PutObjectTagging",
        ]
        Resource = ["${aws_s3_bucket.data.arn}/*"]
      },
    ]
  })
}

# ------------------------------------------------------------ landing watcher
resource "aws_iam_role" "landing_watcher" {
  name        = "${local.name}-landing-watcher"
  description = "HydEn-seq landing watcher Lambda"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "lambda.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

# Matches what the Lambda actually calls: s3:ListBucket (both mates present?),
# three DynamoDB item actions (the claim) and batch:SubmitJob.
resource "aws_iam_role_policy" "landing_watcher" {
  name = "landing-watcher"
  role = aws_iam_role.landing_watcher.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid      = "ReadDataBucket"
        Effect   = "Allow"
        Action   = ["s3:ListBucket", "s3:GetObject"]
        Resource = [aws_s3_bucket.data.arn, "${aws_s3_bucket.data.arn}/*"]
      },
      {
        Sid      = "ClaimSample"
        Effect   = "Allow"
        Action   = ["dynamodb:GetItem", "dynamodb:PutItem", "dynamodb:UpdateItem"]
        Resource = [aws_dynamodb_table.status.arn]
      },
      {
        Sid      = "SubmitPipelineJob"
        Effect   = "Allow"
        Action   = ["batch:SubmitJob"]
        Resource = [aws_batch_job_queue.this.arn, aws_batch_job_definition.pipeline.arn]
      },
      {
        Sid      = "WriteLogs"
        Effect   = "Allow"
        Action   = ["logs:CreateLogStream", "logs:PutLogEvents"]
        Resource = ["${aws_cloudwatch_log_group.landing_watcher.arn}:*"]
      },
    ]
  })
}
