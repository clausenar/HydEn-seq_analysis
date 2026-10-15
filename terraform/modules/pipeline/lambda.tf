data "archive_file" "landing_watcher" {
  type        = "zip"
  source_file = "${path.module}/../../../lambda/landing_watcher_lambda.py"
  output_path = "${path.module}/.build/landing_watcher.zip"
}

resource "aws_cloudwatch_log_group" "landing_watcher" {
  name              = "/aws/lambda/${local.name}-landing-watcher"
  retention_in_days = 30
}

resource "aws_lambda_function" "landing_watcher" {
  function_name    = "${local.name}-landing-watcher"
  description      = "Submits one Batch job per sample once both fastq mates have landed"
  role             = aws_iam_role.landing_watcher.arn
  runtime          = "python3.12"
  handler          = "landing_watcher_lambda.handler"
  filename         = data.archive_file.landing_watcher.output_path
  source_code_hash = data.archive_file.landing_watcher.output_base64sha256
  timeout          = 30
  memory_size      = 256

  logging_config {
    log_format = "Text"
    log_group  = aws_cloudwatch_log_group.landing_watcher.name
  }

  environment {
    variables = {
      STATUS_TABLE_NAME    = aws_dynamodb_table.status.name
      BATCH_JOB_QUEUE      = aws_batch_job_queue.this.arn
      BATCH_JOB_DEFINITION = aws_batch_job_definition.pipeline.arn
    }
  }

  # The role's inline policy must exist before the first invocation.
  depends_on = [aws_iam_role_policy.landing_watcher]
}

resource "aws_lambda_permission" "s3_invoke" {
  statement_id   = "AllowS3Invoke"
  action         = "lambda:InvokeFunction"
  function_name  = aws_lambda_function.landing_watcher.function_name
  principal      = "s3.amazonaws.com"
  source_arn     = aws_s3_bucket.data.arn
  source_account = data.aws_caller_identity.current.account_id
}

# Only uploads under landing/ trigger the watcher; output/ and processed_raw/
# writes by the pipeline itself must not.
resource "aws_s3_bucket_notification" "landing" {
  bucket = aws_s3_bucket.data.id

  lambda_function {
    lambda_function_arn = aws_lambda_function.landing_watcher.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = local.landing_prefix
  }

  depends_on = [aws_lambda_permission.s3_invoke]
}
