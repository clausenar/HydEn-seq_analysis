# Runs the Glue job for each sample the pipeline finishes, so a FASTQ upload
# ends as rows in the Iceberg tables with no manual step:
#
#   Batch job SUCCEEDED ──► EventBridge ──► Step Functions (one execution per sample)
#                                              RunGlueJob (startJobRun.sync, --sample_name)
#                                                 │ retries while 3 runs are already going
#                                                 ├─ ok ─────► Succeeded
#                                                 └─ failed ─► SNS alert naming the sample ─► Failed
#
# The Data Quality rules run inside the Glue job before it writes, so a failed
# rule fails the job and arrives here as an ordinary job failure.

resource "aws_cloudwatch_log_group" "orchestration" {
  name              = "/aws/vendedlogs/states/${local.name}-to-iceberg"
  retention_in_days = 30
}

resource "aws_sfn_state_machine" "to_iceberg" {
  name     = "${local.name}-to-iceberg"
  role_arn = aws_iam_role.orchestration.arn

  definition = jsonencode({
    Comment = "Write one finished sample to the Iceberg tables with the Glue job"
    StartAt = "RunGlueJob"
    States = {
      RunGlueJob = {
        Type     = "Task"
        Resource = "arn:aws:states:::glue:startJobRun.sync"
        Parameters = {
          JobName = aws_glue_job.bedgraph_to_iceberg.name
          Arguments = {
            "--sample_name.$" = "$.sample_id"
          }
        }
        ResultSelector = {
          "JobRunId.$"      = "$.Id"
          "ExecutionTime.$" = "$.ExecutionTime"
        }
        ResultPath = "$.glue"
        Retry = [{
          # max_concurrent_runs on the job; wait for a slot rather than fail.
          ErrorEquals     = ["Glue.ConcurrentRunsExceededException"]
          IntervalSeconds = 60
          BackoffRate     = 2
          MaxAttempts     = 5
        }]
        Catch = [{
          ErrorEquals = ["States.ALL"]
          ResultPath  = "$.error"
          Next        = "AlertFailure"
        }]
        Next = "Succeeded"
      }
      AlertFailure = {
        Type     = "Task"
        Resource = "arn:aws:states:::sns:publish"
        Parameters = {
          TopicArn    = aws_sns_topic.alerts.arn
          Subject     = "HydEn-seq ${var.stage}: Iceberg load failed"
          "Message.$" = "States.Format('Sample {} was processed by Batch (job {}) but the Glue job did not write it to the Iceberg tables. Error: {}. Cause: {}', $.sample_id, $.batch_job_id, $.error.Error, $.error.Cause)"
        }
        ResultPath = null
        Next       = "Failed"
      }
      Failed = {
        Type  = "Fail"
        Error = "GlueJobFailed"
        Cause = "The Glue job failed; an alert was sent."
      }
      Succeeded = {
        Type = "Succeed"
      }
    }
  })

  logging_configuration {
    log_destination        = "${aws_cloudwatch_log_group.orchestration.arn}:*"
    include_execution_data = true
    level                  = "ERROR"
  }
}

# Same event as the landing watcher's job_finished rule, but only successes,
# and only jobs the Lambda submitted (they carry the sample_id parameter).
resource "aws_cloudwatch_event_rule" "job_succeeded" {
  name        = "${local.name}-job-succeeded-to-iceberg"
  description = "Load each sample the pipeline finishes into the Iceberg tables"
  event_pattern = jsonencode({
    source      = ["aws.batch"]
    detail-type = ["Batch Job State Change"]
    detail = {
      status     = ["SUCCEEDED"]
      jobQueue   = [aws_batch_job_queue.this.arn]
      parameters = { sample_id = [{ exists = true }] }
    }
  })
}

resource "aws_cloudwatch_event_target" "job_succeeded_to_state_machine" {
  rule     = aws_cloudwatch_event_rule.job_succeeded.name
  arn      = aws_sfn_state_machine.to_iceberg.arn
  role_arn = aws_iam_role.start_orchestration.arn

  input_transformer {
    input_paths = {
      sample_id    = "$.detail.parameters.sample_id"
      batch_job_id = "$.detail.jobId"
    }
    input_template = <<-EOT
      {"sample_id": <sample_id>, "batch_job_id": <batch_job_id>}
    EOT
  }
}

# ------------------------------------------------------------------- IAM

resource "aws_iam_role" "start_orchestration" {
  name        = "${local.name}-start-to-iceberg"
  description = "EventBridge: start the to-Iceberg state machine"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "events.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy" "start_orchestration" {
  name = "start-to-iceberg"
  role = aws_iam_role.start_orchestration.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Sid      = "StartExecution"
      Effect   = "Allow"
      Action   = ["states:StartExecution"]
      Resource = [aws_sfn_state_machine.to_iceberg.arn]
    }]
  })
}

resource "aws_iam_role" "orchestration" {
  name        = "${local.name}-to-iceberg"
  description = "Step Functions: run the Glue job and alert on failure"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "states.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy" "orchestration" {
  name = "to-iceberg"
  role = aws_iam_role.orchestration.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        # startJobRun.sync starts the run, polls it, and stops it if the
        # execution is stopped.
        Sid      = "RunGlueJob"
        Effect   = "Allow"
        Action   = ["glue:StartJobRun", "glue:GetJobRun", "glue:GetJobRuns", "glue:BatchStopJobRun"]
        Resource = [aws_glue_job.bedgraph_to_iceberg.arn]
      },
      {
        Sid      = "Alert"
        Effect   = "Allow"
        Action   = ["sns:Publish"]
        Resource = [aws_sns_topic.alerts.arn]
      },
      {
        # Step Functions log delivery; these actions have no resource-level
        # permissions.
        Sid    = "DeliverLogs"
        Effect = "Allow"
        Action = [
          "logs:CreateLogDelivery", "logs:GetLogDelivery", "logs:UpdateLogDelivery",
          "logs:DeleteLogDelivery", "logs:ListLogDeliveries", "logs:PutResourcePolicy",
          "logs:DescribeResourcePolicies", "logs:DescribeLogGroups",
        ]
        Resource = ["*"]
      },
    ]
  })
}
