resource "aws_sns_topic" "alerts" {
  name         = "${local.name}-alerts"
  display_name = "HydEn-seq ${var.stage} alerts"
}

resource "aws_sns_topic_subscription" "email" {
  # "" counts as unset: CI passes an empty string when the secret is missing.
  count = var.alert_email == null || var.alert_email == "" ? 0 : 1

  topic_arn = aws_sns_topic.alerts.arn
  protocol  = "email"
  endpoint  = var.alert_email
}

# Batch has no built-in "failed jobs" metric, so failed jobs on this queue are
# routed to the alert topic through EventBridge.
resource "aws_cloudwatch_event_rule" "failed_job" {
  name        = "${local.name}-failed-job"
  description = "Alert when a HydEn-seq Batch job fails"
  event_pattern = jsonencode({
    source      = ["aws.batch"]
    detail-type = ["Batch Job State Change"]
    detail = {
      status   = ["FAILED"]
      jobQueue = [aws_batch_job_queue.this.arn]
    }
  })
}

resource "aws_cloudwatch_event_target" "failed_job_to_sns" {
  rule = aws_cloudwatch_event_rule.failed_job.name
  arn  = aws_sns_topic.alerts.arn
}

# EventBridge may publish to the topic, but only for this rule.
resource "aws_sns_topic_policy" "alerts" {
  arn = aws_sns_topic.alerts.arn
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Sid       = "AllowFailedJobRule"
      Effect    = "Allow"
      Principal = { Service = "events.amazonaws.com" }
      Action    = "sns:Publish"
      Resource  = aws_sns_topic.alerts.arn
      Condition = { ArnEquals = { "aws:SourceArn" = aws_cloudwatch_event_rule.failed_job.arn } }
    }]
  })
}

resource "aws_cloudwatch_metric_alarm" "landing_watcher_errors" {
  alarm_name          = "${local.name}-landing-watcher-errors"
  alarm_description   = "The landing watcher Lambda raised an error"
  namespace           = "AWS/Lambda"
  metric_name         = "Errors"
  dimensions          = { FunctionName = aws_lambda_function.landing_watcher.function_name }
  statistic           = "Sum"
  period              = 300
  evaluation_periods  = 1
  threshold           = 1
  comparison_operator = "GreaterThanOrEqualToThreshold"
  treat_missing_data  = "notBreaching"
  alarm_actions       = [aws_sns_topic.alerts.arn]
}
