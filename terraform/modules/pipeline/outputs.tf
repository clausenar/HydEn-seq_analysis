output "bucket_name" {
  description = "Data bucket. Upload samples to landing/{sample}/; results appear in output/{sample}/."
  value       = aws_s3_bucket.data.bucket
}

output "status_table_name" {
  description = "DynamoDB table with one status row per sample."
  value       = aws_dynamodb_table.status.name
}

output "job_queue_arn" {
  value = aws_batch_job_queue.this.arn
}

output "job_definition_arn" {
  value = aws_batch_job_definition.pipeline.arn
}

output "landing_watcher_function_name" {
  value = aws_lambda_function.landing_watcher.function_name
}

output "alerts_topic_arn" {
  value = aws_sns_topic.alerts.arn
}

output "glue_database_name" {
  description = "Glue database holding the tables over the pipeline's output."
  value       = aws_glue_catalog_database.this.name
}

output "athena_workgroup_name" {
  description = "Athena workgroup to query the Glue tables in (results in athena-results/, 1 GB scan cap)."
  value       = aws_athena_workgroup.this.name
}
