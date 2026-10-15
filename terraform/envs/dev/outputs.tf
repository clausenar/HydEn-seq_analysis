output "bucket_name" {
  value = module.pipeline.bucket_name
}

output "status_table_name" {
  value = module.pipeline.status_table_name
}

output "job_queue_arn" {
  value = module.pipeline.job_queue_arn
}

output "landing_watcher_function_name" {
  value = module.pipeline.landing_watcher_function_name
}

output "alerts_topic_arn" {
  value = module.pipeline.alerts_topic_arn
}
