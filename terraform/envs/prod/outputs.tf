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

output "glue_database_name" {
  value = module.pipeline.glue_database_name
}

output "athena_workgroup_name" {
  value = module.pipeline.athena_workgroup_name
}

output "glue_job_name" {
  value = module.pipeline.glue_job_name
}
