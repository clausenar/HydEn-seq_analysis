# Tests for the pipeline module, the Terraform counterpart of
# infra/tests/test_pipeline_stack.py.
#
# The AWS provider is mocked, so these need no AWS credentials and create
# nothing in any account. Run from terraform/modules/pipeline/:
#
#   terraform init -backend=false && terraform test

mock_provider "aws" {
  # Realistic ARNs/ids so cross-references in policies stay readable.
  mock_data "aws_region" {
    defaults = { region = "us-east-1" }
  }
  mock_data "aws_caller_identity" {
    defaults = { account_id = "123456789012" }
  }
  mock_data "aws_ecr_repository" {
    defaults = {
      arn            = "arn:aws:ecr:us-east-1:123456789012:repository/hyden-seq-pipeline"
      repository_url = "123456789012.dkr.ecr.us-east-1.amazonaws.com/hyden-seq-pipeline"
    }
  }
  mock_resource "aws_s3_bucket" {
    defaults = {
      arn    = "arn:aws:s3:::hyden-seq-test-data"
      id     = "hyden-seq-test-data"
      bucket = "hyden-seq-test-data"
    }
  }
  mock_resource "aws_dynamodb_table" {
    defaults = {
      arn = "arn:aws:dynamodb:us-east-1:123456789012:table/hyden-seq-test-sample-status"
    }
  }
  mock_resource "aws_iam_role" {
    defaults = {
      arn = "arn:aws:iam::123456789012:role/hyden-seq-test"
    }
  }
  mock_resource "aws_batch_compute_environment" {
    defaults = {
      arn = "arn:aws:batch:us-east-1:123456789012:compute-environment/hyden-seq-test-fargate"
    }
  }
  mock_resource "aws_batch_job_queue" {
    defaults = {
      arn = "arn:aws:batch:us-east-1:123456789012:job-queue/hyden-seq-test-queue"
    }
  }
  mock_resource "aws_batch_job_definition" {
    defaults = {
      arn = "arn:aws:batch:us-east-1:123456789012:job-definition/hyden-seq-test-pipeline:1"
    }
  }
  mock_resource "aws_cloudwatch_log_group" {
    defaults = {
      arn = "arn:aws:logs:us-east-1:123456789012:log-group:hyden-seq-test"
    }
  }
  mock_resource "aws_lambda_function" {
    defaults = {
      arn = "arn:aws:lambda:us-east-1:123456789012:function:hyden-seq-test-landing-watcher"
    }
  }
  mock_resource "aws_sns_topic" {
    defaults = {
      arn = "arn:aws:sns:us-east-1:123456789012:hyden-seq-test-alerts"
    }
  }
  mock_resource "aws_cloudwatch_event_rule" {
    defaults = {
      arn = "arn:aws:events:us-east-1:123456789012:rule/hyden-seq-test-failed-job"
    }
  }
}

variables {
  stage = "dev"
}

# ------------------------------------------------------------------ storage

run "bucket_is_private_encrypted_and_https_only" {
  assert {
    condition = alltrue([
      aws_s3_bucket_public_access_block.data.block_public_acls,
      aws_s3_bucket_public_access_block.data.block_public_policy,
      aws_s3_bucket_public_access_block.data.ignore_public_acls,
      aws_s3_bucket_public_access_block.data.restrict_public_buckets,
    ])
    error_message = "All four public access block settings must be on."
  }

  assert {
    condition     = one(one(aws_s3_bucket_server_side_encryption_configuration.data.rule).apply_server_side_encryption_by_default).sse_algorithm == "AES256"
    error_message = "The bucket must be encrypted at rest."
  }

  assert {
    condition = anytrue([
      for s in jsondecode(aws_s3_bucket_policy.https_only.policy).Statement :
      s.Effect == "Deny" && s.Condition.Bool["aws:SecureTransport"] == "false"
    ])
    error_message = "The bucket policy must deny non-HTTPS requests."
  }
}

run "processed_raw_moves_to_cheaper_storage" {
  assert {
    condition = anytrue([
      for r in aws_s3_bucket_lifecycle_configuration.data.rule :
      one(r.filter).prefix == "processed_raw/" &&
      one(r.transition).storage_class == "STANDARD_IA" &&
      one(r.transition).days == 30
    ])
    error_message = "processed_raw/ must move to STANDARD_IA after 30 days."
  }
}

run "status_table_matches_lambda_contract" {
  assert {
    condition     = aws_dynamodb_table.status.hash_key == "sample_id"
    error_message = "The Lambda keys items on sample_id."
  }

  assert {
    condition     = one(aws_dynamodb_table.status.attribute).type == "S"
    error_message = "sample_id must be a string attribute."
  }

  assert {
    condition     = aws_dynamodb_table.status.billing_mode == "PAY_PER_REQUEST"
    error_message = "The table must be on-demand so it costs nothing when idle."
  }
}

run "dev_is_disposable" {
  assert {
    condition     = aws_s3_bucket.data.force_destroy
    error_message = "A dev bucket must be deletable with its contents."
  }

  assert {
    condition     = !aws_dynamodb_table.status.deletion_protection_enabled
    error_message = "A dev table must be deletable."
  }

  assert {
    condition     = !one(aws_dynamodb_table.status.point_in_time_recovery).enabled
    error_message = "Point-in-time recovery is a prod-only cost."
  }
}

run "prod_keeps_data" {
  variables {
    stage = "prod"
  }

  assert {
    condition     = !aws_s3_bucket.data.force_destroy
    error_message = "The prod bucket must refuse deletion while it holds data."
  }

  assert {
    condition     = aws_dynamodb_table.status.deletion_protection_enabled
    error_message = "The prod table must have deletion protection."
  }

  assert {
    condition     = one(aws_dynamodb_table.status.point_in_time_recovery).enabled
    error_message = "The prod table must have point-in-time recovery."
  }
}

# ------------------------------------------------------------------ compute

run "job_definition_runs_on_fargate_with_8_vcpu_16_gb" {
  assert {
    condition     = aws_batch_job_definition.pipeline.platform_capabilities == toset(["FARGATE"])
    error_message = "The job must run on Fargate."
  }

  assert {
    condition = (
      { for r in jsondecode(aws_batch_job_definition.pipeline.container_properties).resourceRequirements : r.type => r.value }
      == { VCPU = "8", MEMORY = "16384" }
    )
    error_message = "The job needs 8 vCPU and 16 GiB."
  }

  assert {
    condition     = jsondecode(aws_batch_job_definition.pipeline.container_properties).image == "123456789012.dkr.ecr.us-east-1.amazonaws.com/hyden-seq-pipeline:latest"
    error_message = "The job must run the image from the existing ECR repository."
  }

  assert {
    condition     = jsondecode(aws_batch_job_definition.pipeline.container_properties).networkConfiguration.assignPublicIp == "ENABLED"
    error_message = "Without a NAT gateway the task needs a public IP to reach ECR."
  }
}

run "network_has_public_subnets_and_s3_endpoint" {
  assert {
    condition     = length(aws_subnet.public) == 2
    error_message = "Two public subnets are expected."
  }

  assert {
    condition     = one(aws_route_table.public.route).gateway_id == aws_internet_gateway.this.id
    error_message = "The default route must go to the internet gateway (no NAT gateway)."
  }

  assert {
    condition     = aws_vpc_endpoint.s3.vpc_endpoint_type == "Gateway"
    error_message = "S3 traffic must use the free gateway endpoint."
  }

  assert {
    condition     = length(aws_security_group.batch.ingress) == 0
    error_message = "Nothing may connect in to the Batch tasks."
  }
}

run "compute_scales_to_zero_with_a_cap" {
  assert {
    condition     = one(aws_batch_compute_environment.fargate.compute_resources).type == "FARGATE"
    error_message = "Fargate compute has no idle instances."
  }

  assert {
    condition     = one(aws_batch_compute_environment.fargate.compute_resources).max_vcpus == 16
    error_message = "Compute must be capped at 16 vCPU."
  }
}

# ---------------------------------------------------------- landing watcher

run "landing_watcher_gets_its_configuration" {
  assert {
    condition = aws_lambda_function.landing_watcher.environment[0].variables == tomap({
      STATUS_TABLE_NAME    = aws_dynamodb_table.status.name
      BATCH_JOB_QUEUE      = aws_batch_job_queue.this.arn
      BATCH_JOB_DEFINITION = aws_batch_job_definition.pipeline.arn
    })
    error_message = "The Lambda must get the table, queue and job definition it reads from os.environ."
  }

  assert {
    condition     = aws_lambda_function.landing_watcher.handler == "landing_watcher_lambda.handler"
    error_message = "Handler must point at lambda/landing_watcher_lambda.py."
  }
}

run "only_landing_uploads_trigger_the_watcher" {
  assert {
    condition     = one(aws_s3_bucket_notification.landing.lambda_function).filter_prefix == "landing/"
    error_message = "Only uploads under landing/ may trigger the watcher, not the pipeline's own output."
  }

  assert {
    condition     = aws_lambda_permission.s3_invoke.source_arn == aws_s3_bucket.data.arn
    error_message = "Only this bucket may invoke the watcher."
  }
}

run "landing_watcher_is_least_privilege" {
  assert {
    condition = toset(flatten([
      for s in jsondecode(aws_iam_role_policy.landing_watcher.policy).Statement : s.Action
      ])) == toset([
      "s3:ListBucket", "s3:GetObject",
      "dynamodb:GetItem", "dynamodb:PutItem", "dynamodb:UpdateItem",
      "batch:SubmitJob",
      "logs:CreateLogStream", "logs:PutLogEvents",
    ])
    error_message = "The watcher must get exactly the actions it calls, nothing more."
  }

  assert {
    condition = alltrue([
      for s in jsondecode(aws_iam_role_policy.landing_watcher.policy).Statement :
      !contains(s.Resource, "*")
    ])
    error_message = "No watcher permission may use a wildcard resource."
  }
}

run "container_can_read_write_its_bucket_only" {
  assert {
    condition = alltrue(flatten([
      for s in jsondecode(aws_iam_role_policy.batch_job.policy).Statement : [
        for r in s.Resource : startswith(r, aws_s3_bucket.data.arn)
      ]
    ]))
    error_message = "The container role may only touch its own data bucket."
  }

  assert {
    condition = contains(flatten([
      for s in jsondecode(aws_iam_role_policy.batch_job.policy).Statement : s.Action
    ]), "s3:DeleteObject")
    error_message = "`aws s3 mv` to processed_raw/ needs DeleteObject."
  }

  assert {
    condition = length(setsubtract(["s3:GetObjectTagging", "s3:PutObjectTagging"], flatten([
      for s in jsondecode(aws_iam_role_policy.batch_job.policy).Statement : s.Action
    ]))) == 0
    error_message = "`aws s3 mv` of a FASTQ over 8 MB copies tags, so it needs Get/PutObjectTagging."
  }
}

# ---------------------------------------------------------------- alerting

run "failed_jobs_are_alerted" {
  assert {
    condition = jsondecode(aws_cloudwatch_event_rule.failed_job.event_pattern) == {
      source      = ["aws.batch"]
      detail-type = ["Batch Job State Change"]
      detail      = { status = ["FAILED"], jobQueue = [aws_batch_job_queue.this.arn] }
    }
    error_message = "Failed jobs on this queue (and only this queue) must be routed to alerts."
  }

  assert {
    condition     = aws_cloudwatch_event_target.failed_job_to_sns.arn == aws_sns_topic.alerts.arn
    error_message = "The failed-job rule must publish to the alert topic."
  }

  assert {
    condition     = jsondecode(aws_sns_topic_policy.alerts.policy).Statement[0].Principal.Service == "events.amazonaws.com"
    error_message = "EventBridge must be allowed to publish to the topic."
  }

  assert {
    condition     = aws_cloudwatch_metric_alarm.landing_watcher_errors.alarm_actions == toset([aws_sns_topic.alerts.arn])
    error_message = "Lambda errors must alarm to the alert topic."
  }
}

run "finished_jobs_update_the_status_table" {
  assert {
    condition = jsondecode(aws_cloudwatch_event_rule.job_finished.event_pattern) == {
      source      = ["aws.batch"]
      detail-type = ["Batch Job State Change"]
      detail      = { status = ["SUCCEEDED", "FAILED"], jobQueue = [aws_batch_job_queue.this.arn] }
    }
    error_message = "Succeeded and failed jobs on this queue must reach the watcher, or the status stays SUBMITTED."
  }

  assert {
    condition     = aws_cloudwatch_event_target.job_finished_to_watcher.arn == aws_lambda_function.landing_watcher.arn
    error_message = "The job-finished rule must invoke the landing watcher."
  }

  assert {
    condition = (
      aws_lambda_permission.events_invoke.principal == "events.amazonaws.com" &&
      aws_lambda_permission.events_invoke.source_arn == aws_cloudwatch_event_rule.job_finished.arn
    )
    error_message = "Only the job-finished rule may invoke the watcher through EventBridge."
  }

  assert {
    condition = contains(flatten([
      for s in jsondecode(aws_iam_role_policy.landing_watcher.policy).Statement : s.Action
    ]), "dynamodb:UpdateItem")
    error_message = "Recording the result needs dynamodb:UpdateItem."
  }
}

run "alert_email_is_optional" {
  assert {
    condition     = length(aws_sns_topic_subscription.email) == 0
    error_message = "No subscription without an email address."
  }
}

run "alert_email_subscribes_when_given" {
  variables {
    alert_email = "alerts@example.com"
  }

  assert {
    condition     = aws_sns_topic_subscription.email[0].endpoint == "alerts@example.com"
    error_message = "The given address must be subscribed."
  }
}

# ------------------------------------------------------------ glue / athena

run "bedgraph_table_reads_only_its_own_prefix" {
  assert {
    condition     = one(aws_glue_catalog_table.bedgraph_signal_raw.storage_descriptor).location == "s3://hyden-seq-test-data/tables/bedgraph/"
    error_message = "The table must point at tables/bedgraph/, not output/ (which also holds BAMs and plots)."
  }

  assert {
    condition = (
      { for e in jsondecode(aws_batch_job_definition.pipeline.container_properties).environment : e.name => e.value }
      ["HYDEN_S3_BEDGRAPH_TABLE_PREFIX"] == one(aws_glue_catalog_table.bedgraph_signal_raw.storage_descriptor).location
    )
    error_message = "The job must publish bedGraphs to the prefix the table reads."
  }

  assert {
    condition = (
      [for c in one(aws_glue_catalog_table.bedgraph_signal_raw.storage_descriptor).columns : "${c.name}:${c.type}"]
      == ["chrom:string", "start_pos:bigint", "end_pos:bigint", "score:double"]
    )
    error_message = "Columns must match the bedGraph layout and the Snowflake table."
  }

  assert {
    condition     = one(one(aws_glue_catalog_table.bedgraph_signal_raw.storage_descriptor).ser_de_info).parameters["field.delim"] == "\t"
    error_message = "bedGraphs are tab-separated."
  }

  assert {
    condition     = aws_glue_catalog_database.this.name == "hyden_seq_dev"
    error_message = "Glue database names use underscores and include the stage."
  }
}

run "athena_workgroup_is_enforced_capped_and_encrypted" {
  assert {
    condition     = one(aws_athena_workgroup.this.configuration).enforce_workgroup_configuration
    error_message = "Clients must not override the workgroup settings."
  }

  assert {
    condition     = one(aws_athena_workgroup.this.configuration).bytes_scanned_cutoff_per_query == 1073741824
    error_message = "Queries must stop after scanning 1 GB."
  }

  assert {
    condition     = one(one(aws_athena_workgroup.this.configuration).result_configuration).output_location == "s3://hyden-seq-test-data/athena-results/"
    error_message = "Results must go to athena-results/ in the stage's own bucket."
  }

  assert {
    condition     = one(one(one(aws_athena_workgroup.this.configuration).result_configuration).encryption_configuration).encryption_option == "SSE_S3"
    error_message = "Query results must be encrypted."
  }

  assert {
    condition = anytrue([
      for r in aws_s3_bucket_lifecycle_configuration.data.rule :
      one(r.filter).prefix == "athena-results/" && one(r.expiration).days == 7
    ])
    error_message = "Athena results must expire after 7 days."
  }
}

run "prod_workgroup_keeps_saved_queries" {
  variables {
    stage = "prod"
  }

  assert {
    condition     = !aws_athena_workgroup.this.force_destroy
    error_message = "Deleting the prod workgroup must not silently delete its saved queries."
  }
}

# --------------------------------------------------------------- variables

run "image_digest_pins_the_job_to_one_image" {
  variables {
    image_digest = "sha256:5e0cf4728bd7dcb68cc11a1a768504c2aac47fa77901d34bdffbaad3ee64f0c6"
  }

  assert {
    condition     = jsondecode(aws_batch_job_definition.pipeline.container_properties).image == "123456789012.dkr.ecr.us-east-1.amazonaws.com/hyden-seq-pipeline@sha256:5e0cf4728bd7dcb68cc11a1a768504c2aac47fa77901d34bdffbaad3ee64f0c6"
    error_message = "With image_digest set, the job must run that exact image, not a tag."
  }
}

run "rejects_an_invalid_image_digest" {
  command = plan

  variables {
    image_digest = "latest"
  }

  expect_failures = [var.image_digest]
}

run "rejects_an_invalid_stage_name" {
  command = plan

  variables {
    stage = "Prod_1"
  }

  expect_failures = [var.stage]
}
