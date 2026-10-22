# Glue Spark job: tables/bedgraph/ text files -> Iceberg tables in the stage's
# Glue database (see glue/bedgraph_to_iceberg.py for what it writes and why).
#
#   tables/bedgraph/*.bedgraph ──► Glue 5.0 job (Spark 3.5, 2 × G.1X) ──► warehouse/ (Iceberg files)
#                                                                    └─► hyden_seq_<stage>.bedgraph_signal
#                                                                        hyden_seq_<stage>.strand_balance
#
# Phase 2 runs it by hand (aws glue start-job-run); phase 3 will start it when
# a Batch job succeeds. The Iceberg tables are created by the job itself
# (CREATE TABLE IF NOT EXISTS), so their schema lives next to the code that
# writes them; Terraform owns the database they go in.

locals {
  glue_script_key  = "glue-scripts/bedgraph_to_iceberg.py"
  warehouse_prefix = "warehouse/"

  # With a security configuration attached, Glue logs to groups named after it
  # and the job role (not the shared /aws-glue/jobs/output and /error).
  glue_job_role_name = "${local.name}-glue-job"
  glue_log_prefix    = "/aws-glue/jobs/${local.name}-role/${local.glue_job_role_name}"
}

resource "aws_s3_object" "glue_script" {
  bucket = aws_s3_bucket.data.id
  key    = local.glue_script_key
  source = "${path.module}/../../../glue/bedgraph_to_iceberg.py"
  etag   = filemd5("${path.module}/../../../glue/bedgraph_to_iceberg.py")
}

# Driver and executor logs (one stream per run and per executor). Glue would
# create these itself with no expiry; created here, they are kept 30 days.
resource "aws_cloudwatch_log_group" "glue_job" {
  for_each          = toset(["output", "error"])
  name              = "${local.glue_log_prefix}/${each.key}"
  retention_in_days = 30
}

resource "aws_glue_security_configuration" "this" {
  name = local.name

  encryption_configuration {
    s3_encryption {
      s3_encryption_mode = "SSE-S3"
    }

    # Both need a customer-managed KMS key (see .checkov.yaml).
    cloudwatch_encryption {
      cloudwatch_encryption_mode = "DISABLED"
    }

    job_bookmarks_encryption {
      job_bookmarks_encryption_mode = "DISABLED"
    }
  }
}

resource "aws_glue_job" "bedgraph_to_iceberg" {
  name                   = "${local.name}-bedgraph-to-iceberg"
  description            = "bedGraphs in tables/bedgraph/ -> Iceberg bedgraph_signal (silver) and strand_balance (gold)"
  role_arn               = aws_iam_role.glue_job.arn
  glue_version           = "5.0"
  worker_type            = "G.1X"
  number_of_workers      = 2
  timeout                = 30 # minutes; a run takes a few
  max_retries            = 0
  security_configuration = aws_glue_security_configuration.this.name

  command {
    name            = "glueetl"
    script_location = "s3://${aws_s3_bucket.data.bucket}/${aws_s3_object.glue_script.key}"
    python_version  = "3"
  }

  # Runs for different samples write different partitions, which Iceberg
  # commits safely side by side (phase 3 starts one run per finished sample).
  execution_property {
    max_concurrent_runs = 3
  }

  default_arguments = {
    "--job-language"     = "python"
    "--datalake-formats" = "iceberg"
    "--enable-metrics"   = "true"
    "--input_path"       = "s3://${aws_s3_bucket.data.bucket}/${local.bedgraph_table_prefix}"
    "--database"         = aws_glue_catalog_database.this.name
    "--warehouse"        = "s3://${aws_s3_bucket.data.bucket}/${local.warehouse_prefix}"
    "--sample_name"      = "ALL"
  }

  # The log groups must exist before the first run, or Glue creates them unmanaged.
  depends_on = [aws_cloudwatch_log_group.glue_job]
}

# ------------------------------------------------------------------- IAM

resource "aws_iam_role" "glue_job" {
  name        = local.glue_job_role_name
  description = "HydEn-seq Glue job: read bedGraphs, write Iceberg tables in its own database"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "glue.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy" "glue_job" {
  name = "bedgraph-to-iceberg"
  role = aws_iam_role.glue_job.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid      = "ListOwnPrefixes"
        Effect   = "Allow"
        Action   = ["s3:ListBucket"]
        Resource = [aws_s3_bucket.data.arn]
        Condition = {
          StringLike = {
            "s3:prefix" = ["${local.bedgraph_table_prefix}*", "${local.warehouse_prefix}*", "glue-scripts/*"]
          }
        }
      },
      {
        Sid    = "ReadInputAndScript"
        Effect = "Allow"
        Action = ["s3:GetObject"]
        Resource = [
          "${aws_s3_bucket.data.arn}/${local.bedgraph_table_prefix}*",
          "${aws_s3_bucket.data.arn}/${local.glue_script_key}",
        ]
      },
      {
        # Iceberg data and metadata files; delete is for expired snapshots.
        Sid      = "ReadWriteWarehouse"
        Effect   = "Allow"
        Action   = ["s3:GetObject", "s3:PutObject", "s3:DeleteObject", "s3:AbortMultipartUpload"]
        Resource = ["${aws_s3_bucket.data.arn}/${local.warehouse_prefix}*"]
      },
      {
        # Iceberg's GlueCatalog: look up the database, create the tables on
        # the first run, and commit each write by updating the table. Glue 5.0
        # also looks up the catalog itself at start-up (GetCatalog).
        Sid    = "IcebergTablesInOwnDatabase"
        Effect = "Allow"
        Action = ["glue:GetCatalog", "glue:GetDatabase", "glue:GetTable", "glue:GetTables", "glue:CreateTable", "glue:UpdateTable"]
        Resource = [
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:catalog",
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:database/${aws_glue_catalog_database.this.name}",
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:table/${aws_glue_catalog_database.this.name}/*",
        ]
      },
      {
        # Only into the two groups above: no CreateLogGroup, so Glue can't
        # fall back to creating unmanaged ones.
        Sid      = "WriteLogs"
        Effect   = "Allow"
        Action   = ["logs:CreateLogStream", "logs:PutLogEvents"]
        Resource = [for g in aws_cloudwatch_log_group.glue_job : "${g.arn}:*"]
      },
      {
        Sid       = "JobMetrics"
        Effect    = "Allow"
        Action    = ["cloudwatch:PutMetricData"]
        Resource  = ["*"] # this action has no resource-level permissions
        Condition = { StringEquals = { "cloudwatch:namespace" = ["Glue", "Glue Data Quality"] } }
      },
      {
        # The job's Data Quality results, shown on the job's Data quality tab.
        # Glue checks this against a ruleset named after each run's evaluation
        # context (data_quality_context() in the script), not the job.
        Sid      = "PublishDataQualityResults"
        Effect   = "Allow"
        Action   = ["glue:PublishDataQuality"]
        Resource = ["arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:dataQualityRuleset/bedgraph_signal_*"]
      },
    ]
  })
}
