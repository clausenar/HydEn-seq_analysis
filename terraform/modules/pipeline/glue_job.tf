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
}

resource "aws_s3_object" "glue_script" {
  bucket = aws_s3_bucket.data.id
  key    = local.glue_script_key
  source = "${path.module}/../../../glue/bedgraph_to_iceberg.py"
  etag   = filemd5("${path.module}/../../../glue/bedgraph_to_iceberg.py")
}

# Glue's own output and error logs go to the account-wide /aws-glue/jobs/*
# groups; the job's continuous logs (driver and executors) go here, kept 30 days.
resource "aws_cloudwatch_log_group" "glue_job" {
  name              = "/aws-glue/${local.name}"
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
    "--job-language"                     = "python"
    "--datalake-formats"                 = "iceberg"
    "--enable-metrics"                   = "true"
    "--enable-continuous-cloudwatch-log" = "true"
    "--continuous-log-logGroup"          = aws_cloudwatch_log_group.glue_job.name
    "--input_path"                       = "s3://${aws_s3_bucket.data.bucket}/${local.bedgraph_table_prefix}"
    "--database"                         = aws_glue_catalog_database.this.name
    "--warehouse"                        = "s3://${aws_s3_bucket.data.bucket}/${local.warehouse_prefix}"
    "--sample_name"                      = "ALL"
  }
}

# ------------------------------------------------------------------- IAM

resource "aws_iam_role" "glue_job" {
  name        = "${local.name}-glue-job"
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
        # the first run, and commit each write by updating the table.
        Sid    = "IcebergTablesInOwnDatabase"
        Effect = "Allow"
        Action = ["glue:GetDatabase", "glue:GetTable", "glue:GetTables", "glue:CreateTable", "glue:UpdateTable"]
        Resource = [
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:catalog",
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:database/${aws_glue_catalog_database.this.name}",
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:table/${aws_glue_catalog_database.this.name}/*",
        ]
      },
      {
        Sid    = "WriteLogs"
        Effect = "Allow"
        Action = ["logs:CreateLogGroup", "logs:CreateLogStream", "logs:PutLogEvents"]
        Resource = [
          "${aws_cloudwatch_log_group.glue_job.arn}:*",
          "arn:aws:logs:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:log-group:/aws-glue/jobs/*",
        ]
      },
      {
        Sid       = "JobMetrics"
        Effect    = "Allow"
        Action    = ["cloudwatch:PutMetricData"]
        Resource  = ["*"] # this action has no resource-level permissions
        Condition = { StringEquals = { "cloudwatch:namespace" = "Glue" } }
      },
    ]
  })
}
