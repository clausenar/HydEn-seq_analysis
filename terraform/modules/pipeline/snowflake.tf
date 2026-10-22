# Read-only access for Snowflake to the Iceberg tables the Glue job writes, so
# Snowflake queries them in place (no copy, no COPY INTO):
#
#   Snowflake catalog integration (CATALOG_SOURCE = GLUE) ──► Glue Data Catalog: hyden_seq_<stage>
#   Snowflake external volume (ALLOW_WRITES = FALSE)      ──► s3://<bucket>/warehouse/ (Iceberg files)
#
# Both Snowflake objects assume this one role. Linking them takes two steps,
# because each side needs an ID the other one generates:
#   1. apply with snowflake_iceberg_trust = null: the role exists but its trust
#      only matches a placeholder external ID, so nobody can assume it;
#   2. create the two Snowflake objects with this role's ARN (snowflake/iceberg_glue.sql),
#      copy the IAM user ARN and the two external IDs from DESC into
#      snowflake_iceberg_trust, and apply again.
# The external IDs are not secrets (AWS: they prevent the confused-deputy
# problem, not authenticate), so they can live in envs/<stage>/main.tf.

locals {
  snowflake_trust_linked = var.snowflake_iceberg_trust != null
}

resource "aws_iam_role" "snowflake_iceberg" {
  name        = "${local.name}-snowflake-iceberg"
  description = "Snowflake: read the Iceberg tables in hyden_seq_${var.stage} (Glue catalog + warehouse/)"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect = "Allow"
      Principal = {
        AWS = local.snowflake_trust_linked ? var.snowflake_iceberg_trust.iam_user_arns : ["arn:aws:iam::${data.aws_caller_identity.current.account_id}:root"]
      }
      Action = "sts:AssumeRole"
      Condition = {
        StringEquals = {
          "sts:ExternalId" = local.snowflake_trust_linked ? var.snowflake_iceberg_trust.external_ids : ["not-linked-to-snowflake-yet"]
        }
      }
    }]
  })
}

resource "aws_iam_role_policy" "snowflake_iceberg" {
  name = "read-iceberg-tables"
  role = aws_iam_role.snowflake_iceberg.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        # The catalog integration reads table metadata locations; it never
        # changes the catalog.
        Sid    = "ReadGlueCatalog"
        Effect = "Allow"
        Action = ["glue:GetCatalog", "glue:GetDatabase", "glue:GetDatabases", "glue:GetTable", "glue:GetTables"]
        Resource = [
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:catalog",
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:database/${aws_glue_catalog_database.this.name}",
          "arn:aws:glue:${data.aws_region.current.region}:${data.aws_caller_identity.current.account_id}:table/${aws_glue_catalog_database.this.name}/*",
        ]
      },
      {
        Sid      = "BucketLocation"
        Effect   = "Allow"
        Action   = ["s3:GetBucketLocation"]
        Resource = [aws_s3_bucket.data.arn]
      },
      {
        Sid       = "ListWarehouse"
        Effect    = "Allow"
        Action    = ["s3:ListBucket"]
        Resource  = [aws_s3_bucket.data.arn]
        Condition = { StringLike = { "s3:prefix" = ["${local.warehouse_prefix}*"] } }
      },
      {
        # Iceberg metadata and data files; read-only (ALLOW_WRITES = FALSE).
        Sid      = "ReadWarehouse"
        Effect   = "Allow"
        Action   = ["s3:GetObject", "s3:GetObjectVersion"]
        Resource = ["${aws_s3_bucket.data.arn}/${local.warehouse_prefix}*"]
      },
    ]
  })
}
