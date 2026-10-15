resource "aws_s3_bucket" "data" {
  bucket_prefix = "${local.name}-data-"
  force_destroy = !local.is_prod
}

resource "aws_s3_bucket_public_access_block" "data" {
  bucket = aws_s3_bucket.data.id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_s3_bucket_server_side_encryption_configuration" "data" {
  bucket = aws_s3_bucket.data.id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

resource "aws_s3_bucket_policy" "https_only" {
  bucket = aws_s3_bucket.data.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Sid       = "DenyInsecureTransport"
      Effect    = "Deny"
      Principal = "*"
      Action    = ["s3:*"]
      Resource  = [aws_s3_bucket.data.arn, "${aws_s3_bucket.data.arn}/*"]
      Condition = { Bool = { "aws:SecureTransport" = "false" } }
    }]
  })

  # S3 rejects some policies while the public access block is being applied.
  depends_on = [aws_s3_bucket_public_access_block.data]
}

resource "aws_s3_bucket_lifecycle_configuration" "data" {
  bucket = aws_s3_bucket.data.id

  # Raw fastqs are only re-read for reprocessing, so move them to cheaper
  # storage once they have been processed for a month.
  rule {
    id     = "ArchiveProcessedRaw"
    status = "Enabled"

    filter {
      prefix = local.processed_raw_prefix
    }

    transition {
      days          = 30
      storage_class = "STANDARD_IA"
    }
  }

  rule {
    id     = "AbortIncompleteUploads"
    status = "Enabled"

    filter {}

    abort_incomplete_multipart_upload {
      days_after_initiation = 7
    }
  }
}

# Exactly-once job submission + job status per sample (see the Lambda docstring).
resource "aws_dynamodb_table" "status" {
  name                        = "${local.name}-sample-status"
  hash_key                    = "sample_id"
  billing_mode                = "PAY_PER_REQUEST"
  deletion_protection_enabled = local.is_prod

  attribute {
    name = "sample_id"
    type = "S"
  }

  point_in_time_recovery {
    enabled = local.is_prod
  }
}
