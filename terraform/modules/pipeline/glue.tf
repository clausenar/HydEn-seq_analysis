# Glue Data Catalog table over the pipeline's bedGraphs, queried with Athena.
#
#   Batch job ──► output/{sample}/ (BAMs, plots, bedGraphs ...)
#             └─► tables/bedgraph/{sample}__{strand}.bedgraph ◄── Glue table ◄── Athena
#
# The table points at its own prefix rather than output/, because Athena reads
# every file under a table's location and output/ also holds BAMs and PNGs.
# There are no partitions: a handful of samples is a few tens of MB, and the
# sample and strand come from each row's file name ("$path"), the same way
# stg_bedgraph_signal parses them in Snowflake.

resource "aws_glue_catalog_database" "this" {
  # Glue and Athena names are lowercase with underscores.
  name        = replace(local.name, "-", "_")
  description = "HydEn-seq ${var.stage}: tables over the pipeline's S3 output"
}

resource "aws_glue_catalog_table" "bedgraph_signal_raw" {
  name          = "bedgraph_signal_raw"
  database_name = aws_glue_catalog_database.this.name
  description   = "One row per bedGraph interval; sample and strand are in \"$path\" ({sample}__{forward|reverse}.bedgraph)"
  table_type    = "EXTERNAL_TABLE"

  parameters = {
    EXTERNAL       = "TRUE"
    classification = "csv"
  }

  storage_descriptor {
    location      = "s3://${aws_s3_bucket.data.bucket}/${local.bedgraph_table_prefix}"
    input_format  = "org.apache.hadoop.mapred.TextInputFormat"
    output_format = "org.apache.hadoop.hive.ql.io.HiveIgnoreKeyTextOutputFormat"

    ser_de_info {
      serialization_library = "org.apache.hadoop.hive.serde2.lazy.LazySimpleSerDe"
      parameters            = { "field.delim" = "\t" }
    }

    # Same names and types as BEDGRAPH_SIGNAL_RAW / stg_bedgraph_signal.
    columns {
      name    = "chrom"
      type    = "string"
      comment = "Chromosome, sacCer3 roman-numeral names (chrI ...)"
    }
    columns {
      name    = "start_pos"
      type    = "bigint"
      comment = "0-based interval start"
    }
    columns {
      name    = "end_pos"
      type    = "bigint"
      comment = "0-based interval end (exclusive)"
    }
    columns {
      name    = "score"
      type    = "double"
      comment = "Ribonucleotide pileup count"
    }
  }
}

# ------------------------------------------------------------------- Athena

resource "aws_athena_workgroup" "this" {
  name        = local.name
  description = "HydEn-seq ${var.stage} queries"

  # Deleting a non-prod workgroup also deletes its saved queries.
  force_destroy = !local.is_prod

  configuration {
    # Clients can't override the results location, encryption or the scan cap.
    enforce_workgroup_configuration    = true
    publish_cloudwatch_metrics_enabled = true

    # Athena bills $5 per TB scanned; every table here is MBs, so a query
    # that scans a GB has gone wrong.
    bytes_scanned_cutoff_per_query = 1024 * 1024 * 1024

    engine_version {
      selected_engine_version = "Athena engine version 3"
    }

    result_configuration {
      output_location = "s3://${aws_s3_bucket.data.bucket}/${local.athena_results_prefix}"

      encryption_configuration {
        encryption_option = "SSE_S3"
      }
    }
  }
}

# Saved queries, shown under "Saved queries" in the Athena console.
resource "aws_athena_named_query" "signal_summary" {
  name        = "bedgraph_signal_summary"
  description = "Positions and signal per sample, strand and chromosome (matches the dbt mart in Snowflake)"
  workgroup   = aws_athena_workgroup.this.id
  database    = aws_glue_catalog_database.this.name
  query       = file("${path.module}/athena/bedgraph_signal_summary.sql")
}

resource "aws_athena_named_query" "strand_share" {
  name        = "bedgraph_strand_share"
  description = "Share of each sample's signal on the forward strand"
  workgroup   = aws_athena_workgroup.this.id
  database    = aws_glue_catalog_database.this.name
  query       = file("${path.module}/athena/bedgraph_strand_share.sql")
}
