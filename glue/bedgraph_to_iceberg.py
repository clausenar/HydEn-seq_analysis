"""AWS Glue Spark job: the pipeline's bedGraphs -> Iceberg tables in the Glue Data Catalog.

    tables/bedgraph/{sample}__{strand}.bedgraph      text, written by the Batch job
        -> {database}.bedgraph_signal   silver: one typed row per interval,
                                        sample and strand as columns
        -> {database}.strand_balance    gold: per sample and chromosome, forward
                                        and reverse signal and the forward share

Before anything is written, the silver rows must pass the Glue Data Quality
ruleset below (inside Glue) - a failed rule fails the job and nothing lands.

Both tables are partitioned by sample_name and written with overwritePartitions(),
so rerunning a sample replaces its rows instead of adding to them. A run without
a sample (--sample_name ALL) rewrites every sample found in the input; a sample
whose files were deleted keeps its old rows until it is deleted from the tables.

Job arguments (set by Terraform, see terraform/modules/pipeline/glue_job.tf):
    --input_path   s3://<bucket>/tables/bedgraph/
    --database     Glue database, e.g. hyden_seq_dev
    --warehouse    s3://<bucket>/warehouse/   (Iceberg data and metadata files)
    --sample_name  one sample, or ALL

Everything except main() runs on a plain local SparkSession, which is how
tests/test_bedgraph_to_iceberg.py tests it.
"""

import re
import sys
from datetime import datetime, timedelta, timezone
from typing import Callable

from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as F
from pyspark.sql.types import DoubleType, LongType, StringType, StructField, StructType

CATALOG = "glue_catalog"
ALL_SAMPLES = "ALL"
FILE_NAME = r"^(.+)__(forward|reverse)\.bedgraph$"

# Old snapshots keep replaced data files alive; keep a week of history (and at
# least the last 5 snapshots) for time travel, and let Iceberg delete the rest.
SNAPSHOT_MAX_AGE = timedelta(days=7)
SNAPSHOTS_TO_KEEP = 5

BEDGRAPH_SCHEMA = StructType([
    StructField("chrom", StringType(), nullable=False),
    StructField("start_pos", LongType(), nullable=False),
    StructField("end_pos", LongType(), nullable=False),
    StructField("score", DoubleType(), nullable=False),
])

SACCER3_CHROMS = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
                  "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrM"]

# Glue Data Quality (DQDL) rules for the silver rows of the samples in this run.
# Run by glue_data_quality_gate() inside Glue; results go to the Glue console
# (Data quality tab of the job) and to CloudWatch metrics.
SILVER_RULESET = f"""Rules = [
    RowCount > 0,
    IsComplete "chrom",
    IsComplete "start_pos",
    IsComplete "end_pos",
    IsComplete "score",
    ColumnValues "strand" in ["forward", "reverse"],
    DistinctValuesCount "strand" = 2,
    ColumnValues "chrom" in [{", ".join(f'"{c}"' for c in SACCER3_CHROMS)}],
    ColumnValues "start_pos" >= 0,
    ColumnValues "score" > 0,
    CustomSql "select count(*) from primary where end_pos <= start_pos" = 0
]"""

QualityGate = Callable[[DataFrame], None]

SILVER_DDL = """
    sample_name  string    COMMENT 'From the file name {sample}__{strand}.bedgraph',
    strand       string    COMMENT 'forward or reverse',
    chrom        string    COMMENT 'sacCer3 roman-numeral names (chrI ...)',
    start_pos    bigint    COMMENT '0-based interval start',
    end_pos      bigint    COMMENT '0-based interval end (exclusive)',
    score        double    COMMENT 'Ribonucleotide pileup count',
    source_file  string    COMMENT 'bedGraph file the row came from',
    processed_at timestamp COMMENT 'When this job wrote the row'
"""

GOLD_DDL = """
    sample_name   string,
    chrom         string,
    n_positions   bigint    COMMENT 'Intervals with signal, both strands',
    forward_score double,
    reverse_score double,
    forward_share double    COMMENT 'forward_score / (forward_score + reverse_score)',
    processed_at  timestamp
"""


def glue_catalog_conf(warehouse: str) -> dict:
    """Spark settings for an Iceberg catalog backed by the Glue Data Catalog."""
    return {
        # processed_at and the snapshot cut-off are naive UTC times.
        "spark.sql.session.timeZone": "UTC",
        "spark.sql.extensions": "org.apache.iceberg.spark.extensions.IcebergSparkSessionExtensions",
        f"spark.sql.catalog.{CATALOG}": "org.apache.iceberg.spark.SparkCatalog",
        f"spark.sql.catalog.{CATALOG}.catalog-impl": "org.apache.iceberg.aws.glue.GlueCatalog",
        f"spark.sql.catalog.{CATALOG}.io-impl": "org.apache.iceberg.aws.s3.S3FileIO",
        f"spark.sql.catalog.{CATALOG}.warehouse": warehouse,
    }


def read_bedgraphs(spark: SparkSession, input_path: str, sample_name: str = ALL_SAMPLES) -> DataFrame:
    """Read the bedGraph files for one sample (or all), with sample and strand from each file name.

    Fails on any line that doesn't parse (FAILFAST) rather than loading it as nulls.
    """
    pattern = "*" if sample_name == ALL_SAMPLES else sample_name
    path = f"{input_path.rstrip('/')}/{pattern}__*.bedgraph"
    try:
        df = (spark.read
              .option("sep", "\t")
              .option("mode", "FAILFAST")
              .schema(BEDGRAPH_SCHEMA)
              .csv(path))
    except Exception as exc:  # AnalysisException: no file matches the path
        if "PATH_NOT_FOUND" in str(exc) or "Path does not exist" in str(exc):
            raise ValueError(f"No bedGraph files match {path}") from exc
        raise

    file_name = F.col("_metadata.file_name")
    return (df
            .withColumn("sample_name", F.regexp_extract(file_name, FILE_NAME, 1))
            .withColumn("strand", F.regexp_extract(file_name, FILE_NAME, 2))
            .withColumn("source_file", F.col("_metadata.file_path")))


def to_silver(raw: DataFrame, processed_at: datetime) -> DataFrame:
    return raw.select(
        "sample_name", "strand", "chrom", "start_pos", "end_pos", "score", "source_file",
        F.lit(processed_at).cast("timestamp").alias("processed_at"),
    )


def to_gold(silver: DataFrame) -> DataFrame:
    forward = F.sum(F.when(F.col("strand") == "forward", F.col("score")).otherwise(0.0))
    reverse = F.sum(F.when(F.col("strand") == "reverse", F.col("score")).otherwise(0.0))
    return (silver
            .groupBy("sample_name", "chrom")
            .agg(F.count("*").alias("n_positions"),
                 forward.alias("forward_score"),
                 reverse.alias("reverse_score"),
                 F.max("processed_at").alias("processed_at"))
            .withColumn("forward_share",
                        F.col("forward_score") / (F.col("forward_score") + F.col("reverse_score")))
            .select("sample_name", "chrom", "n_positions", "forward_score", "reverse_score",
                    "forward_share", "processed_at"))


def check_silver(silver: DataFrame) -> int:
    """Stop before writing anything if the input is empty or a file name didn't parse."""
    counts = silver.agg(
        F.count("*").alias("rows"),
        F.sum(F.when(F.col("sample_name") == "", 1).otherwise(0)).alias("unparsed"),
    ).first()
    if counts["rows"] == 0:
        raise ValueError("The bedGraph files hold no rows")
    if counts["unparsed"]:
        raise ValueError(f"{counts['unparsed']} rows come from files not named {{sample}}__{{strand}}.bedgraph")
    return counts["rows"]


def write_partitions(df: DataFrame, table: str, ddl: str) -> None:
    """Create the table if needed, then replace only the sample partitions present in df."""
    spark = df.sparkSession
    spark.sql(f"""
        CREATE TABLE IF NOT EXISTS {table} ({ddl})
        USING iceberg
        PARTITIONED BY (sample_name)
        TBLPROPERTIES ('format-version' = '2')
    """)
    df.writeTo(table).overwritePartitions()


def expire_snapshots(spark: SparkSession, table: str, now: datetime) -> None:
    older_than = (now - SNAPSHOT_MAX_AGE).strftime("%Y-%m-%d %H:%M:%S")
    spark.sql(f"""
        CALL {CATALOG}.system.expire_snapshots(
            table => '{table.split('.', 1)[1]}',
            older_than => TIMESTAMP '{older_than}',
            retain_last => {SNAPSHOTS_TO_KEEP})
    """)


def no_quality_gate(silver: DataFrame) -> None:
    """Outside Glue (tests, local runs) there is no Glue Data Quality to call."""


def data_quality_context(sample_name: str) -> str:
    """Name for the Data Quality results: Glue allows only letters, digits, - and _."""
    return "bedgraph_signal_" + re.sub(r"[^A-Za-z0-9_-]", "_", sample_name)


def glue_data_quality_gate(spark: SparkSession, context: str) -> QualityGate:
    """Evaluate SILVER_RULESET with Glue Data Quality; raise if any rule fails."""
    from awsglue.context import GlueContext  # only exist inside Glue
    from awsglue.dynamicframe import DynamicFrame
    from awsgluedq.transforms import EvaluateDataQuality

    glue_context = GlueContext(spark.sparkContext)

    def gate(silver: DataFrame) -> None:
        outcomes = EvaluateDataQuality.apply(
            frame=DynamicFrame.fromDF(silver, glue_context, "silver"),
            ruleset=SILVER_RULESET,
            publishing_options={
                "dataQualityEvaluationContext": context,
                "enableDataQualityCloudWatchMetrics": True,
                "enableDataQualityResultsPublishing": True,
            },
        ).toDF()
        failed = outcomes.filter(F.col("Outcome") != "Passed").select("Rule", "FailureReason").collect()
        for rule in outcomes.select("Rule", "Outcome").collect():
            print(f"data quality: {rule['Outcome']}: {rule['Rule']}")
        if failed:
            reasons = "; ".join(f"{r['Rule']} ({r['FailureReason']})" for r in failed)
            raise ValueError(f"{len(failed)} data quality rule(s) failed, nothing written: {reasons}")

    return gate


def run(spark: SparkSession, input_path: str, database: str, sample_name: str = ALL_SAMPLES,
        now: datetime | None = None, quality_gate: QualityGate = no_quality_gate) -> dict:
    """The whole job. Returns row counts for the log (and the tests)."""
    now = now or datetime.now(timezone.utc).replace(tzinfo=None)
    silver = to_silver(read_bedgraphs(spark, input_path, sample_name), now).cache()
    rows = check_silver(silver)
    quality_gate(silver)  # before any write: a failed rule leaves the tables as they were
    gold = to_gold(silver)

    silver_table = f"{CATALOG}.{database}.bedgraph_signal"
    gold_table = f"{CATALOG}.{database}.strand_balance"
    write_partitions(silver, silver_table, SILVER_DDL)
    write_partitions(gold, gold_table, GOLD_DDL)
    for table in (silver_table, gold_table):
        expire_snapshots(spark, table, now)

    samples = sorted(r["sample_name"] for r in silver.select("sample_name").distinct().collect())
    return {"silver_rows": rows, "gold_rows": gold.count(), "samples": samples}


def main() -> None:
    from awsglue.utils import getResolvedOptions  # only exists inside Glue

    args = getResolvedOptions(sys.argv, ["input_path", "database", "warehouse", "sample_name"])
    builder = SparkSession.builder.appName("bedgraph_to_iceberg")
    for key, value in glue_catalog_conf(args["warehouse"]).items():
        builder = builder.config(key, value)
    spark = builder.getOrCreate()

    gate = glue_data_quality_gate(spark, context=data_quality_context(args["sample_name"]))
    result = run(spark, args["input_path"], args["database"], args["sample_name"], quality_gate=gate)
    print(f"bedgraph_to_iceberg: {result}")


if __name__ == "__main__":
    main()
