"""Tests for the Glue job, on a local SparkSession with a file-based Iceberg catalog.

Same Spark (3.5) and Iceberg (1.7.1) versions as Glue 5.0. Only the catalog
differs: a Hadoop catalog in a temp folder stands in for the Glue Data Catalog.
Needs Java 17 and network access the first time (Spark downloads the Iceberg jar):

    pip install -r glue/requirements-dev.txt
    pytest glue/tests
"""

import sys
from datetime import datetime, timedelta
from pathlib import Path

import pytest
from pyspark.sql import SparkSession

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bedgraph_to_iceberg as job  # noqa: E402

ICEBERG_PACKAGE = "org.apache.iceberg:iceberg-spark-runtime-3.5_2.12:1.7.1"
DB = "testdb"
NOW = datetime(2026, 10, 21, 12, 0, 0)


@pytest.fixture(scope="session")
def spark(tmp_path_factory):
    warehouse = tmp_path_factory.mktemp("warehouse")
    spark = (SparkSession.builder
             .master("local[2]")
             .appName("test_bedgraph_to_iceberg")
             # Stay on loopback: the machine's network address may be unreachable (VPNs, sandboxes, CI).
             .config("spark.driver.host", "127.0.0.1")
             .config("spark.driver.bindAddress", "127.0.0.1")
             .config("spark.jars.packages", ICEBERG_PACKAGE)
             .config("spark.sql.extensions", "org.apache.iceberg.spark.extensions.IcebergSparkSessionExtensions")
             .config(f"spark.sql.catalog.{job.CATALOG}", "org.apache.iceberg.spark.SparkCatalog")
             .config(f"spark.sql.catalog.{job.CATALOG}.type", "hadoop")
             .config(f"spark.sql.catalog.{job.CATALOG}.warehouse", str(warehouse))
             .config("spark.sql.session.timeZone", "UTC")
             .config("spark.sql.shuffle.partitions", "2")
             .config("spark.ui.enabled", "false")
             .getOrCreate())
    spark.sql(f"CREATE NAMESPACE IF NOT EXISTS {job.CATALOG}.{DB}")
    yield spark
    spark.stop()


@pytest.fixture(autouse=True)
def empty_tables(spark):
    for table in ("bedgraph_signal", "strand_balance"):
        spark.sql(f"DROP TABLE IF EXISTS {job.CATALOG}.{DB}.{table} PURGE")


def write_bedgraph(folder: Path, name: str, rows: list[tuple]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / name).write_text("".join("\t".join(str(v) for v in row) + "\n" for row in rows))


@pytest.fixture
def two_samples(tmp_path):
    """Sample A (with dots in its name, like the Kunkel sample) and sample B."""
    folder = tmp_path / "bedgraph"
    write_bedgraph(folder, "Kunkel.1b.1__forward.bedgraph", [("chrI", 10, 11, 3), ("chrI", 20, 21, 1), ("chrII", 5, 6, 2)])
    write_bedgraph(folder, "Kunkel.1b.1__reverse.bedgraph", [("chrI", 12, 13, 4), ("chrII", 7, 8, 6)])
    write_bedgraph(folder, "SRR1__forward.bedgraph", [("chrI", 1, 2, 5)])
    write_bedgraph(folder, "SRR1__reverse.bedgraph", [("chrI", 3, 4, 5)])
    return f"{folder}/"


def table(spark, name):
    return spark.table(f"{job.CATALOG}.{DB}.{name}")


# ------------------------------------------------------------ transforms

def test_sample_and_strand_come_from_the_file_name(spark, two_samples):
    raw = job.read_bedgraphs(spark, two_samples)
    got = {(r.sample_name, r.strand) for r in raw.select("sample_name", "strand").distinct().collect()}
    assert got == {("Kunkel.1b.1", "forward"), ("Kunkel.1b.1", "reverse"), ("SRR1", "forward"), ("SRR1", "reverse")}


def test_one_sample_reads_only_its_own_files(spark, two_samples):
    raw = job.read_bedgraphs(spark, two_samples, "SRR1")
    assert {r.sample_name for r in raw.collect()} == {"SRR1"}
    assert raw.count() == 2


def test_silver_has_typed_columns(spark, two_samples):
    silver = job.to_silver(job.read_bedgraphs(spark, two_samples), NOW)
    assert [(f.name, f.dataType.simpleString()) for f in silver.schema.fields] == [
        ("sample_name", "string"), ("strand", "string"), ("chrom", "string"),
        ("start_pos", "bigint"), ("end_pos", "bigint"), ("score", "double"),
        ("source_file", "string"), ("processed_at", "timestamp"),
    ]


def test_gold_sums_each_strand_and_the_forward_share(spark, two_samples):
    silver = job.to_silver(job.read_bedgraphs(spark, two_samples, "Kunkel.1b.1"), NOW)
    gold = {r.chrom: r for r in job.to_gold(silver).collect()}
    assert (gold["chrI"].n_positions, gold["chrI"].forward_score, gold["chrI"].reverse_score) == (3, 4.0, 4.0)
    assert gold["chrI"].forward_share == 0.5
    assert gold["chrII"].forward_share == 0.25  # 2 / (2 + 6)


# ------------------------------------------------------------ bad input

def test_a_malformed_line_fails_the_job(spark, tmp_path):
    folder = tmp_path / "bad"
    write_bedgraph(folder, "S__forward.bedgraph", [("chrI", 1, 2, 3), ("chrI", "x", 2, 3)])
    with pytest.raises(Exception, match="(?i)malformed|FAILFAST"):
        job.run(spark, f"{folder}/", DB, now=NOW)
    assert not spark.catalog.tableExists(f"{job.CATALOG}.{DB}.bedgraph_signal")


def test_no_input_files_is_a_clear_error(spark, tmp_path):
    with pytest.raises(ValueError, match="No bedGraph files match"):
        job.run(spark, f"{tmp_path}/", DB, "missing", now=NOW)


def test_a_file_with_an_unknown_strand_stops_the_job(spark, tmp_path):
    folder = tmp_path / "strand"
    write_bedgraph(folder, "S__plus.bedgraph", [("chrI", 1, 2, 3)])
    with pytest.raises(ValueError, match="not named"):
        job.run(spark, f"{folder}/", DB, now=NOW)


# ------------------------------------------------------------ Iceberg writes

def test_run_writes_silver_and_gold(spark, two_samples):
    result = job.run(spark, two_samples, DB, now=NOW)
    assert result == {"silver_rows": 7, "gold_rows": 3, "samples": ["Kunkel.1b.1", "SRR1"]}
    assert table(spark, "bedgraph_signal").count() == 7
    assert table(spark, "strand_balance").count() == 3


def test_tables_are_partitioned_by_sample(spark, two_samples):
    job.run(spark, two_samples, DB, now=NOW)
    for name in ("bedgraph_signal", "strand_balance"):
        partitions = spark.sql(f"SELECT partition.sample_name FROM {job.CATALOG}.{DB}.{name}.partitions").collect()
        assert {p.sample_name for p in partitions} == {"Kunkel.1b.1", "SRR1"}


def test_rerunning_a_sample_replaces_it_and_leaves_the_others(spark, two_samples, tmp_path):
    job.run(spark, two_samples, DB, now=NOW)

    # SRR1 is reprocessed with different output; Kunkel must not change.
    folder = Path(two_samples)
    write_bedgraph(folder, "SRR1__forward.bedgraph", [("chrI", 1, 2, 7), ("chrII", 1, 2, 1)])
    job.run(spark, two_samples, DB, "SRR1", now=NOW)
    job.run(spark, two_samples, DB, "SRR1", now=NOW)  # and again: still no duplicates

    silver = table(spark, "bedgraph_signal")
    assert silver.filter("sample_name = 'SRR1'").count() == 3
    assert silver.filter("sample_name = 'Kunkel.1b.1'").count() == 5
    srr1 = {r.chrom: r.forward_score for r in table(spark, "strand_balance").filter("sample_name = 'SRR1'").collect()}
    assert srr1 == {"chrI": 7.0, "chrII": 1.0}


def test_old_snapshots_are_expired_but_the_last_five_kept(spark, two_samples):
    def snapshots(name):
        return spark.sql(f"SELECT count(*) AS n FROM {job.CATALOG}.{DB}.{name}.snapshots").first().n

    for _ in range(7):
        job.run(spark, two_samples, DB, "SRR1", now=datetime.utcnow())
    assert snapshots("bedgraph_signal") == 7  # all younger than a week: none expired

    # A run 8 days later: everything before it is over a week old, so only
    # the last 5 snapshots (including the new one) survive.
    job.run(spark, two_samples, DB, "SRR1", now=datetime.utcnow() + timedelta(days=8))
    assert snapshots("bedgraph_signal") == job.SNAPSHOTS_TO_KEEP
    assert snapshots("strand_balance") == job.SNAPSHOTS_TO_KEEP
