"""Unit tests for landing_watcher_lambda.py, run entirely against moto-mocked
AWS services (no network access, no real AWS account touched).

Run with:
    pip install -r requirements-dev.txt
    pytest lambda/test_landing_watcher_lambda.py -v
"""

import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import MagicMock, patch

import boto3
import pytest
from moto import mock_aws

# Module-level code in landing_watcher_lambda.py reads these at import time.
os.environ.setdefault("AWS_DEFAULT_REGION", "us-east-1")
os.environ.setdefault("AWS_ACCESS_KEY_ID", "testing")
os.environ.setdefault("AWS_SECRET_ACCESS_KEY", "testing")
os.environ.setdefault("STATUS_TABLE_NAME", "hyden-seq-sample-status")
os.environ.setdefault("BATCH_JOB_QUEUE", "hyden-seq-job-queue")
os.environ.setdefault("BATCH_JOB_DEFINITION", "hyden-seq-pipeline-job")

sys.path.insert(0, os.path.dirname(__file__))  # "lambda" isn't importable as a package name
import landing_watcher_lambda as landing_watcher  # noqa: E402

BUCKET = "kunkel-ribo-data-2026"


def s3_event(bucket, key):
    return {"Records": [{"s3": {"bucket": {"name": bucket}, "object": {"key": key}}}]}


@pytest.fixture(autouse=True)
def aws_mocks():
    with mock_aws():
        boto3.client("s3").create_bucket(Bucket=BUCKET)
        boto3.resource("dynamodb").create_table(
            TableName=os.environ["STATUS_TABLE_NAME"],
            KeySchema=[{"AttributeName": "sample_id", "KeyType": "HASH"}],
            AttributeDefinitions=[{"AttributeName": "sample_id", "AttributeType": "S"}],
            BillingMode="PAY_PER_REQUEST",
        )
        # landing_watcher.py creates its S3/DynamoDB clients at import time,
        # before this mock context exists - point them at the mocked ones.
        landing_watcher.S3 = boto3.client("s3")
        landing_watcher.TABLE = boto3.resource("dynamodb").Table(os.environ["STATUS_TABLE_NAME"])
        yield


class TestParseSample:
    def test_valid_end1_fastq(self):
        assert landing_watcher._parse_sample("landing/sample1/sample1_end1.fastq") == "sample1"

    def test_valid_end2_fastq_gz(self):
        assert landing_watcher._parse_sample("landing/sample1/sample1_end2.fastq.gz") == "sample1"

    def test_sample_name_with_dots(self):
        key = "landing/my.sample.v2/my.sample.v2_end1.fastq.gz"
        assert landing_watcher._parse_sample(key) == "my.sample.v2"

    def test_rejects_mismatched_directory_and_filename(self):
        # defensive: a stray/misnamed upload shouldn't be treated as sample "sample2"
        assert landing_watcher._parse_sample("landing/sample1/sample2_end1.fastq") is None

    def test_rejects_keys_outside_landing_prefix(self):
        assert landing_watcher._parse_sample("output/sample1/sample1_end1.fastq") is None

    def test_rejects_non_fastq_files(self):
        assert landing_watcher._parse_sample("landing/sample1/sample1_end1.bam") is None

    def test_rejects_filename_without_end_marker(self):
        assert landing_watcher._parse_sample("landing/sample1/sample1.fastq") is None


class TestBothMatesPresent:
    def test_false_when_only_one_mate_uploaded(self):
        landing_watcher.S3.put_object(Bucket=BUCKET, Key="landing/s1/s1_end1.fastq", Body=b"x")
        assert landing_watcher._both_mates_present(BUCKET, "s1") is False

    def test_true_once_both_mates_uploaded(self):
        landing_watcher.S3.put_object(Bucket=BUCKET, Key="landing/s1/s1_end1.fastq", Body=b"x")
        landing_watcher.S3.put_object(Bucket=BUCKET, Key="landing/s1/s1_end2.fastq", Body=b"x")
        assert landing_watcher._both_mates_present(BUCKET, "s1") is True

    def test_false_when_sample_has_no_uploads_yet(self):
        assert landing_watcher._both_mates_present(BUCKET, "nonexistent") is False


class TestTryClaim:
    def test_first_claim_succeeds(self):
        assert landing_watcher._try_claim("s1") is True

    def test_second_claim_on_same_sample_fails(self):
        assert landing_watcher._try_claim("s1") is True
        assert landing_watcher._try_claim("s1") is False

    def test_claims_on_different_samples_are_independent(self):
        assert landing_watcher._try_claim("s1") is True
        assert landing_watcher._try_claim("s2") is True

    def test_concurrent_claims_on_same_sample_only_one_wins(self):
        # Exercises the actual conditional-update race the design relies on,
        # rather than just two sequential calls.
        with ThreadPoolExecutor(max_workers=8) as pool:
            results = list(pool.map(landing_watcher._try_claim, ["race-sample"] * 8))
        assert results.count(True) == 1
        assert results.count(False) == 7


class TestBatchJobName:
    # Regression coverage for a real production failure: AWS Batch rejected
    # jobName="hyden-seq-Kunkel_Ribo-seq_Pol2MGrnh201.1b.1" with
    # "Job name should match valid pattern" because of the dots - a sample
    # name that's perfectly valid as a DynamoDB key and S3 path segment.
    VALID_BATCH_JOB_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,127}$")

    def test_simple_sample_name_is_unchanged(self):
        assert landing_watcher._batch_job_name("sample1") == "hyden-seq-sample1"

    def test_dotted_sample_name_matches_batch_naming_pattern(self):
        name = landing_watcher._batch_job_name("Kunkel_Ribo-seq_Pol2MGrnh201.1b.1")
        assert self.VALID_BATCH_JOB_NAME.match(name)

    def test_dotted_sample_name_is_sanitized_predictably(self):
        assert landing_watcher._batch_job_name("Kunkel_Ribo-seq_Pol2MGrnh201.1b.1") == \
            "hyden-seq-Kunkel_Ribo-seq_Pol2MGrnh201-1b-1"


class TestHandler:
    def _put_both_mates(self, sample):
        landing_watcher.S3.put_object(Bucket=BUCKET, Key=f"landing/{sample}/{sample}_end1.fastq", Body=b"x")
        landing_watcher.S3.put_object(Bucket=BUCKET, Key=f"landing/{sample}/{sample}_end2.fastq", Body=b"x")

    @patch("landing_watcher_lambda.BATCH")
    def test_submits_job_once_both_mates_present(self, mock_batch):
        mock_batch.submit_job.return_value = {"jobId": "job-123"}
        self._put_both_mates("sampleA")

        landing_watcher.handler(s3_event(BUCKET, "landing/sampleA/sampleA_end2.fastq"), None)

        mock_batch.submit_job.assert_called_once()
        call_kwargs = mock_batch.submit_job.call_args.kwargs
        assert call_kwargs["jobQueue"] == os.environ["BATCH_JOB_QUEUE"]
        assert call_kwargs["jobDefinition"] == os.environ["BATCH_JOB_DEFINITION"]
        command = call_kwargs["containerOverrides"]["command"]
        assert command == [
            f"s3://{BUCKET}/landing/sampleA/",
            f"s3://{BUCKET}/output/sampleA/",
            f"s3://{BUCKET}/reference/",
            f"s3://{BUCKET}/processed_raw/sampleA/",
        ]
        # Batch's submit_job API rejects an empty-string command element outright.
        assert all(part for part in command)

        item = landing_watcher.TABLE.get_item(Key={"sample_id": "sampleA"})["Item"]
        assert item["status"] == "SUBMITTED"
        assert item["job_id"] == "job-123"

    @patch("landing_watcher_lambda.BATCH")
    def test_submits_job_for_real_dotted_sample_name(self, mock_batch):
        # End-to-end regression test for the dotted-sample-name Batch naming bug.
        sample = "Kunkel_Ribo-seq_Pol2MGrnh201.1b.1"
        mock_batch.submit_job.return_value = {"jobId": "job-789"}
        self._put_both_mates(sample)

        landing_watcher.handler(s3_event(BUCKET, f"landing/{sample}/{sample}_end2.fastq"), None)

        mock_batch.submit_job.assert_called_once()
        job_name = mock_batch.submit_job.call_args.kwargs["jobName"]
        assert TestBatchJobName.VALID_BATCH_JOB_NAME.match(job_name)

    @patch("landing_watcher_lambda.BATCH")
    def test_does_not_submit_when_only_one_mate_present(self, mock_batch):
        landing_watcher.S3.put_object(Bucket=BUCKET, Key="landing/sampleB/sampleB_end1.fastq", Body=b"x")

        landing_watcher.handler(s3_event(BUCKET, "landing/sampleB/sampleB_end1.fastq"), None)

        mock_batch.submit_job.assert_not_called()

    @patch("landing_watcher_lambda.BATCH")
    def test_ignores_non_landing_keys(self, mock_batch):
        landing_watcher.handler(s3_event(BUCKET, "output/sampleC/sampleC__forward.bedgraph"), None)
        mock_batch.submit_job.assert_not_called()

    @patch("landing_watcher_lambda.BATCH")
    def test_second_event_for_already_submitted_sample_is_a_noop(self, mock_batch):
        # Simulates the real-world case this whole design exists for: both
        # mates' S3 events each invoke the Lambda and both see "both mates
        # present" - only the first invocation to reach batch.submit_job
        # should actually call it.
        mock_batch.submit_job.return_value = {"jobId": "job-456"}
        self._put_both_mates("sampleD")

        landing_watcher.handler(s3_event(BUCKET, "landing/sampleD/sampleD_end1.fastq"), None)
        landing_watcher.handler(s3_event(BUCKET, "landing/sampleD/sampleD_end2.fastq"), None)

        mock_batch.submit_job.assert_called_once()
