"""hyden-seq-landing-watcher: submits a Batch job once both mates land.

Triggered by S3 ObjectCreated events under landing/{sample}/. The two mates
of a pair (_end1/_end2) upload as separate S3 events, often within
milliseconds of each other, so both invocations may see "both mates present"
at the same time. Exactly-once submission is enforced with a conditional
DynamoDB update (table STATUS_TABLE_NAME, partition key sample_id): only the
invocation that successfully flips status to SUBMITTED calls batch.submit_job.
The other invocation's condition check fails and it exits quietly. This
replaces an earlier S3 conditional-write lock file (locks/{sample}.lock)
with the DynamoDB-native equivalent of the same claim pattern.
"""

import os
import re
from datetime import datetime, timezone

import boto3
from botocore.exceptions import ClientError

S3 = boto3.client("s3")
BATCH = boto3.client("batch")
TABLE = boto3.resource("dynamodb").Table(os.environ["STATUS_TABLE_NAME"])

JOB_QUEUE = os.environ["BATCH_JOB_QUEUE"]
JOB_DEFINITION = os.environ["BATCH_JOB_DEFINITION"]

LANDING_PREFIX = "landing/"
FASTQ_EXTENSIONS = (".fastq", ".fastq.gz")


def _now():
    return datetime.now(timezone.utc).isoformat()


def _parse_sample(key):
    """Mirrors the Snakefile's own SAMPLES parsing (i.split("_end")[0]) so a
    key is only treated as a landing-zone upload if it would also resolve to
    a valid sample there."""
    if not key.startswith(LANDING_PREFIX):
        return None
    rel = key[len(LANDING_PREFIX):]
    if "/" not in rel:
        return None
    sample_dir, filename = rel.split("/", 1)
    if not filename.endswith(FASTQ_EXTENSIONS) or "_end" not in filename:
        return None
    sample = filename.split("_end")[0]
    return sample if sample == sample_dir else None


def _both_mates_present(bucket, sample):
    resp = S3.list_objects_v2(Bucket=bucket, Prefix=f"{LANDING_PREFIX}{sample}/{sample}_end")
    keys = [obj["Key"] for obj in resp.get("Contents", [])]
    has_end1 = any("_end1" in k for k in keys)
    has_end2 = any("_end2" in k for k in keys)
    return has_end1 and has_end2


def _try_claim(sample_id):
    """Atomically flips status to SUBMITTED. Returns True iff this invocation won the race."""
    try:
        TABLE.update_item(
            Key={"sample_id": sample_id},
            UpdateExpression="SET #s = :submitted, updated_at = :now",
            ConditionExpression="attribute_not_exists(sample_id) OR #s <> :submitted",
            ExpressionAttributeNames={"#s": "status"},
            ExpressionAttributeValues={":submitted": "SUBMITTED", ":now": _now()},
        )
        return True
    except ClientError as e:
        if e.response["Error"]["Code"] == "ConditionalCheckFailedException":
            return False
        raise


def _batch_job_name(sample):
    # Batch job names must match ^[a-zA-Z0-9][a-zA-Z0-9_-]{0,127}$ - sample
    # names (e.g. "Kunkel_Ribo-seq_Pol2MGrnh201.1b.1") may contain dots or
    # other characters that aren't valid there, even though they're fine in
    # DynamoDB keys and S3 paths.
    return re.sub(r"[^A-Za-z0-9_-]", "-", f"hyden-seq-{sample}")


def _submit_batch_job(bucket, sample):
    input_uri = f"s3://{bucket}/landing/{sample}/"
    output_uri = f"s3://{bucket}/output/{sample}/"
    # batch_entrypoint.sh's own default; passed explicitly rather than as ""
    # because Batch's submit_job API rejects an empty string command element
    # outright, so the shell script's own "${3:-default}" fallback never
    # gets a chance to run.
    reference_uri = f"s3://{bucket}/reference/"
    processed_raw_uri = f"s3://{bucket}/processed_raw/{sample}/"
    resp = BATCH.submit_job(
        jobName=_batch_job_name(sample),
        jobQueue=JOB_QUEUE,
        jobDefinition=JOB_DEFINITION,
        # positional command = [input_uri, output_uri, reference_uri, processed_raw_uri]
        containerOverrides={"command": [input_uri, output_uri, reference_uri, processed_raw_uri]},
    )
    return resp["jobId"]


def handler(event, context):
    for record in event["Records"]:
        bucket = record["s3"]["bucket"]["name"]
        key = record["s3"]["object"]["key"]

        sample = _parse_sample(key)
        if sample is None:
            continue

        if not _both_mates_present(bucket, sample):
            continue  # still waiting on the other mate

        if not _try_claim(sample):
            continue  # the other mate's invocation already submitted this job

        job_id = _submit_batch_job(bucket, sample)
        TABLE.update_item(
            Key={"sample_id": sample},
            UpdateExpression="SET job_id = :j",
            ExpressionAttributeValues={":j": job_id},
        )
