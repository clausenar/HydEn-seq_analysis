#!/bin/bash
set -euo pipefail
export PATH="/opt/conda/bin:$PATH"

S3_INPUT_PREFIX="${1:?Usage: batch_entrypoint.sh s3://input-prefix/ s3://output-prefix/ [s3://reference-prefix/] [s3://processed-raw-prefix/]}"
S3_OUTPUT_PREFIX="${2:?Usage: batch_entrypoint.sh s3://input-prefix/ s3://output-prefix/ [s3://reference-prefix/] [s3://processed-raw-prefix/]}"
S3_REFERENCE_PREFIX="${3:-s3://kunkel-ribo-data-2026/reference/}"
S3_PROCESSED_RAW_PREFIX="${4:-}"

# Fixed container-internal layout; the Snakefile/scripts pick these up via
# the HYDEN_* environment variables added for exactly this purpose, so this
# is the *unmodified* git-tracked pipeline code, not a patched fork of it.
export HYDEN_RAW_DIR=/data/raw
export HYDEN_GENOME=/data/reference/genome/sacCer3
export HYDEN_OLIGOS=/data/reference/oligo/oligo
export HYDEN_OUT_DIR=/data/output
export HYDEN_ORIGINS_FILE=/data/reference/or200.txt

echo "Syncing reference data from ${S3_REFERENCE_PREFIX} ..."
mkdir -p /data/reference /data/raw /data/output
aws s3 sync "${S3_REFERENCE_PREFIX}" /data/reference

echo "Syncing input from ${S3_INPUT_PREFIX} ..."
aws s3 sync "${S3_INPUT_PREFIX}" /data/raw

echo "Running pipeline ..."
cd /app
conda run --no-capture-output -n hyden_pipeline snakemake --cores 8 -p

echo "Syncing output to ${S3_OUTPUT_PREFIX} ..."
aws s3 sync /data/output "${S3_OUTPUT_PREFIX}"

if [ -n "${S3_PROCESSED_RAW_PREFIX}" ]; then
  echo "Moving raw input from ${S3_INPUT_PREFIX} to ${S3_PROCESSED_RAW_PREFIX} ..."
  aws s3 mv "${S3_INPUT_PREFIX}" "${S3_PROCESSED_RAW_PREFIX}" --recursive
fi

echo "Done."
