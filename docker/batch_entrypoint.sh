#!/bin/bash
set -euo pipefail
export PATH="/opt/conda/bin:$PATH"

USAGE="Usage: batch_entrypoint.sh s3://input-prefix/ s3://output-prefix/ s3://reference-prefix/ [s3://processed-raw-prefix/]"
S3_INPUT_PREFIX="${1:?$USAGE}"
S3_OUTPUT_PREFIX="${2:?$USAGE}"
# Required, with no default: each stage has its own bucket, and a default
# naming one of them would quietly read another stage's reference data.
S3_REFERENCE_PREFIX="${3:?$USAGE}"
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

# The Glue/Athena table reads every file under its S3 location, and output/
# also holds BAMs and plots, so the two bedGraphs get a prefix of their own.
# Set on the job definition by Terraform; the file names already carry the
# sample and strand, so a rerun of a sample overwrites rather than duplicates.
if [ -n "${HYDEN_S3_BEDGRAPH_TABLE_PREFIX:-}" ]; then
  echo "Publishing bedGraphs to ${HYDEN_S3_BEDGRAPH_TABLE_PREFIX} ..."
  for f in /data/output/*__forward.bedgraph /data/output/*__reverse.bedgraph; do
    aws s3 cp "$f" "${HYDEN_S3_BEDGRAPH_TABLE_PREFIX}$(basename "$f")"
  done
fi

if [ -n "${S3_PROCESSED_RAW_PREFIX}" ]; then
  echo "Moving raw input from ${S3_INPUT_PREFIX} to ${S3_PROCESSED_RAW_PREFIX} ..."
  aws s3 mv "${S3_INPUT_PREFIX}" "${S3_PROCESSED_RAW_PREFIX}" --recursive
fi

echo "Done."
