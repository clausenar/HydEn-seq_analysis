-- Snowflake reads the dev Iceberg tables in place, through the Glue Data Catalog.
-- Run in Snowsight as ACCOUNTADMIN, one statement at a time, in order.
-- The AWS side is terraform/modules/pipeline/snowflake.tf (role
-- hyden-seq-dev-snowflake-iceberg). Prod: iceberg_glue_prod.sql, same steps.
--
-- Never CREATE OR REPLACE the integration or the volume once linked: that
-- generates new external IDs and breaks the role's trust. Use ALTER.

-- ---------------------------------------------------------------- step 1
-- Before: terraform apply in envs/dev (creates the role, not yet trusted).

CREATE CATALOG INTEGRATION HYDEN_GLUE_DEV_INT
  CATALOG_SOURCE = GLUE
  CATALOG_NAMESPACE = 'hyden_seq_dev'
  TABLE_FORMAT = ICEBERG
  GLUE_AWS_ROLE_ARN = 'arn:aws:iam::069509443906:role/hyden-seq-dev-snowflake-iceberg'
  GLUE_CATALOG_ID = '069509443906'
  GLUE_REGION = 'us-east-1'
  REFRESH_INTERVAL_SECONDS = 60
  ENABLED = TRUE;

CREATE EXTERNAL VOLUME HYDEN_ICEBERG_DEV_VOL
  STORAGE_LOCATIONS = ((
    NAME = 'hyden-seq-dev-warehouse'
    STORAGE_PROVIDER = 'S3'
    STORAGE_BASE_URL = 's3://hyden-seq-dev-data-90322281e63f81697d01cae60c/warehouse/'
    STORAGE_AWS_ROLE_ARN = 'arn:aws:iam::069509443906:role/hyden-seq-dev-snowflake-iceberg'
  ))
  ALLOW_WRITES = FALSE;

-- Copy GLUE_AWS_IAM_USER_ARN and GLUE_AWS_EXTERNAL_ID from this:
DESC CATALOG INTEGRATION HYDEN_GLUE_DEV_INT;

-- and STORAGE_AWS_IAM_USER_ARN and STORAGE_AWS_EXTERNAL_ID from the
-- STORAGE_LOCATION_1 row of this:
DESC EXTERNAL VOLUME HYDEN_ICEBERG_DEV_VOL;

-- ---------------------------------------------------------------- step 2
-- Before: those four values go into snowflake_iceberg_trust in
-- terraform/envs/dev/main.tf, and terraform apply.

-- Should return "success": true.
SELECT SYSTEM$VERIFY_EXTERNAL_VOLUME('HYDEN_ICEBERG_DEV_VOL');

CREATE SCHEMA IF NOT EXISTS HYDEN_SEQ.ICEBERG_DEV
  COMMENT = 'Iceberg tables written by the dev Glue job, read through the Glue Data Catalog';

-- AUTO_REFRESH: Snowflake polls the Glue catalog (REFRESH_INTERVAL_SECONDS
-- above) and picks up each new snapshot the Glue job commits.
CREATE ICEBERG TABLE HYDEN_SEQ.ICEBERG_DEV.BEDGRAPH_SIGNAL
  EXTERNAL_VOLUME = 'HYDEN_ICEBERG_DEV_VOL'
  CATALOG = 'HYDEN_GLUE_DEV_INT'
  CATALOG_TABLE_NAME = 'bedgraph_signal'
  AUTO_REFRESH = TRUE;

CREATE ICEBERG TABLE HYDEN_SEQ.ICEBERG_DEV.STRAND_BALANCE
  EXTERNAL_VOLUME = 'HYDEN_ICEBERG_DEV_VOL'
  CATALOG = 'HYDEN_GLUE_DEV_INT'
  CATALOG_TABLE_NAME = 'strand_balance'
  AUTO_REFRESH = TRUE;

-- ---------------------------------------------------------------- step 3: check
-- Expect SRR1609190sub5 and SRR1609190sub6, each 99,204 forward / 99,729
-- reverse rows, score 138,880 / 143,728 (the same as Athena).
SELECT sample_name, strand, COUNT(*) AS n_positions, SUM(score) AS total_score
FROM HYDEN_SEQ.ICEBERG_DEV.BEDGRAPH_SIGNAL
GROUP BY 1, 2
ORDER BY 1, 2;

SELECT sample_name, SUM(forward_score) / (SUM(forward_score) + SUM(reverse_score)) AS forward_share
FROM HYDEN_SEQ.ICEBERG_DEV.STRAND_BALANCE
GROUP BY 1
ORDER BY 1;
