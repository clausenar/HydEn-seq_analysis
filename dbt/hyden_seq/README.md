# hyden_seq dbt project

Transforms the raw bedGraph data landed in Snowflake (see the main repo's
[README, "Snowflake warehouse"](../../README.md#snowflake-warehouse) section
for how it gets there) into a small staging → marts model, replacing a single
hand-written view with a real, tested dbt project.

## Layout

```
models/
  staging/
    _sources.yml            declares HYDEN_SEQ.PIPELINE.BEDGRAPH_SIGNAL_RAW
                             as a dbt source (landed outside dbt, via COPY INTO)
    stg_bedgraph_signal.sql  parses sample_name/strand out of the S3 path,
                             casts types - the dbt-managed replacement for
                             the hand-written BEDGRAPH_SIGNAL view
    _staging.yml             not_null/accepted_values/uniqueness/score>=0 tests
  marts/
    bedgraph_signal_summary.sql  per-sample/strand/chrom totals
    braid_ratio_binned.sql       SQL-native log2(Watson/Crick), binned
                                  genome-wide - the warehouse-side equivalent
                                  of braid_plot.py's core computation
    _marts.yml                    tests on both marts
```

Staging models materialize as views (schema `ANALYTICS_staging`); marts
materialize as tables (schema `ANALYTICS_marts`) - the default schema
(`ANALYTICS`) plus dbt's own `<custom_schema>` suffixing, set in
`dbt_project.yml`.

## Running it

Needs a `~/.dbt/profiles.yml` with key-pair auth configured (private key at
`~/.snowflake/rsa_key.p8`, registered against the Snowflake user via
`ALTER USER ... SET RSA_PUBLIC_KEY=...` - see the main README for the
account-level setup). Then, from this directory:

```bash
dbt deps    # installs dbt_utils (used for the composite-uniqueness and
            # score>=0 tests)
dbt debug   # verify the Snowflake connection
dbt build   # runs every model and every test
```

`braid_ratio_binned`'s bin size defaults to 500bp (matching `braid_plot.py`'s
own default) and can be overridden per run:

```bash
dbt run --select braid_ratio_binned --vars '{bin_size: 50}'
```

## Why this exists

The pipeline itself doesn't depend on any of this - it's a practice/portfolio
layer on top of the Snowflake warehouse setup, demonstrating a proper
staging/marts split with source declarations and automated tests instead of
one hand-written `CREATE VIEW`.
