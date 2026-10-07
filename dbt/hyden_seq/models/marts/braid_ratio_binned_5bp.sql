-- Same computation as braid_ratio_binned.sql, fixed at a 5bp bin size -
-- effectively per-base resolution given the source bedGraphs are
-- themselves ~1bp intervals, so this is about as fine-grained as the data
-- actually supports. Materialized as a view rather than a table: at 5bp
-- this is ~100x braid_ratio_binned's row count, and it's a resolution
-- you'd query occasionally to zoom into one region, not scan wholesale -
-- not worth paying to store precomputed.

{{ config(materialized='view') }}

{{ binned_watson_crick_ratio(5) }}
