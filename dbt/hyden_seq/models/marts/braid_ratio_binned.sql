-- SQL-native version of braid_plot.py's core computation: bins the genome
-- into fixed-size windows and computes log2(Watson/Crick) per bin, per
-- sample. Bin size defaults to 500bp, matching the script's own default -
-- override at build time with `dbt run --vars '{bin_size: 50}'`.
--
-- For a standing fine-resolution cut, see braid_ratio_binned_5bp.sql
-- instead (a view, fixed at 5bp, doesn't get overwritten by --vars runs
-- against this model).

{{ binned_watson_crick_ratio(var('bin_size', 500)) }}
