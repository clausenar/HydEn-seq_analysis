-- Positions and signal per sample, strand and chromosome. The same columns as
-- the dbt mart bedgraph_signal_summary, so the two can be compared row by row.
-- Sample and strand come from the file name, as in stg_bedgraph_signal.
SELECT
    regexp_extract("$path", '([^/]+)__(forward|reverse)\.bedgraph$', 1) AS sample_name,
    regexp_extract("$path", '([^/]+)__(forward|reverse)\.bedgraph$', 2) AS strand,
    chrom,
    count(*) AS n_positions,
    sum(score) AS total_score,
    avg(score) AS mean_score,
    max(score) AS max_score
FROM bedgraph_signal_raw
GROUP BY 1, 2, 3
ORDER BY 1, 2, 3
