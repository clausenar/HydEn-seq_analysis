-- Share of each sample's signal on the forward strand, over the whole genome.
WITH signal AS (
    SELECT
        regexp_extract("$path", '([^/]+)__(forward|reverse)\.bedgraph$', 1) AS sample_name,
        regexp_extract("$path", '([^/]+)__(forward|reverse)\.bedgraph$', 2) AS strand,
        score
    FROM bedgraph_signal_raw
)
SELECT
    sample_name,
    sum(score) FILTER (WHERE strand = 'forward') AS forward_score,
    sum(score) FILTER (WHERE strand = 'reverse') AS reverse_score,
    sum(score) FILTER (WHERE strand = 'forward') / sum(score) AS forward_share
FROM signal
GROUP BY sample_name
ORDER BY sample_name
