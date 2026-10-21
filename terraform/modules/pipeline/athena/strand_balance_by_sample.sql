-- Gold table written by the Glue job (Iceberg): forward share per sample over
-- the whole genome, from the per-chromosome rows.
SELECT
    sample_name,
    sum(n_positions) AS n_positions,
    sum(forward_score) AS forward_score,
    sum(reverse_score) AS reverse_score,
    sum(forward_score) / (sum(forward_score) + sum(reverse_score)) AS forward_share,
    max(processed_at) AS processed_at
FROM strand_balance
GROUP BY sample_name
ORDER BY sample_name
