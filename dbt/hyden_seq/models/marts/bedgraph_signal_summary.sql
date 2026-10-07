-- Per-sample, per-strand, per-chromosome totals - the aggregate view you'd
-- pull for a quick sanity check of a new sample's overall coverage.

select
    sample_name,
    strand,
    chrom,
    count(*) as n_positions,
    sum(score) as total_score,
    avg(score) as mean_score,
    max(score) as max_score
from {{ ref('stg_bedgraph_signal') }}
group by 1, 2, 3
