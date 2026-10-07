-- SQL-native version of origin_metaplot.py's core windowing step: for every
-- origin, bins the +/-window bp around its midpoint into fixed-size bins
-- (position expressed *relative to the origin*, so bin_start=0 is always
-- the origin midpoint itself, same as the script's BIN_CENTERS), and sums
-- each strand's signal per origin per bin.
--
-- One row per (origin, sample, strand, relative bin) - the long-form
-- equivalent of the script's watson/crick (origins x bins) matrices.
-- Final mean-across-origins aggregation (the actual metaplot line) is left
-- to the querying app/notebook, same division of labor as the rest of this
-- project: the warehouse does the expensive per-position join, small
-- aggregations happen client-side.

{% set window = var('origin_window', 2000) %}
{% set bin_size = var('origin_bin_size', 50) %}

with origins as (
    select
        row_number() over (order by chrom, midpoint) as origin_id,
        chrom,
        midpoint
    from {{ source('pipeline', 'origins') }}
),

signal as (
    select * from {{ ref('stg_bedgraph_signal') }}
),

joined as (
    select
        o.origin_id,
        o.chrom,
        o.midpoint,
        s.sample_name,
        s.strand,
        s.start_pos - o.midpoint as relative_pos,
        s.score
    from origins o
    inner join signal s
        on s.chrom = o.chrom
       and s.start_pos between o.midpoint - {{ window }} and o.midpoint + {{ window }}
)

select
    origin_id,
    chrom,
    midpoint,
    sample_name,
    strand,
    floor(relative_pos / {{ bin_size }}) * {{ bin_size }} as bin_start,
    sum(score) as score_sum
from joined
group by 1, 2, 3, 4, 5, 6
