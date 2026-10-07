{% macro binned_watson_crick_ratio(bin_size) %}

with binned as (
    select
        sample_name,
        strand,
        chrom,
        floor(start_pos / {{ bin_size }}) * {{ bin_size }} as bin_start,
        score
    from {{ ref('stg_bedgraph_signal') }}
),

per_bin as (
    select
        sample_name,
        chrom,
        bin_start,
        sum(case when strand = 'forward' then score else 0 end) as watson,
        sum(case when strand = 'reverse' then score else 0 end) as crick
    from binned
    group by 1, 2, 3
)

select
    sample_name,
    chrom,
    bin_start,
    watson,
    crick,
    log(2, (watson + 1) / (crick + 1)) as log2_watson_crick
from per_bin

{% endmacro %}
