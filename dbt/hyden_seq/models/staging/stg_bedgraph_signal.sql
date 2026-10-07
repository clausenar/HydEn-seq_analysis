-- Cleans and types the raw landing table: parses sample_name/strand out of
-- the S3 path (they're not columns in the source bedGraph files themselves),
-- and casts positions/score explicitly rather than relying on COPY INTO's
-- inferred types.

with source as (
    select * from {{ source('pipeline', 'bedgraph_signal_raw') }}
)

select
    regexp_substr(src_file, '([^/]+)__(forward|reverse)\\.bedgraph$', 1, 1, 'e', 1) as sample_name,
    regexp_substr(src_file, '([^/]+)__(forward|reverse)\\.bedgraph$', 1, 1, 'e', 2) as strand,
    chrom,
    start_pos::number as start_pos,
    end_pos::number as end_pos,
    score::float as score,
    loaded_at
from source
