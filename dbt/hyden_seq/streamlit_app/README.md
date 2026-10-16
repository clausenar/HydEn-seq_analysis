# HydEn-seq Streamlit app (Streamlit in Snowflake)

Braid plots and origin metaplots drawn live from the Snowflake tables that dbt
builds, instead of from local bedGraph files.

| Page | What it shows | Source |
|---|---|---|
| Braid plot | log2(Watson / Crick) genome-wide at 500bp, or one chromosome at 5bp | `ANALYTICS_MARTS.BRAID_RATIO_BINNED`, `BRAID_RATIO_BINNED_5BP` |
| Region braid plot | one window of one chromosome at any bin size, one or more samples stacked, origins marked; PNG and CSV download | `ANALYTICS_STAGING.STG_BEDGRAPH_SIGNAL`, binned at query time |
| Origin metaplot | mean signal and strand ratio around replication origins | `ANALYTICS_MARTS.ORIGIN_SIGNAL_BINNED` |

The region page takes a centre position, a window (± bp) and a bin size. The
window's start is rounded down to a multiple of the bin size, so its bins line
up with the genome-wide ones. One plot is limited to 20,000 bins per sample
(for example 1bp bins over a 20kb window). Binning matches `braid_plot.py`: each
1bp bedGraph interval counts in the bin holding its start, empty bins are zero
on both strands, and the ratio uses a pseudocount of 1. Checked on 2026-10-16:
every bin of three windows (chrIV middle at 100bp, chrI start at 1bp, chrXII
end at 500bp, both samples) matched `braid_plot.py`'s own functions run on the
bedGraph files.

## Files

- `app.py`: the pages. All queries use bound parameters (`params=[...]`).
- `braid_region.py`: the region page's query, bin filling and plot, with no
  Streamlit or Snowflake dependency.
- `test_braid_region.py`: unit tests for `braid_region.py`
  (`python -m pytest -q test_braid_region.py`, needs pandas, numpy, matplotlib).
- `environment.yml`: packages the app gets in Snowflake.

## Deploying a change

There is no Snowflake CLI on the development Mac, so changes go in through
Snowsight: open the app under **Projects → Streamlit**, choose **Edit**, replace
`app.py`, add or replace `braid_region.py` in the file list, then **Run**. The
test file doesn't need to be uploaded. The app needs an active (paid)
Snowflake account, and its warehouse `HYDEN_WH` runs while the app is open.
