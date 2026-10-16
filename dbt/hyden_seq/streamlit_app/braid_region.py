"""Region braid plot: log2(Watson / Crick) for one window of one chromosome,
with a bin size chosen at query time.

Kept apart from app.py so the query, the binning and the plot can be tested
without Streamlit or a live Snowflake session (see test_braid_region.py).
The binning matches braid_plot.py and the dbt binned_watson_crick_ratio
macro: each bedGraph interval is counted in the bin holding its start
position, empty bins count as zero on both strands, and the ratio uses a
pseudocount of 1.
"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIGNAL_TABLE = "HYDEN_SEQ.ANALYTICS_STAGING.STG_BEDGRAPH_SIGNAL"

# Above this many bins per sample a plot is unreadable and the query slow;
# e.g. 1bp bins are fine for a 20kb window but not for a whole chromosome.
MAX_BINS = 20_000


def region_bounds(centre, half_window, bin_size, chrom_len):
    """Returns (start, end) of the window around centre, clipped to the
    chromosome. start is rounded down to a multiple of bin_size, so the bins
    line up with the genome-wide 500bp grid and with braid_plot.py's."""
    if bin_size < 1:
        raise ValueError("Bin size must be at least 1 bp.")
    if half_window < 1:
        raise ValueError("Window must be at least 1 bp.")
    if not 0 <= centre < chrom_len:
        raise ValueError(f"Centre must be between 0 and {chrom_len - 1:,} on this chromosome.")
    start = max(0, centre - half_window) // bin_size * bin_size
    end = min(chrom_len, centre + half_window)
    n_bins = -(-(end - start) // bin_size)  # ceil division
    if n_bins > MAX_BINS:
        raise ValueError(
            f"{n_bins:,} bins is too many to plot (limit {MAX_BINS:,}). "
            "Use a larger bin size or a smaller window."
        )
    return start, end


def region_query(n_samples):
    """SQL with ? placeholders, bound in this order:
    bin_size, bin_size, chrom, start, end, then one sample name per sample."""
    sample_marks = ", ".join(["?"] * n_samples)
    return f"""
        SELECT
            sample_name,
            FLOOR(start_pos / ?) * ? AS bin_start,
            SUM(IFF(strand = 'forward', score, 0)) AS watson,
            SUM(IFF(strand = 'reverse', score, 0)) AS crick
        FROM {SIGNAL_TABLE}
        WHERE chrom = ?
          AND start_pos >= ? AND start_pos < ?
          AND sample_name IN ({sample_marks})
        GROUP BY 1, 2
    """


def region_params(samples, chrom, start, end, bin_size):
    return [bin_size, bin_size, chrom, start, end, *samples]


def dense_bins(rows, samples, start, end, bin_size):
    """One row per (sample, bin) across the whole window, empty bins as zero.

    rows: DataFrame with columns SAMPLE_NAME, BIN_START, WATSON, CRICK, holding
    only the bins that had signal (what region_query returns).
    """
    bin_starts = np.arange(start, end, bin_size)
    full = pd.MultiIndex.from_product([samples, bin_starts], names=["SAMPLE_NAME", "BIN_START"])
    if len(rows):
        rows = rows.astype({"BIN_START": "int64", "WATSON": "float64", "CRICK": "float64"})
        dense = rows.set_index(["SAMPLE_NAME", "BIN_START"])[["WATSON", "CRICK"]]
        dense = dense.reindex(full, fill_value=0.0)
    else:
        dense = pd.DataFrame(0.0, index=full, columns=["WATSON", "CRICK"])
    dense = dense.reset_index()
    dense["LOG2_WATSON_CRICK"] = np.log2((dense["WATSON"] + 1) / (dense["CRICK"] + 1))
    return dense


def plot_region(dense, samples, chrom, start, end, bin_size, origins=()):
    """One braid panel per sample, sharing the x axis; red where Watson > Crick,
    blue where Crick > Watson, dashed lines at origin midpoints in the window."""
    fig, axes = plt.subplots(
        len(samples), 1, figsize=(14, max(2.2 * len(samples), 3)), sharex=True, squeeze=False
    )
    axes = axes[:, 0]
    limit = max(1.0, float(np.abs(dense["LOG2_WATSON_CRICK"]).max())) if len(dense) else 1.0
    for ax, sample in zip(axes, samples):
        sub = dense[dense["SAMPLE_NAME"] == sample]
        pos = sub["BIN_START"].to_numpy()
        ratio = sub["LOG2_WATSON_CRICK"].to_numpy()
        ax.fill_between(pos, ratio, 0, where=(ratio >= 0), color="tab:red", linewidth=0, step="post")
        ax.fill_between(pos, ratio, 0, where=(ratio < 0), color="tab:blue", linewidth=0, step="post")
        ax.axhline(0, color="black", linewidth=0.5)
        for mid in origins:
            if start <= mid < end:
                ax.axvline(mid, color="black", linewidth=0.7, alpha=0.6, linestyle="--")
        ax.set_xlim(start, end)
        ax.set_ylim(-limit * 1.05, limit * 1.05)
        ax.set_ylabel(f"{sample}\nlog2(W/C)", rotation=0, labelpad=45, va="center", fontsize=8)
    axes[-1].set_xlabel(f"{chrom} position (bp)")
    axes[-1].ticklabel_format(axis="x", style="plain")
    fig.suptitle(f"{chrom}:{start:,}-{end:,}, {bin_size:,}bp bins")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig
