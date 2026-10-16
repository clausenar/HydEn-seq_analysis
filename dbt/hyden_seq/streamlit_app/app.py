"""Braid plot + origin metaplot, live from Snowflake - Streamlit-in-Snowflake
reproduction of this repo's braid_plot.py and origin_metaplot.py, querying
HYDEN_SEQ.ANALYTICS_MARTS directly instead of reading local bedGraph files.
The region braid plot computes its bins at query time from the staging view,
so the window and bin size can be anything (see braid_region.py).
"""

import io

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st
from snowflake.snowpark.context import get_active_session

import braid_region

session = get_active_session()

st.set_page_config(layout="wide")

samples = [
    r[0]
    for r in session.sql(
        "SELECT DISTINCT sample_name FROM HYDEN_SEQ.ANALYTICS_MARTS.BRAID_RATIO_BINNED ORDER BY 1"
    ).collect()
]
sample = st.sidebar.selectbox("Sample", samples)
page = st.sidebar.radio("View", ["Braid plot", "Region braid plot", "Origin metaplot"])

chrom_sizes_df = session.sql(
    "SELECT chrom, length FROM HYDEN_SEQ.PIPELINE.CHROM_SIZES"
).to_pandas()
chrom_sizes = dict(zip(chrom_sizes_df["CHROM"], chrom_sizes_df["LENGTH"]))


def braid_plot_page():
    st.title("Braid plot: log2(Watson / Crick) ribonucleotide ratio")

    show_origins = st.sidebar.checkbox("Show origins", value=True)
    origins_df = session.sql("SELECT chrom, midpoint FROM HYDEN_SEQ.PIPELINE.ORIGINS").to_pandas()
    origins = {}
    for _, row in origins_df.iterrows():
        origins.setdefault(row["CHROM"], []).append(row["MIDPOINT"])

    resolution = st.sidebar.radio("Resolution", ["500bp (genome-wide)", "5bp (single chromosome)"])

    if resolution.startswith("500bp"):
        data = session.sql("""
            SELECT chrom, bin_start, log2_watson_crick
            FROM HYDEN_SEQ.ANALYTICS_MARTS.BRAID_RATIO_BINNED
            WHERE sample_name = ?
            ORDER BY chrom, bin_start
        """, params=[sample]).to_pandas()
        chroms = [c for c in chrom_sizes if c in set(data["CHROM"])]
        bin_size_label = "500bp"
    else:
        chrom = st.sidebar.selectbox("Chromosome", list(chrom_sizes.keys()))
        data = session.sql("""
            SELECT chrom, bin_start, log2_watson_crick
            FROM HYDEN_SEQ.ANALYTICS_MARTS.BRAID_RATIO_BINNED_5BP
            WHERE sample_name = ? AND chrom = ?
            ORDER BY bin_start
        """, params=[sample, chrom]).to_pandas()
        chroms = [chrom]
        bin_size_label = "5bp"

    st.caption(f"{sample} - {bin_size_label} bins - {len(data):,} rows")

    max_len = max(chrom_sizes[c] for c in chroms)
    fig, axes = plt.subplots(len(chroms), 1, figsize=(16, max(1.4 * len(chroms), 3)), sharex=False)
    axes = [axes] if len(chroms) == 1 else axes

    for ax, c in zip(axes, chroms):
        sub = data[data["CHROM"] == c]
        positions = sub["BIN_START"].values
        ratio = sub["LOG2_WATSON_CRICK"].values

        ax.fill_between(positions, ratio, 0, where=(ratio >= 0), color="tab:red", linewidth=0)
        ax.fill_between(positions, ratio, 0, where=(ratio < 0), color="tab:blue", linewidth=0)
        ax.axhline(0, color="black", linewidth=0.5)
        if show_origins:
            for mid in origins.get(c, []):
                ax.axvline(mid, color="black", linewidth=0.5, alpha=0.5, linestyle="--")
        ax.set_xlim(0, chrom_sizes[c])
        ax.set_ylabel(c.replace("chr", ""), rotation=0, labelpad=20, va="center")
        ax.set_yticks([])

    axes[-1].set_xlabel("Position (bp)")
    fig.suptitle(f"{sample} - {bin_size_label} bins")
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    if len(chroms) > 1:
        for ax, c in zip(axes, chroms):
            pos = ax.get_position()
            width = pos.width * chrom_sizes[c] / max_len
            ax.set_position([pos.x0, pos.y0, width, pos.height])

    st.pyplot(fig)

    with st.expander("Sample totals"):
        st.dataframe(
            session.sql("""
                SELECT strand, chrom, n_positions, total_score
                FROM HYDEN_SEQ.ANALYTICS_MARTS.BEDGRAPH_SIGNAL_SUMMARY
                WHERE sample_name = ?
                ORDER BY strand, chrom
            """, params=[sample]).to_pandas()
        )


@st.cache_data(ttl=600, show_spinner=False)
def region_rows(region_samples, chrom, start, end, bin_size):
    """Bins with signal in the window, one row per (sample, bin). Cached for
    10 minutes so moving a widget back and forth doesn't re-run the query."""
    return session.sql(
        braid_region.region_query(len(region_samples)),
        params=braid_region.region_params(list(region_samples), chrom, start, end, bin_size),
    ).to_pandas()


def region_braid_page():
    st.title("Region braid plot: log2(Watson / Crick) in a chosen window")

    region_samples = st.sidebar.multiselect(
        "Samples (stacked)", samples, default=[sample],
        help="Pick several to compare them on the same axis.",
    )
    chrom = st.sidebar.selectbox("Chromosome", list(chrom_sizes.keys()))
    chrom_len = int(chrom_sizes[chrom])
    centre = st.sidebar.number_input(
        "Centre position (bp)", min_value=0, max_value=chrom_len - 1,
        value=min(chrom_len // 2, chrom_len - 1), step=1000,
    )
    half_window = st.sidebar.number_input(
        "Window (+/- bp around centre)", min_value=1, max_value=chrom_len, value=20000, step=1000,
    )
    bin_size = st.sidebar.number_input("Bin size (bp)", min_value=1, value=100, step=50)
    show_origins = st.sidebar.checkbox("Show origins", value=True)

    if not region_samples:
        st.info("Pick at least one sample.")
        return
    try:
        start, end = braid_region.region_bounds(int(centre), int(half_window), int(bin_size), chrom_len)
    except ValueError as e:
        st.warning(str(e))
        return

    with st.spinner("Querying Snowflake..."):
        rows = region_rows(tuple(region_samples), chrom, start, end, int(bin_size))
    dense = braid_region.dense_bins(rows, region_samples, start, end, int(bin_size))

    origins = []
    if show_origins:
        origins = [
            r[0]
            for r in session.sql(
                "SELECT midpoint FROM HYDEN_SEQ.PIPELINE.ORIGINS WHERE chrom = ? AND midpoint >= ? AND midpoint < ?",
                params=[chrom, start, end],
            ).collect()
        ]

    n_bins = len(dense) // len(region_samples)
    st.caption(
        f"{chrom}:{start:,}-{end:,} - {int(bin_size):,}bp bins - {n_bins:,} bins per sample - "
        f"{len(rows):,} bins with signal - {len(origins)} origins in window"
    )
    fig = braid_region.plot_region(dense, region_samples, chrom, start, end, int(bin_size), origins)
    st.pyplot(fig)

    png = io.BytesIO()
    fig.savefig(png, format="png", dpi=200)
    plt.close(fig)
    name = f"braid_{chrom}_{start}-{end}_{int(bin_size)}bp"
    col1, col2 = st.columns(2)
    col1.download_button("Download plot (PNG)", png.getvalue(), f"{name}.png", "image/png")
    col2.download_button(
        "Download bins (CSV)",
        dense.rename(columns=str.lower).to_csv(index=False).encode(),
        f"{name}.csv",
        "text/csv",
    )

    with st.expander("Binned values"):
        st.dataframe(dense)


def origin_metaplot_page():
    st.title("Origin metaplot: signal and strand ratio around replication origins")

    window = st.sidebar.number_input("Window (+/- bp around origin)", value=2000, step=500)
    bin_size = st.sidebar.number_input("Bin size (bp)", value=50, step=10)

    long_df = session.sql("""
        SELECT origin_id, strand, bin_start, score_sum
        FROM HYDEN_SEQ.ANALYTICS_MARTS.ORIGIN_SIGNAL_BINNED
        WHERE sample_name = ?
    """, params=[sample]).to_pandas()
    st.caption(f"{sample} - {window}bp window, {bin_size}bp bins - {long_df['ORIGIN_ID'].nunique():,} origins")

    # Dense-fill every (origin, bin) combination with 0, matching the script's
    # zero-initialized numpy arrays - origin_signal_binned only has rows where
    # a strand actually had reads in that bin, so without this a bin with true
    # zero signal would be silently dropped from the average instead of
    # counted as zero.
    bin_starts = np.arange(-window, window, bin_size)
    origin_ids = long_df["ORIGIN_ID"].unique()
    full_index = pd.MultiIndex.from_product([origin_ids, bin_starts], names=["ORIGIN_ID", "BIN_START"])

    pivot = long_df.pivot_table(
        index=["ORIGIN_ID", "BIN_START"], columns="STRAND", values="SCORE_SUM", aggfunc="sum"
    ).reindex(full_index, fill_value=0)
    for col in ("forward", "reverse"):
        if col not in pivot.columns:
            pivot[col] = 0
    pivot = pivot.fillna(0).reset_index()
    pivot["log2_ratio"] = np.log2((pivot["forward"] + 1) / (pivot["reverse"] + 1))

    def mean_sem(df, col):
        g = df.groupby("BIN_START")[col]
        n = g.count()
        return g.mean(), g.std() / np.sqrt(n)

    watson_mean, watson_sem = mean_sem(pivot, "forward")
    crick_mean, crick_sem = mean_sem(pivot, "reverse")
    ratio_mean, ratio_sem = mean_sem(pivot, "log2_ratio")

    col1, col2 = st.columns(2)

    with col1:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(watson_mean.index, watson_mean.values, color="tab:blue", label="Watson (+)")
        ax.fill_between(watson_mean.index, watson_mean - watson_sem, watson_mean + watson_sem, color="tab:blue", alpha=0.2)
        ax.plot(crick_mean.index, -crick_mean.values, color="tab:red", label="Crick (-)")
        ax.fill_between(crick_mean.index, -crick_mean - crick_sem, -crick_mean + crick_sem, color="tab:red", alpha=0.2)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_xlabel("Position relative to origin midpoint (bp)")
        ax.set_ylabel("Mean rNMP hit count per bin")
        ax.set_title("Metaplot")
        ax.legend()
        st.pyplot(fig)

    with col2:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(ratio_mean.index, ratio_mean.values, color="tab:purple", label="log2(Watson / Crick)")
        ax.fill_between(ratio_mean.index, ratio_mean - ratio_sem, ratio_mean + ratio_sem, color="tab:purple", alpha=0.2)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_xlabel("Position relative to origin midpoint (bp)")
        ax.set_ylabel("Mean log2(Watson / Crick)")
        ax.set_title("Strand ratio metaplot")
        ax.legend()
        st.pyplot(fig)


if page == "Braid plot":
    braid_plot_page()
elif page == "Region braid plot":
    region_braid_page()
else:
    origin_metaplot_page()
