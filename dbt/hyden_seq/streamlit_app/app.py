"""Braid plot + origin metaplot, live from Snowflake - Streamlit-in-Snowflake
reproduction of this repo's braid_plot.py and origin_metaplot.py, querying
HYDEN_SEQ.ANALYTICS_MARTS directly instead of reading local bedGraph files.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st
from snowflake.snowpark.context import get_active_session

session = get_active_session()

st.set_page_config(layout="wide")

samples = [
    r[0]
    for r in session.sql(
        "SELECT DISTINCT sample_name FROM HYDEN_SEQ.ANALYTICS_MARTS.BRAID_RATIO_BINNED ORDER BY 1"
    ).collect()
]
sample = st.sidebar.selectbox("Sample", samples)
page = st.sidebar.radio("View", ["Braid plot", "Origin metaplot"])

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
        data = session.sql(f"""
            SELECT chrom, bin_start, log2_watson_crick
            FROM HYDEN_SEQ.ANALYTICS_MARTS.BRAID_RATIO_BINNED
            WHERE sample_name = '{sample}'
            ORDER BY chrom, bin_start
        """).to_pandas()
        chroms = [c for c in chrom_sizes if c in set(data["CHROM"])]
        bin_size_label = "500bp"
    else:
        chrom = st.sidebar.selectbox("Chromosome", list(chrom_sizes.keys()))
        data = session.sql(f"""
            SELECT chrom, bin_start, log2_watson_crick
            FROM HYDEN_SEQ.ANALYTICS_MARTS.BRAID_RATIO_BINNED_5BP
            WHERE sample_name = '{sample}' AND chrom = '{chrom}'
            ORDER BY bin_start
        """).to_pandas()
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
            session.sql(f"""
                SELECT strand, chrom, n_positions, total_score
                FROM HYDEN_SEQ.ANALYTICS_MARTS.BEDGRAPH_SIGNAL_SUMMARY
                WHERE sample_name = '{sample}'
                ORDER BY strand, chrom
            """).to_pandas()
        )


def origin_metaplot_page():
    st.title("Origin metaplot: signal and strand ratio around replication origins")

    window = st.sidebar.number_input("Window (+/- bp around origin)", value=2000, step=500)
    bin_size = st.sidebar.number_input("Bin size (bp)", value=50, step=10)

    long_df = session.sql(f"""
        SELECT origin_id, strand, bin_start, score_sum
        FROM HYDEN_SEQ.ANALYTICS_MARTS.ORIGIN_SIGNAL_BINNED
        WHERE sample_name = '{sample}'
    """).to_pandas()
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
else:
    origin_metaplot_page()
