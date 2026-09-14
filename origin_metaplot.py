import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt

# --- CONFIGURATION (defaults; overridable via CLI args, see main()) ---
ORIGINS_FILE = "/Users/xranea/git/HydEn-seq_analysis/or200.txt"
FORWARD_BEDGRAPH = "/Users/xranea/bedgraphs/Kunkel_Ribo-seq_Pol2MGrnh201.1b.1__forward.bedgraph"
REVERSE_BEDGRAPH = "/Users/xranea/bedgraphs/Kunkel_Ribo-seq_Pol2MGrnh201.1b.1__reverse.bedgraph"
OUTPUT_DIR = "/Users/xranea/bedgraphs/processed_results"

WINDOW = 2000   # bp on each side of the origin midpoint
BIN_SIZE = 50   # bp per bin
# ---------------------

N_BINS = (2 * WINDOW) // BIN_SIZE
BIN_EDGES = np.arange(-WINDOW, WINDOW + BIN_SIZE, BIN_SIZE)
BIN_CENTERS = (BIN_EDGES[:-1] + BIN_EDGES[1:]) / 2

# or200.txt uses arabic chromosome numbers (chr1..chr16); the bedgraphs use
# roman numerals (chrI..chrXVI), matching the sacCer3 reference naming.
ROMAN = [
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
    "XI", "XII", "XIII", "XIV", "XV", "XVI",
]
CHROM_MAP = {f"chr{i+1}": f"chr{roman}" for i, roman in enumerate(ROMAN)}


def load_origins(path, exclude_names=None):
    """Reads the origin BED-like file and returns a list of (chrom, midpoint).

    If exclude_names is given, rows whose name column (4th column) is in it
    are skipped (e.g. {"null"} to drop origins with no assigned ARS name).
    """
    exclude_names = exclude_names or set()
    origins = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("track") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[:3]
            name = fields[3] if len(fields) > 3 else ""
            if name in exclude_names:
                continue
            chrom = CHROM_MAP.get(chrom, chrom)
            midpoint = (int(start) + int(end)) // 2
            origins.append((chrom, midpoint))
    return origins


def load_bedgraph(path):
    """Reads a bedgraph into {chrom: [(start, end, score), ...]}, sorted by start."""
    by_chrom = {}
    with open(path, "r") as f:
        for line in f:
            if line.startswith("track") or not line.strip():
                continue
            chrom, start, end, score = line.strip().split("\t")
            by_chrom.setdefault(chrom, []).append(
                (int(start), int(end), float(score))
            )
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x[0])
    return by_chrom


def bin_signal_for_window(intervals, window_start, window_end):
    """Accumulates per-bin totals for one origin window from overlapping bedgraph intervals."""
    row = np.zeros(N_BINS)
    if not intervals:
        return row

    starts = [iv[0] for iv in intervals]
    # Narrow down to intervals that could overlap the window.
    lo = max(0, np.searchsorted(starts, window_start - BIN_SIZE, side="left") - 1)
    hi = np.searchsorted(starts, window_end, side="right")

    for start, end, score in intervals[lo:hi]:
        if end <= window_start or start >= window_end:
            continue
        ov_start = max(start, window_start)
        ov_end = min(end, window_end)
        first_bin = (ov_start - window_start) // BIN_SIZE
        last_bin = (ov_end - 1 - window_start) // BIN_SIZE
        for b in range(int(first_bin), int(last_bin) + 1):
            bin_start = window_start + b * BIN_SIZE
            bin_end = bin_start + BIN_SIZE
            overlap = min(ov_end, bin_end) - max(ov_start, bin_start)
            if overlap > 0:
                row[b] += score * overlap
    return row


def build_matrices(origins, fwd_bedgraph, rev_bedgraph):
    watson = np.zeros((len(origins), N_BINS))
    crick = np.zeros((len(origins), N_BINS))

    for i, (chrom, mid) in enumerate(origins):
        window_start = mid - WINDOW
        window_end = mid + WINDOW
        watson[i] = bin_signal_for_window(
            fwd_bedgraph.get(chrom, []), window_start, window_end
        )
        crick[i] = bin_signal_for_window(
            rev_bedgraph.get(chrom, []), window_start, window_end
        )
    return watson, crick


HEATMAP_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "blue_red", ["blue", "red"]
)
HEATMAP_CMAP_R = HEATMAP_CMAP.reversed()
SPLIT_BIN = N_BINS // 2  # bin index of the origin midpoint (position 0)


def plot_heatmap(watson, crick, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 8), sharey=True)
    combined_nonzero = np.concatenate([watson[watson > 0], crick[crick > 0]])
    vmax = max(np.percentile(combined_nonzero, 99), 1) if combined_nonzero.size else 1

    strands = (watson, crick)
    titles = ("Watson (+)", "Crick (-)")
    # Watson: blue=low/red=high upstream, red=low/blue=high downstream (flips
    # at the origin). Crick uses the exact opposite of Watson on each side,
    # so the two panels are color-inverted mirrors of one another.
    left_cmaps = (HEATMAP_CMAP, HEATMAP_CMAP_R)
    right_cmaps = (HEATMAP_CMAP_R, HEATMAP_CMAP)
    # Canonical handles for the colorbars: Watson's images, since Watson
    # always uses HEATMAP_CMAP on the left and HEATMAP_CMAP_R on the right
    # regardless of the per-strand swap above.
    im_normal = im_reversed = None
    for i, (ax, matrix, title, left_cmap, right_cmap) in enumerate(
        zip(axes, strands, titles, left_cmaps, right_cmaps)
    ):
        n_rows = len(matrix)
        im_left = ax.imshow(
            matrix[:, :SPLIT_BIN],
            aspect="auto",
            cmap=left_cmap,
            vmin=0,
            vmax=vmax,
            extent=[-WINDOW, 0, n_rows, 0],
        )
        im_right = ax.imshow(
            matrix[:, SPLIT_BIN:],
            aspect="auto",
            cmap=right_cmap,
            vmin=0,
            vmax=vmax,
            extent=[0, WINDOW, n_rows, 0],
        )
        if i == 0:  # Watson
            im_normal, im_reversed = im_left, im_right
        ax.set_xlim(-WINDOW, WINDOW)
        ax.set_title(title)
        ax.set_xlabel("Position relative to origin midpoint (bp)")
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")

    axes[0].set_ylabel("Origins")
    fig.colorbar(
        im_normal, ax=axes, label="rNMP hit count per bin", shrink=0.8, extend="max",
    )
    fig.suptitle("Ribonucleotide incorporation around replication origins")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def compute_log2_ratio(watson, crick, pseudocount=1.0):
    """log2(Watson/Crick) per bin, with a pseudocount to avoid division by zero."""
    return np.log2((watson + pseudocount) / (crick + pseudocount))


def plot_ratio_heatmap(ratio, out_path):
    fig, ax = plt.subplots(figsize=(6, 8))
    vabs = max(np.percentile(np.abs(ratio), 99), 0.1)

    im = ax.imshow(
        ratio,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vabs,
        vmax=vabs,
        extent=[-WINDOW, WINDOW, len(ratio), 0],
    )
    ax.set_title("Watson / Crick strand ratio")
    ax.set_xlabel("Position relative to origin midpoint (bp)")
    ax.set_ylabel("Origins")
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
    fig.colorbar(im, ax=ax, label="log2(Watson / Crick)", extend="both")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_ratio_metaplot(ratio, out_path):
    ratio_mean = ratio.mean(axis=0)
    ratio_sem = ratio.std(axis=0) / np.sqrt(len(ratio))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(BIN_CENTERS, ratio_mean, color="tab:purple", label="log2(Watson / Crick)")
    ax.fill_between(
        BIN_CENTERS,
        ratio_mean - ratio_sem,
        ratio_mean + ratio_sem,
        color="tab:purple",
        alpha=0.2,
    )
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Position relative to origin midpoint (bp)")
    ax.set_ylabel("Mean log2(Watson / Crick)")
    ax.set_title("Strand ratio around origins")
    ax.legend()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_metaplot(watson, crick, out_path):
    watson_mean = watson.mean(axis=0)
    crick_mean = crick.mean(axis=0)
    watson_sem = watson.std(axis=0) / np.sqrt(len(watson))
    crick_sem = crick.std(axis=0) / np.sqrt(len(crick))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(BIN_CENTERS, watson_mean, color="tab:blue", label="Watson (+)")
    ax.fill_between(
        BIN_CENTERS,
        watson_mean - watson_sem,
        watson_mean + watson_sem,
        color="tab:blue",
        alpha=0.2,
    )
    ax.plot(BIN_CENTERS, -crick_mean, color="tab:red", label="Crick (-)")
    ax.fill_between(
        BIN_CENTERS,
        -crick_mean - crick_sem,
        -crick_mean + crick_sem,
        color="tab:red",
        alpha=0.2,
    )
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Position relative to origin midpoint (bp)")
    ax.set_ylabel("Mean rNMP hit count per bin")
    ax.set_title("Metaplot of ribonucleotide incorporation around origins")
    ax.legend()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def generate_outputs(origins, fwd_bedgraph, rev_bedgraph, suffix=""):
    """Builds matrices for the given origins and writes all plots/matrices,
    with output filenames tagged by `suffix` (e.g. "_null")."""
    watson, crick = build_matrices(origins, fwd_bedgraph, rev_bedgraph)

    heatmap_path = os.path.join(OUTPUT_DIR, f"origin_heatmap{suffix}.png")
    metaplot_path = os.path.join(OUTPUT_DIR, f"origin_metaplot{suffix}.png")
    plot_heatmap(watson, crick, heatmap_path)
    plot_metaplot(watson, crick, metaplot_path)

    ratio = compute_log2_ratio(watson, crick)
    ratio_heatmap_path = os.path.join(OUTPUT_DIR, f"origin_ratio_heatmap{suffix}.png")
    ratio_metaplot_path = os.path.join(OUTPUT_DIR, f"origin_ratio_metaplot{suffix}.png")
    plot_ratio_heatmap(ratio, ratio_heatmap_path)
    plot_ratio_metaplot(ratio, ratio_metaplot_path)

    np.savetxt(
        os.path.join(OUTPUT_DIR, f"origin_matrix_watson{suffix}.tsv"),
        watson,
        delimiter="\t",
        header="\t".join(str(int(c)) for c in BIN_CENTERS),
        comments="",
    )
    np.savetxt(
        os.path.join(OUTPUT_DIR, f"origin_matrix_crick{suffix}.tsv"),
        crick,
        delimiter="\t",
        header="\t".join(str(int(c)) for c in BIN_CENTERS),
        comments="",
    )

    print(f"  Heatmap saved to: {heatmap_path}")
    print(f"  Metaplot saved to: {metaplot_path}")
    print(f"  Ratio heatmap saved to: {ratio_heatmap_path}")
    print(f"  Ratio metaplot saved to: {ratio_metaplot_path}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--forward-bedgraph", default=FORWARD_BEDGRAPH)
    parser.add_argument("--reverse-bedgraph", default=REVERSE_BEDGRAPH)
    parser.add_argument("--output-dir", default=OUTPUT_DIR)
    parser.add_argument("--origins-file", default=ORIGINS_FILE)
    return parser.parse_args()


def main():
    args = parse_args()
    global OUTPUT_DIR
    OUTPUT_DIR = args.output_dir
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    fwd_bedgraph = load_bedgraph(args.forward_bedgraph)
    rev_bedgraph = load_bedgraph(args.reverse_bedgraph)

    origins = load_origins(args.origins_file)
    print(f"Loaded {len(origins)} origins from {args.origins_file}")
    generate_outputs(origins, fwd_bedgraph, rev_bedgraph, suffix="")

    named_origins = load_origins(args.origins_file, exclude_names={"null"})
    print(
        f"Loaded {len(named_origins)} named origins (excluding 'null') "
        f"from {args.origins_file}"
    )
    generate_outputs(named_origins, fwd_bedgraph, rev_bedgraph, suffix="_named")

    print(f"All outputs saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
