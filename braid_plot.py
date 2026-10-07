"""braid_plot.py - genome-wide Watson/Crick strand ratio ("braid plot").

Bins the whole genome into fixed-size windows, sums each strand's
ribonucleotide hits per bin from the sample's forward (Watson, +) and
reverse (Crick, -) bedGraphs, and plots log2(Watson/Crick) as one panel per
chromosome, windowed to that chromosome's own size, with each origin from
the origins file marked as a vertical line. The ratio's sign alternates
across replication origins and termini, giving the characteristic "braid"
pattern this plot is named for. Matches origin_metaplot.py's own
Watson/Crick convention.
"""

import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- CONFIGURATION (defaults; overridable via CLI args, see main()) ---
# or200.txt ships alongside this script, so default to it there rather than
# an absolute path tied to one specific checkout location.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FORWARD_BEDGRAPH = "/Users/xranea/bedgraphs/Kunkel_Ribo-seq_Pol2MGrnh201.1b.1__forward.bedgraph"
REVERSE_BEDGRAPH = "/Users/xranea/bedgraphs/Kunkel_Ribo-seq_Pol2MGrnh201.1b.1__reverse.bedgraph"
GENOME_FAI = os.environ.get("HYDEN_GENOME", "/Users/xranea/genome/sacCer3") + ".fa.fai"
ORIGINS_FILE = os.environ.get("HYDEN_ORIGINS_FILE", os.path.join(_SCRIPT_DIR, "or200.txt"))
OUTPUT_DIR = os.path.join(os.environ.get("HYDEN_OUT_DIR", "/Users/xranea/bedgraphs"), "processed_results")

BIN_SIZE = 500   # bp per bin, genome-wide
# ---------------------

# or200.txt uses arabic chromosome numbers (chr1..chr16); the bedgraphs use
# roman numerals (chrI..chrXVI), matching the sacCer3 reference naming.
ROMAN = [
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
    "XI", "XII", "XIII", "XIV", "XV", "XVI",
]
CHROM_MAP = {f"chr{i+1}": f"chr{roman}" for i, roman in enumerate(ROMAN)}


def load_origins(path):
    """Reads the origin BED-like file into {chrom: [midpoint, ...]}."""
    by_chrom = {}
    with open(path, "r") as f:
        for line in f:
            if line.startswith("track") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[:3]
            chrom = CHROM_MAP.get(chrom, chrom)
            midpoint = (int(start) + int(end)) // 2
            by_chrom.setdefault(chrom, []).append(midpoint)
    return by_chrom


def load_chrom_sizes(fai_path):
    """Reads a .fa.fai index into {chrom: length}, in file order."""
    sizes = {}
    with open(fai_path, "r") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            sizes[fields[0]] = int(fields[1])
    return sizes


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


def bin_chrom_signal(intervals, chrom_len, bin_size):
    """Sums per-bin totals across a whole chromosome.

    Every interval in these bedGraphs is exactly 1bp wide (bedtools
    genomecov -bg -5 pileups), so each one falls entirely inside a single
    bin - no need for origin_metaplot's overlap-splitting logic here, just
    a bincount on each interval's start position.
    """
    n_bins = -(-chrom_len // bin_size)  # ceil division
    bins = np.zeros(n_bins)
    if not intervals:
        return bins
    starts = np.array([iv[0] for iv in intervals])
    scores = np.array([iv[2] for iv in intervals])
    bin_idx = starts // bin_size
    np.add.at(bins, bin_idx, scores)
    return bins


def compute_log2_ratio(watson, crick, pseudocount=1.0):
    """log2(Watson/Crick) per bin, with a pseudocount to avoid division by zero."""
    return np.log2((watson + pseudocount) / (crick + pseudocount))


def plot_braid(chrom_sizes, fwd_bedgraph, rev_bedgraph, origins, bin_size, out_path):
    chroms = list(chrom_sizes.keys())
    max_len = max(chrom_sizes.values())
    fig, axes = plt.subplots(len(chroms), 1, figsize=(16, 1.4 * len(chroms)), sharex=False)
    if len(chroms) == 1:
        axes = [axes]

    for ax, chrom in zip(axes, chroms):
        chrom_len = chrom_sizes[chrom]
        watson = bin_chrom_signal(fwd_bedgraph.get(chrom, []), chrom_len, bin_size)
        crick = bin_chrom_signal(rev_bedgraph.get(chrom, []), chrom_len, bin_size)
        ratio = compute_log2_ratio(watson, crick)
        positions = np.arange(len(ratio)) * bin_size

        ax.fill_between(positions, ratio, 0, where=(ratio >= 0), color="tab:red", linewidth=0)
        ax.fill_between(positions, ratio, 0, where=(ratio < 0), color="tab:blue", linewidth=0)
        ax.axhline(0, color="black", linewidth=0.5)
        for mid in origins.get(chrom, []):
            ax.axvline(mid, color="black", linewidth=0.5, alpha=0.5, linestyle="--")
        # Window this panel to its own chromosome's size, not a shared
        # genome-wide axis, so short chromosomes aren't mostly blank space.
        ax.set_xlim(0, chrom_len)
        ax.set_ylabel(chrom.replace("chr", ""), rotation=0, labelpad=20, va="center")
        ax.set_yticks([])

    axes[-1].set_xlabel("Position (bp)")
    fig.suptitle(f"Braid plot: log2(Watson / Crick) ribonucleotide ratio, {bin_size}bp bins")
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    # Narrow each panel's on-page width to its chromosome's share of the
    # longest chromosome, now that tight_layout has settled everyone's
    # position - so panel width visually reflects relative chromosome length
    # on top of each panel's own x-axis being windowed to its own size.
    for ax, chrom in zip(axes, chroms):
        pos = ax.get_position()
        width = pos.width * chrom_sizes[chrom] / max_len
        ax.set_position([pos.x0, pos.y0, width, pos.height])

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_matrix(chrom_sizes, fwd_bedgraph, rev_bedgraph, bin_size, out_path):
    """Writes one row per genome bin: chrom, bin_start, watson, crick, log2_ratio."""
    with open(out_path, "w") as f:
        f.write("chrom\tbin_start\twatson\tcrick\tlog2_watson_crick\n")
        for chrom, chrom_len in chrom_sizes.items():
            watson = bin_chrom_signal(fwd_bedgraph.get(chrom, []), chrom_len, bin_size)
            crick = bin_chrom_signal(rev_bedgraph.get(chrom, []), chrom_len, bin_size)
            ratio = compute_log2_ratio(watson, crick)
            for i in range(len(ratio)):
                f.write(f"{chrom}\t{i * bin_size}\t{watson[i]}\t{crick[i]}\t{ratio[i]:.4f}\n")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--forward-bedgraph", default=FORWARD_BEDGRAPH)
    parser.add_argument("--reverse-bedgraph", default=REVERSE_BEDGRAPH)
    parser.add_argument("--genome-fai", default=GENOME_FAI)
    parser.add_argument("--origins-file", default=ORIGINS_FILE)
    parser.add_argument("--output-dir", default=OUTPUT_DIR)
    parser.add_argument("--bin-size", type=int, default=BIN_SIZE)
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    chrom_sizes = load_chrom_sizes(args.genome_fai)
    fwd_bedgraph = load_bedgraph(args.forward_bedgraph)
    rev_bedgraph = load_bedgraph(args.reverse_bedgraph)
    origins = load_origins(args.origins_file)
    print(f"Loaded {sum(len(v) for v in origins.values())} origins from {args.origins_file}")

    plot_path = os.path.join(args.output_dir, "braid_plot.png")
    matrix_path = os.path.join(args.output_dir, "braid_matrix.tsv")

    plot_braid(chrom_sizes, fwd_bedgraph, rev_bedgraph, origins, args.bin_size, plot_path)
    save_matrix(chrom_sizes, fwd_bedgraph, rev_bedgraph, args.bin_size, matrix_path)

    print(f"Braid plot saved to: {plot_path}")
    print(f"Binned matrix saved to: {matrix_path}")


if __name__ == "__main__":
    main()
