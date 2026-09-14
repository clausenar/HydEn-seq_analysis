import collections
import glob
import os
from pybedtools import BedTool
import pysam

# --- CONFIGURATION ---
REFERENCE_FASTA = "/Users/xranea/genome/sacCer3.fa"
INPUT_DIR = "/Users/xranea/bedgraphs"
OUTPUT_DIR = "/Users/xranea/bedgraphs/processed_results"
TOTALS_FILE = "base_count_totals.txt"
# ---------------------

# Complement dictionary for reverse strand conversion
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def get_complement(sequence):
    """Returns the complementary base sequence for the reverse strand."""
    return "".join([COMPLEMENT.get(base, "N") for base in sequence.upper()])


def process_single_strand(bg_path, fasta_path, strand):
    """Processes a single bedGraph file for a given strand and returns its text data and base counts."""
    bg = BedTool(bg_path)

    # Build a coordinate -> score lookup from the original bedGraph, since
    # bedtools' sequence extraction does not carry the score column through.
    scores = {}
    with open(bg_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("track"):
                continue
            chrom, start_str, end_str, score = line.split("\t")
            scores[f"{chrom}:{start_str}-{end_str}"] = score

    # Fetch sequences using bedtools
    bg_with_seq = bg.sequence(fi=fasta_path, tab=True)

    local_counts = collections.Counter()
    output_lines = []

    with open(bg_with_seq.seqfn, "r") as f:
        for line in f:
            coords, sequence = line.strip().split("\t")
            chrom, intervals = coords.split(":")
            start_str, end_str = intervals.split("-")
            start = int(start_str)
            score = scores[coords]

            sequence = sequence.upper()
            # If it's the reverse strand, flip the bases to their complements
            if strand == "-":
                sequence = get_complement(sequence)

            for offset, base in enumerate(sequence):
                exact_position = start + offset
                output_lines.append(
                    f"{chrom}\t{exact_position}\t{strand}\t{score}\t{base}\n"
                )
                local_counts[base] += 1

    return output_lines, local_counts


def batch_process_directory():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Find all bedGraph files in the input directory
    all_files = glob.glob(os.path.join(INPUT_DIR, "*.bedGraph")) + glob.glob(
        os.path.join(INPUT_DIR, "*.bedgraph")
    )

    # Group files by sample prefix to link forward/reverse strands
    # Expects naming like: sample_name_forward.bedGraph / sample_name_reverse.bedGraph
    samples = {}
    for f in all_files:
        filename = os.path.basename(f)
        if "forward" in filename.lower():
            prefix = (
                filename.lower().split("forward")[0].rstrip("_.").rstrip("-")
            )
            samples.setdefault(prefix, {})["+"] = f
        elif "reverse" in filename.lower() or "rev" in filename.lower():
            prefix = (
                filename.lower()
                .split("reverse")[0]
                .split("rev")[0]
                .rstrip("_.")
                .rstrip("-")
            )
            samples.setdefault(prefix, {})["-"] = f
        else:
            # If no strand is specified in the name, default it to the forward strand
            prefix = os.path.splitext(filename)[0]
            samples.setdefault(prefix, {})["+"] = f

    if not samples:
        print(f"❌ No bedGraph files found in {INPUT_DIR}.")
        return

    print(f"📂 Found {len(samples)} distinct sample(s) to process.\n")

    all_sample_counts = {}

    # Process each identified sample
    for sample_name, strands in samples.items():
        print(f"🧬 Processing Sample: {sample_name}")
        global_counts = collections.Counter()
        out_file_path = os.path.join(OUTPUT_DIR, f"{sample_name}_expanded.txt")

        with open(out_file_path, "w") as out_file:
            out_file.write(
                "Chromosome\tPosition\tStrand\tBedGraph_Value\tBase\n"
            )

            # Process forward Strand if available
            if "+" in strands:
                print("  (+) Processing forward strand...")
                lines, counts = process_single_strand(
                    strands["+"], REFERENCE_FASTA, "+"
                )
                out_file.writelines(lines)
                global_counts.update(counts)

            # Process reverse Strand if available (will trigger complement logic)
            if "-" in strands:
                print("  (-) Processing reverse strand (complementing bases)...")
                lines, counts = process_single_strand(
                    strands["-"], REFERENCE_FASTA, "-"
                )
                out_file.writelines(lines)
                global_counts.update(counts)

        # Print the summary metrics for this specific sample
        print(f"📊 --- TOTAL BASE COUNTS FOR {sample_name.upper()} ---")
        for base, count in sorted(global_counts.items()):
            print(f"  Base {base}: {count:,}")
        print(f"💾 Per-position details saved to: {out_file_path}\n")

        all_sample_counts[sample_name] = global_counts

    # Write a combined totals summary for all samples
    totals_path = os.path.join(OUTPUT_DIR, TOTALS_FILE)
    with open(totals_path, "w") as totals_file:
        totals_file.write("Sample\tBase\tCount\n")
        for sample_name, counts in all_sample_counts.items():
            for base, count in sorted(counts.items()):
                totals_file.write(f"{sample_name}\t{base}\t{count}\n")
    print(f"💾 Base count totals saved to: {totals_path}\n")


if __name__ == "__main__":
    batch_process_directory()
    print("✅ Batch processing completed successfully!")
