import os

RAW_DIR = "/Users/xranea/raw"
OLIGOS = "/Users/xranea/git/HydEn-seq_analysis/bin/oligo/oligo"
GENOME = "/Users/xranea/genome/sacCer3"
GENOME_FAI = GENOME + ".fa.fai"

# count_bases.py and origin_metaplot.py need pybedtools/matplotlib, which
# live in bio_env, not the hyden_pipeline env this Snakefile itself runs in.
# Must go through `conda run` (not a direct path to bio_env's python binary):
# pybedtools shells out to the `bedtools` CLI, which needs bio_env's PATH set
# up, and a direct binary path skips that activation step entirely.
BIO_ENV_PYTHON = "conda run -n bio_env python"

# Denna mapp kommer nu att hålla ALLA genererade filer (både temporära och slutliga)
OUT_DIR = "/Users/xranea/bedgraphs"

# FIX: Lägg till [0] efter split för att få ut strängen, vilket gör den "hashable" för set()
# Accept both plain and gzipped fastq; cutadapt/bowtie read .fastq.gz natively.
SAMPLES = list(set([
	i.split("_end")[0]
	for i in os.listdir(RAW_DIR)
	if i.endswith("fastq") or i.endswith("fastq.gz")
]))

print(SAMPLES)


def raw_fastq(sample, end):
	"""Resolves a raw input fastq path, preferring the gzipped form if both exist."""
	for ext in ("fastq.gz", "fastq"):
		path = os.path.join(RAW_DIR, f"{sample}_{end}.{ext}")
		if os.path.exists(path):
			return path
	raise FileNotFoundError(f"No {end} fastq found for sample {sample} in {RAW_DIR}")

rule all:
	input:
		expand([os.path.join(OUT_DIR, "{sample}__forward.bedgraph"), os.path.join(OUT_DIR, "{sample}__reverse.bedgraph")], sample=SAMPLES),
		os.path.join(OUT_DIR, "processed_results", "base_count_totals.txt"),
		expand(os.path.join(OUT_DIR, "processed_results", "{sample}", "origin_metaplot.png"), sample=SAMPLES)


rule cut_adapt_pair:
	input:
		read1=lambda wc: raw_fastq(wc.sample, "end1"),
		read2=lambda wc: raw_fastq(wc.sample, "end2")
	output:
		read1=temp(os.path.join(OUT_DIR, "{sample}_R1_cut")),
		read2=temp(os.path.join(OUT_DIR, "{sample}_R2_cut"))
	threads: 8
	# The old pipeline had a rectify_trimmed_pairs.pl step here that (a)
	# re-synced R1/R2 by read ID - redundant, cutadapt's paired mode already
	# keeps them in lockstep (verified by diffing read IDs) - and (b) blindly
	# trimmed the last base off every read. That trim is unneeded: the mapped
	# ribonucleotide position comes from mate1's 5' end (see the bedgraph
	# rule), which a 3'-end trim doesn't affect.
	shell:"cutadapt -j {threads} -a file:./bin/for.txt -A file:./bin/for.txt --match-read-wildcards --quiet -m 15 -q 10 -o {output.read1} -p {output.read2} {input.read1} {input.read2}"

rule map_toward_oligo:
	input:
		os.path.join(OUT_DIR, "{sample}_R1_cut")
	output:
		temp(os.path.join(OUT_DIR, "{sample}_R1_cut.unhit"))
	params:
		oligo_max=os.path.join(OUT_DIR, "{sample}_R1.cutadapt.paired.oligo_max"),
		dump=os.path.join(OUT_DIR, "{sample}_dump.tmp")
	# FIX: Lagt till {OUT_DIR}/ framför oligo_max-filen så att den hamnar utanför din repo
	# oligo_max/dump are wildcarded by {sample} so concurrent samples don't race on the same file.
	shell: "bowtie -m1 -v2 --max {params.oligo_max} --un {output} -x {OLIGOS} {input} {params.dump}"


rule extract_pair:
	input:
		i1=os.path.join(OUT_DIR, "{sample}_R1_cut.unhit"),
		i2=os.path.join(OUT_DIR, "{sample}_R2_cut")
	output:
		temp(os.path.join(OUT_DIR, "{sample}_R2_cut.unhit"))
	params:
		tmpdir=os.path.join(OUT_DIR, "{sample}_pairtmp")
	# seqkit pair's default output naming breaks on our dotted sample names
	# (it mis-detects an "extension" and inserts ".paired" mid-name), so we
	# give it an explicit scratch dir (-O keeps filenames unchanged there)
	# and move the R2 result out to our declared output path.
	shell:
		"seqkit pair -1 {input.i1} -2 {input.i2} -O {params.tmpdir} -f && "
		"mv {params.tmpdir}/$(basename {input.i2}) {output} && "
		"rm -rf {params.tmpdir}"

rule align:
	input:
		i1=os.path.join(OUT_DIR, "{sample}_R1_cut.unhit"),
		i2=os.path.join(OUT_DIR, "{sample}_R2_cut.unhit")
	output:
		temp(os.path.join(OUT_DIR, "{sample}_pair.sam"))
	threads: 8
	shell: "bowtie -S -m1 -v2 -p{threads} -X2000 -x {GENOME} -1 {input.i1} -2 {input.i2} {output}"

rule sort_bam:
	input:
		os.path.join(OUT_DIR, "{sample}_pair.sam")
	output:
		temp(os.path.join(OUT_DIR, "{sample}_pair.sorted.bam"))
	shell: "samtools sort -o {output} {input}"

rule extract_mate1:
	input:
		os.path.join(OUT_DIR, "{sample}_pair.sorted.bam")
	output:
		os.path.join(OUT_DIR, "{sample}_mate1.bam")
	# -f 64: first-in-pair (mate1) reads only, -F 4: drop unmapped reads
	shell: "samtools view -b -f 64 -F 4 {input} > {output}"

rule bedgraph:
	input:
		os.path.join(OUT_DIR, "{sample}_mate1.bam")
	output:
		fw = os.path.join(OUT_DIR, "{sample}__forward.bedgraph"),
		rv = os.path.join(OUT_DIR, "{sample}__reverse.bedgraph")
	# -5: each read's 5' base, -strand: split by mapped strand.
	# The incorporated ribonucleotide sits one base 5' of the read start, so
	# shift forward-strand hits left by 1bp and reverse-strand hits right by 1bp.
	shell:
		"bedtools genomecov -ibam {input} -bg -5 -strand + | "
		"awk -F'\\t' 'BEGIN{{OFS=\"\\t\"}} {{s=$2-1; e=$3-1; if(s>=0){{print $1,s,e,$4}}}}' > {output.fw} && "
		"bedtools genomecov -ibam {input} -bg -5 -strand - | "
		"awk -F'\\t' 'BEGIN{{OFS=\"\\t\"}} "
		"NR==FNR{{len[$1]=$2; next}} "
		"{{e=$3+1; if(e<=len[$1]){{print $1,$2+1,e,$4}}}}' {GENOME_FAI} - > {output.rv}"

rule count_bases:
	input:
		expand([os.path.join(OUT_DIR, "{sample}__forward.bedgraph"), os.path.join(OUT_DIR, "{sample}__reverse.bedgraph")], sample=SAMPLES)
	output:
		os.path.join(OUT_DIR, "processed_results", "base_count_totals.txt")
	# count_bases.py batch-processes every bedgraph pair it finds in its own
	# INPUT_DIR (which defaults to OUT_DIR), so it needs no per-sample args.
	shell: "{BIO_ENV_PYTHON} count_bases.py"

rule origin_metaplot:
	input:
		fw=os.path.join(OUT_DIR, "{sample}__forward.bedgraph"),
		rv=os.path.join(OUT_DIR, "{sample}__reverse.bedgraph")
	output:
		os.path.join(OUT_DIR, "processed_results", "{sample}", "origin_metaplot.png")
	params:
		outdir=os.path.join(OUT_DIR, "processed_results", "{sample}")
	shell:
		"{BIO_ENV_PYTHON} origin_metaplot.py "
		"--forward-bedgraph {input.fw} --reverse-bedgraph {input.rv} "
		"--output-dir {params.outdir}"
