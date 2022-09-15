samples=["Pol2MG.rnh201.1b.GSM1521150.SRR1609186"]

rule all:
	input:"output_bedtools_bg_plus_strand.txt"


rule cut_adapt_pair:
	input:
		read1="data/samples/Pol2MG.rnh201.2a.GSM1521151.SRR1609188_R1_001.fastq.gz",
		read2="data/samples/Pol2MG.rnh201.2a.GSM1521151.SRR1609188_R2_001.fastq.gz"
	output:
		read1="SRR1609188_R1_001001.out",
		read2="SRR1609188_R2_001_001.out"
	shell:"cutadapt -a file:./bin/for.txt -A file:./bin/for.txt -f fastq --match-read-wildcards --quiet -m 15 -q 10 -o {output.read1} -p {output.read2} {input.read1} {input.read2}"

rule rectify:
	input:"SRR1609188_R1_001001.out","SRR1609188_R2_001_001.out"
	output:"SRR1609188_R1_001001.out.paired.fastq","SRR1609188_R2_001_001.out.paired.fastq"
	shell:"perl rectify_trimmed_pairs.pl {input}"

rule map_toward_oligo:
	input:"SRR1609188_R1_001001.out.paired.fastq"
	output:"R1.cutadapt.paired_1.oligo_unhit"
	shell:"bowtie -m1 -v2 --max R1.cutadapt.paired.oligo_max --un {output} ./oligo/oligo -1 {input.read1} -2 {input.read2} R1.cutadapt.paired.oligo_unhit_0913_V6.out"

rule extract_pair:
	input:i1="R1.cutadapt.paired.oligo_unhit",i2="SRR1609188_R1_001001.out.paired.fastq"
	output:"unhit_oligo-KOH_6"
	shell:"python extract.py {input.i1} {input.i2} {output}"

rule align:
	input:i1="R1.cutadapt.paired_1.oligo_unhit",i2="unhit_oligo-KOH_6"
	output:"sample1152"
	shell:"bowtie -m1 -v2 -p8 -X2000  /data4/clausenLab_repository/programs/fastq_pipeline_genomes/yeast/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BowtieIndex/genome -1 {input.i1} -2 {input.i2} {output}"

rule extract_mate1:
	input:"sample1152"
	output:"mate1.out.txt"
	shell:"perl -ane 'print $_ if $F[0]=/\/1/;'  {input} > {output}"

#rule realign_unhit_to_organism:
#	input:
#	output:
#	shell:"bowtie $m -v2 -p8 $index{$org} $tmp/"$sample"_1.unhit $tmp/"$sample"_1.unhit.out 2>>  $stats.3.tmp"

rule sam2bam:
	input:"sample_sam.pair"
	output:"sample_sam.pair.bam"
	shell:"samtools view -S {input} > {output}"

rule sort:
	input:"sample_sam.pair.bam"
	output:"sample_sam.pair.bam.sorted"
	shell:"samtools sort {input} -o {output}"

rule bedtools:
	input:"sample_sam.pair.bam.sorted"
	output:"output_bedtools_bg_plus_strand.txt"
	shell:"bedtools genomecov -d -5 -ibam {input} -g genome.fa > {output}"
