import os

SAMPLES=[i.split("_")[0] for i in os.listdir() if i.endswith("gz")] 

rule all:
	input:expand(["{sample}__forward.bedgraph","{sample}__reverse.bedgraph"],sample=SAMPLES)

rule cut_adapt_pair:
	input:
		read1="{sample}_R1_001.fastq.gz",
		read2="{sample}_R2_001.fastq.gz"
	output:
		read1="{sample}_R1_cut",
		read2="{sample}_R2_cut"
	shell:"cutadapt -a file:./bin/for.txt -A file:./bin/for.txt -f fastq --match-read-wildcards --quiet -m 15 -q 10 -o {output.read1} -p {output.read2} {input.read1} {input.read2}"

#rule cut_adapt_R1:
#	input:"{sample}_R1_001.fastq.gz"
#	output:"{sample}_R1_cut"
#	shell:"cutadapt -f fastq --match-read-wildcards -m 15 -q 10 -a file:./bin/for.txt {input} > {output}"

#rule cut_adapt_R2:
#        input:"{sample}_R2_001.fastq.gz"
#        output:"{sample}_R2_cut"
#        shell:"cutadapt -f fastq --match-read-wildcards -m 15 -q 10 -a file:./bin/for.txt {input} > {output}"

rule rectify:
	input:"{sample}_R1_cut","{sample}_R2_cut"
	output:"{sample}_R1_cut.paired.fastq","{sample}_R2_cut.paired.fastq"
	shell:"perl rectify_trimmed_pairs.pl {input}"

rule map_toward_oligo:
	input:"{sample}_R1_cut.paired.fastq"
	output:"{sample}_R1_cut.paired.fastq.unhit"
	shell:"bowtie -m1 -v2 --max R1.cutadapt.paired.oligo_max --un {output} ./oligo/oligo {input} dump.tmp"

rule extract_pair:
	input:i1="{sample}_R1_cut.paired.fastq.unhit",i2="{sample}_R2_cut.paired.fastq"
	output:"{sample}_R2_cut.paired.fastq.unhit"
	shell:"python extract.py {input.i1} {input.i2} {output}"

rule align:
	input:i1="{sample}_R1_cut.paired.fastq.unhit",i2="{sample}_R2_cut.paired.fastq.unhit"
	output:"{sample}_pair"
	#input:"Pol2MG.rnh201.1b.GSM1521150.SRR1609186_1122"
	shell:"bowtie -m1 -v2 -p8 -X2000  /data4/clausenLab_repository/programs/fastq_pipeline_genomes/yeast/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BowtieIndex/genome -1 {input.i1} -2 {input.i2} {output}"

rule extract_mate1:
	input:"{sample}_pair"
	output:"{sample}_"
	shell:"perl -ane 'print $_ if $F[0]=/\/1/;'  {input} > {output}"

rule bedgraph:
	input:"{sample}_"
	output:"{sample}__forward.bedgraph","{sample}__reverse.bedgraph"
	shell:"perl bowtie2bedgraph_new.pl {input} {input}"
