import os
import time
import sys
import csv
import re

human_fasta="/data4/clausenLab_repository/programs/fastq_pipeline_genomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
yeast_fasta="/data4/clausenLab_repository/programs/fastq_pipeline_genomes/yeast/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa"

fasta=sys.argv[1]
seqs=sys.argv[2]

if fasta=="human":
    fasta=human_fasta
    print ("using human genome")

if seqs=="hincii":
    seqs=['GTCAAC','GTCGAC','GTTAAC','GTTGAC']
    print ("Looking for hincII")

from collections import defaultdict
chromosomes=defaultdict(list)
a=time.time()
for line in open(fasta).readlines():
    line=line.strip()
    if line.startswith(">"):
        chrm=line.replace(">","").split()[0]
        #chromosomes2[chrm]=""
        print (chrm)
    else:
        #chromosomes[chrm].append(list(line.strip()))
        line=[line]
        chromosomes[chrm].extend(line)
        
for key,value in chromosomes.items():
    chromosomes[key]="".join(value)

b=time.time()
print (b-a)

def sites(seqs,offset):
    sites=defaultdict(list)
    
    for i in chromosomes:
        cut_sites=[]
        print (i)
        for m in re.finditer(r'|'.join(seqs),chromosomes[i]):
            cut_sites.extend([m.start()+offset])
        sites[i]=cut_sites
    return sites

def make_csv(sites,out):
    with open (out,"w") as f:
        writer = csv.writer(f, delimiter='\t')
        for i in sites:
            for pos in sites[i]:
                output=i,pos
                writer.writerow(output)
                
def total_sites(sites):
    total_sites=0
    for i in sites:
        total_sites+=len(sites[i])
    print ("The total number of sites was:",total_sites)

forward=sites(seqs,3)
make_csv(forward,sys.argv[1]+"_"+sys.argv[2]+"_forward.txt")

reverse=sites(seqs,3)
make_csv(reverse,sys.argv[1]+"_"+sys.argv[2]+"_reverse.txt")

total_sites(forward)

