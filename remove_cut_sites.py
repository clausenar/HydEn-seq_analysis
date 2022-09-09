
# coding: utf-8



import os
import pandas as pd
import csv

print (os.getcwd())

bed_file="KCl-uniq.human.SRR5035909_forward.bedgraph"
cut_file="forward_hincII_human_new.txt"
bed_df=pd.read_csv(bed_file,sep="\t",header=None,skiprows=1,nrows=5)
cut_df=pd.read_csv("../"+cut_file,sep="\t",header=None,skiprows=1,nrows=4)

def make_mask(file):
    output=open("mask"+file,"w+")
    with open("../"+file) as csv_file:
        data=csv.reader(csv_file,delimiter="\t")
        for row in data:
            for i in range(-5,6):
                chrm=row[0]
                pos=int(row[1])+i
                #line=chrm+","+str(pos)
                #line=",".join(line)
                line=row[0]+","+str(pos)
                #print (row[0],pos)
                #print(line)
                output.write(line)
    output.close()
                
make_mask(cut_file)

cut_file.head()


new_file=pd.merge(bed_file,cut_file)

mask=pd.read_csv("maskforward_hincII_human_new.txt",sep="\t")

no_cut = bed_df.merge(mask.drop_duplicates(), on=[0,1], 
                   how='left', indicator=True)

print (no_cut.shape)

