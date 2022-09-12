import os
import pandas as pd

if not os.path.exists('./masked'):
    os.makedirs('./masked')


def mask_file():
    files=[i for i in os.listdir("./bg/test_bg") if i.endswith("bedgraph")]
    #files=["KCl-uniq.human.SRR5035928_forward.bedgraph"]
    cut_file="./human_hincii_forward.txt"
    print (files)
    print ("Cut sites-file",cut_file)
    path="./bg/test_bg/"
    for file in files:
        print ("Removing cut sites from",file)
        bed_df=pd.read_csv(path+file,sep="\t",skiprows=1,header=None)
        cut_df=pd.read_csv(cut_file,sep="\t",header=None)
        cut_df[1]=cut_df[1]-5
        merged_df=bed_df.copy()
        for i in range(10):
            cut_df[1]=cut_df[1]+1
            merged_df = merged_df.merge(cut_df,indicator='i', how='outer').query('i == "left_only"').drop('i', axis=1)
        print ("I have removed",bed_df.shape[0]-merged_df.shape[0],"sites")
        print ("I have removed",bed_df[3].sum()-merged_df[3].sum(),"reads")
        bed_df.to_csv("./masked/cut_removed"+file,index=None)

mask_file()
