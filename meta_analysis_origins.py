#!/usr/bin/env python
# coding: utf-8




import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os





trans={
        "chr1":"chrI",
        "chr2":"chrII",
        "chr3":"chrIII",
        "chr4":"chrIV",
        "chr5":"chrV",
        "chr6":"chrVI",
        "chr7":"chrVII",
        "chr8":"chrVIII",
        "chr9":"chrIX",
        "chr10":"chrX",
        "chr11":"chrXI",
        "chr12":"chrXII",
        "chr13":"chrXII",
        "chr14":"chrXIV",
        "chr15":"chrXV",
        "chr16":"chrXVI"}





def make_meta(file):
    meta=pd.DataFrame(pd.interval_range(start=-2000, freq=50, end=2000))
    meta=pd.DataFrame()
    print (meta.shape)
    origins=open('saccer1_origins_liftover_saccer3_redone.tsv')
    
    for line in origins.readlines():
        chrm=line.split()[2]
        chrm=trans[chrm]
        pos=int(line.split()[3])
        print (chrm,pos)
        data=pd.read_csv(file,sep="\t",header=None,skiprows=1)
        interval_range = pd.interval_range(start=-2000+pos, freq=50, end=2000+pos)
        data['bin']=pd.cut(data[1],bins=interval_range)
        data=data[data[0]==chrm]
        #data.dropna()
        data['bin'].dropna(inplace=True)
        data=data.groupby(['bin']).sum()
        data.reset_index(drop=True,inplace=True)
        data=data.set_index(np.arange(-2000,2000,50))
        meta=pd.concat([meta,data[3]],axis=1)
        #print (meta)
        #data[3].plot()
        #plt.show()
    #sns.heatmap(data.T)
    #plt.savefig("sns"+file+".png")
    meta['total']=meta.sum(axis = 1, skipna = True)
    #data['total'].plot()
    #plt.show()
    return meta['total']



def total_counts(file):
    data=pd.read_csv(file,sep="\t",header=None,skiprows=1)
    return data[3].sum()

files=[i for i in os.listdir() if i.endswith("bedgraph")]
for file in files:
    print (file,total_counts(file))





def make_plots():
    files=[i for i in os.listdir() if i.endswith("forward.bedgraph")]
    files=['uniq.yeast.Pol2MG.rnh201.1b.GSM1521150.SRR1609185_forward.bedgraph']
    control="uniq.yeast.Rnh201.1.a.GSM1521143.SRR1609175_forward.bedgraph"
    control_forward=make_meta(control)
    control_reverse_file=control.replace("forward","reverse")
    control_reverse=make_meta(control_reverse_file)
    for file in files:
        forward_file=file
        reverse_file=file.replace("forward","reverse")
        forward=make_meta(forward_file)
        reverse=make_meta(reverse_file)
        full=pd.DataFrame()
        full=pd.concat([full,forward,reverse,control_forward,control_reverse],axis=1,names=["Top","Bottom"])
        full.columns = ['Top','Bottom','Control_Top','Control_Bottom']
        full['graded_top']=full['Top']/full['Control_Top']*total_counts(control)/total_counts(forward_file)
        full['graded_bottom']=full['Bottom']/full['Control_Bottom']*total_counts(control_reverse_file)/total_counts(reverse_file)
        
        #full=full.set_index(np.arange(-1950,2050,50))
        full[['Top','Bottom','Control_Top','Control_Bottom']].plot()
        full[['graded_top','graded_bottom']].plot()
        plt.title(file)
        plt.savefig(file+".png",facecolor='auto')
        plt.show()
    return full










total=make_plots()








