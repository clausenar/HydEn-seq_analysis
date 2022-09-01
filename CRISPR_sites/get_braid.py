import pandas as pd
import numpy as np
import os
import math

from pathlib import Path
Path("./img").mkdir(parents=True, exist_ok=True)


bins=2000
bins2=np.arange(0,200000,bins)

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

path=os.getcwd()

print ("Starting analysis...")

def braid(f1,f2,f3,f4,f5,f6,f7,f8,bins,chrm,pos_begin,pos_end,origin,crisp):


    file1=pd.read_csv(path+"/"+f1,sep="\t",header=None,skiprows=1)
    file2=pd.read_csv(path+"/"+f2,sep="\t",header=None,skiprows=1)
    file3=pd.read_csv(path+"/"+f3,sep="\t",header=None,skiprows=1)
    file4=pd.read_csv(path+"/"+f4,sep="\t",header=None,skiprows=1)
    file5=pd.read_csv(path+"/"+f5,sep="\t",header=None,skiprows=1)
    file6=pd.read_csv(path+"/"+f6,sep="\t",header=None,skiprows=1)
    file7=pd.read_csv(path+"/"+f7,sep="\t",header=None,skiprows=1)
    file8=pd.read_csv(path+"/"+f8,sep="\t",header=None,skiprows=1)

    files2=[file1,file2,file3,file4,file5,file6,file7,file8]

    data=[]
    appended_data2 = []
    for file in files2:
        data=file[file.iloc[:,0]==chrm]

        chr_bins=np.arange(int(pos_begin),int(pos_end),bins)
        file['bin']=pd.cut(x=data.iloc[:,2],bins=chr_bins,labels=chr_bins[:-1])
        groups = file.groupby(['bin']).sum()
        appended_data2.append(groups)
        data2=pd.concat(appended_data2,axis=1)

    data2.drop([1,2], axis=1,inplace=True)
    data2.columns=['file1','file2','file3','file4','file5','file6','file7','file8']

    ###crispr
    data2['crisp_mut_ratio']=data2['file1']/data2['file2']
    #data2['mut_norm']=1/(1+1/data2['crisp_mut_ratio'])

    #crispr_background
    data2['background_ratio']=data2['file3']/data2['file4']
    #data2['background_ratio_norm']=1/(1+1/data2['background_ratio'])

    #data2['crisp_braid']=1/(1+1/data2['crisp_mut_ratio']/data2['crisp_mut_ratio'])

    ##non-transformed
    data2['ctrl_ratio']=data2['file5']/data2['file6']
    #data2['ctrl_norm']=1/(1+1/data2['ctrl_ratio'])

    #non-transformed_background
    data2['ctrl_background_ratio']=data2['file7']/data2['file8']
    #data2['ctrl_background_ratio_norm']=1/(1+1/data2['ctrl_background_ratio'])

    ###new
    data2['crisp_mut_ratio_background_ratio_sqrt']=(data2['crisp_mut_ratio']/data2['background_ratio']).pow(1./2)
    data2['ctrl_ratio_bctrl_background_ratio_sqrt']=(data2['ctrl_ratio']/data2['ctrl_background_ratio']).pow(1./2)

    data2['crisp_mut_ratio_background_ratio_sqrt_norm']=1/(1+1/data2['crisp_mut_ratio_background_ratio_sqrt'])
    data2['ctrl_ratio_bctrl_background_ratio_sqrt_norm']=1/(1+1/data2['ctrl_ratio_bctrl_background_ratio_sqrt'])

    plt.clf()
    ax = plt.gca()

    data2.plot(kind='line',y='crisp_mut_ratio_background_ratio_sqrt_norm',ax=ax,figsize=(30,10))
    data2.plot(kind='line',y='ctrl_ratio_bctrl_background_ratio_sqrt_norm',ax=ax,figsize=(30,10))

    plt.axvline(x = int(origin)/2000, color = 'b', label = 'axvline - full height')
    plt.axvline(x = int(crisp)/2000, color = 'g', label = 'axvline - full height')
    file_name=chrm+origin+f1+f3+".png"
    print (file_name)

    plt.savefig(path+"/img/"+file_name)
    plt.show()

def read_files():
    f = open(path+"/master.csv", "r")
    f.readline()
    for line in f:
        #print (line)
        f1=line.split(';')[0]
        f2=f1.replace("forward.bedgraph","reverse.bedgraph")
        f3=line.split(';')[1]
        f4=f3.replace("forward.bedgraph","reverse.bedgraph")
        f5=line.split(';')[2]
        f6=f5.replace("forward.bedgraph","reverse.bedgraph").rstrip()
        f7=line.split(';')[3].rstrip()
        f8=f7.replace("forward.bedgraph","reverse.bedgraph")
        bins3=line.split(';')[4]
        chrm=line.split(';')[5]
        pos_begin=line.split(';')[6]
        pos_end=line.split(';')[7]
        origin=line.split(';')[8]
        crisp=line.split(';')[9]
        print (pos_end)

        braid (f1,f2,f3,f4,f5,f6,f7,f8,bins,chrm,pos_begin,pos_end,origin,crisp)


read_files()
