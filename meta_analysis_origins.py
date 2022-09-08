#!/usr/bin/env python
# coding: utf-8

# In[15]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os


# In[37]:


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


# In[61]:


def make_meta(file):
    meta=pd.DataFrame(pd.interval_range(start=-2000, freq=50, end=2000))
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
        data=data.groupby(['bin']).sum()
        data.reset_index(drop=True,inplace=True)
        data=data.set_index(np.arange(-1950,2050,50))
    sns.heatmap(data.T)
    plt.savefig(sns+file+".png")
    data['total']=data.sum(axis = 1, skipna = True)
    return data['total']


# In[62]:


def make_plots():
    files=[i for i in os.listdir() if i.endswith("forward.bedgraph")]
    for file in files:
        forward_file=file
        reverse_file=file.replace("forward","reverse")
        forward=make_meta(forward_file)
        reverse=make_meta(reverse_file)
        full=pd.DataFrame()
        full=pd.concat([full,forward,reverse],axis=1,names=["Top","Bottom"])
        full.columns = ['Top','Bottom']
        full=full.set_index(np.arange(-1950,2050,50))
        full.plot(use_index=True)
        plt.title(file)
        plt.savefig(file+".png",facecolor='auto')
        plt.show()


# In[ ]:


make_plots()


# In[514]:





# In[473]:





# In[474]:





# In[ ]:




