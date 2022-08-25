
# coding: utf-8

# In[60]:

import pybedtools
import os
import pandas as pd
import shutil
import numpy as np

import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
#import matplotlib as plt
sns.set_theme()


# In[41]:

def relative_ratios(df):
    df.drop(["Unnamed: 0"],axis=1,inplace=True)
    df.set_index(bases['file'],inplace=True)
    df.drop(["file"],axis=1,inplace=True)
    df['total']=bases.sum(axis=1)
    for i in ['A','G','C','T']:
        df["%"+i]=bases[i]/bases['total']
    return (df)


# In[42]:

def filter_bases(df):
    new_data=pd.DataFrame()
    roi={1:[(300,5746),(5848,7850),(7860,16199)]}
    for chrm,positions in roi.items():
        
        for pos in positions:
            
            roi=df[df[chrm].between(pos[0],pos[1])]
            #print (roi.head())
            new_data=new_data.append(roi)
    return (new_data)


# In[43]:

def mapping_to_genotype(df):
    names=pd.read_csv("SraRunTable-3.txt",sep=",")
    #print (names.head())
    print (names.columns)
    names['Genotype']=names['Genotype']+names['restriction_digest']
    #names['Genotype2'] = names.Genotype.str.cat(names['restriction digest'], sep='-')
    convert=names[['Run','Genotype']]
    #print (convert.head())
    convert=dict(zip(convert["Run"], convert["Genotype"]))
    df['Genotype']=df['SRR'].map(convert)
    #df=df.set_index(['Genotype'])
    return df


# In[44]:

def find_bases_forward(files):
    fasta = pybedtools.BedTool('/data4/clausenLab_repository/programs/fastq_pipeline_genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
    print ("There are",len(files),"libraries")
    bases=pd.DataFrame()
    for file in files:
        
        path="/home/anders/get_base/bg/"
        data=pd.read_csv(path+file,sep="\t",header=None,skiprows=1)
        data = data.rename({3: file}, axis=1)
        data[1]=data[1]-1
        data=data[data[1]>0]
        data=data[data[0].str.match('chrM')].reset_index(drop=True)
        data.to_csv(path+"skip"+file,sep="\t",index=False,header=None)
        C=pybedtools.BedTool(path+"skip"+file)
        C.sequence(fi=fasta,bedOut=True,fo="fito.tmp",name=False)
        df_base=pd.read_csv("fito.tmp",header=None)
        data.insert(4,'base',df_base[0])
        data=filter_bases(data)
        b1=data.drop([1,2],axis=1).groupby([data.iloc[:,4].str.upper()]).sum().T
        b1['total']=b1.sum(axis=1)
        for i in ['A','G','C','T']:
            b1["%"+i]=b1[i]/b1['total']
        bases=bases.append(b1,sort=True)
    
    bases['SRR']=bases.index.str[15:25]
    bases=mapping_to_genotype(bases)
    bases.to_csv("bases_forward.csv",sep="\t",index=None)
    return (bases)


# In[45]:

def find_bases_reverse(files):
    fasta = pybedtools.BedTool('/data4/clausenLab_repository/programs/fastq_pipeline_genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
    #files=[i for i in os.listdir("./bg") if i.endswith("reverse.bedgraph") and i.startswith("KCl-u") and int(i[24])%2!=0]
    print ("There are",len(files),"libraries")
    bases=pd.DataFrame()
    for file in files:
        print(file)
        path="/home/anders/get_base/bg/"
        data=pd.read_csv(path+file,sep="\t",header=None,skiprows=1)
        data = data.rename({3: file}, axis=1)
        
        data[1]=data[1]-1
        data=data[data[1]>0]
        data=data[data[0].str.match('chrM')].reset_index(drop=True)
        data.to_csv(path+"skip"+file,sep="\t",index=False,header=None)
        C=pybedtools.BedTool(path+"skip"+file)
        C.sequence(fi=fasta,bedOut=True,fo="fito.tmp",name=False,s=False)
        os.remove(path+"skip"+file)
        df_base=pd.read_csv("fito.tmp",header=None)
        data.insert(4,'base',df_base[0])
        data=filter_bases(data)
        b1=data.drop([1,2],axis=1).groupby([data.iloc[:,4].str.upper()]).sum().T
        b1['total']=b1.sum(axis=1)
        b1=b1.rename(columns={"A":"rU","G":"rC","C":"rG","T":"rA"})
        for i in ['rA','rG','rC','rU']:
            b1["%"+i]=b1[i]/b1['total']
        bases=bases.append(b1,sort=True)
    
    bases['SRR']=bases.index.str[15:25]
    bases=mapping_to_genotype(bases)
    bases.to_csv("bases_reverse.csv",sep="\t",index=None)
    return (bases)


# In[46]:

files_forward=[i for i in os.listdir("./bg/KOH") if i.endswith("forward.bedgraph")]

bases_forward=find_bases_forward(files_forward)


# In[47]:

files_reverse=[i for i in os.listdir("./bg/KOH") if i.endswith("reverse.bedgraph") ]
bases_reverse=find_bases_reverse(files_reverse)


# In[48]:

#total_df=bases_forward+bases_reverse
#total_df=pd.merge(bases_forward,bases_reverse,left_index=True, right_index=True)
total_df=pd.merge(bases_forward,bases_reverse,on='SRR')


# In[49]:

total_df.set_index('Genotype_x',inplace=True)


# In[50]:

total_df.to_csv("final.csv",sep="\t",index=None)


# In[51]:

total_df=pd.read_csv("final.csv",sep="\t")
total_df=total_df.set_index('Genotype_y')


# In[52]:

unscaled=total_df.iloc[:,np.r_[0:4,11:15]]


# In[53]:

with open('chrM.fa') as f:
    text=f.read()
    for i in ['A','T','G','C']:
        print (i,text.count(i))

col=unscaled.columns
bases=[4990,5068,2095,3997,3997,2095,5068,4990]
print (col)


for i in range(len(col)):
    unscaled['s'+col[i]]=unscaled[col[i]]/(bases[i]/16568)


# In[58]:

#sns.clustermap(unscaled[:,np.r_[0:4,11:15]])
sns.clustermap(unscaled.iloc[:,8:16], vmin=0, vmax=2)


# In[61]:

plt.savefig('output.png')


# In[ ]:



