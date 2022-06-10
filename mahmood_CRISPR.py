import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

df=pd.read_csv('uniq.human.190216-MF-AR-58-human-i11-i157_reverse.bedgraph',sep="\t",header=None,skiprows=1)

chrms=df[0].drop_duplicates()



def get_files():
    files=os.listdir()
    bedgraphs=[str(i) for i in files if i.endswith('bedgraph')]
    return (bedgraphs)


def make_plots():

    for file in get_files():
        print ("this is first file")
        print (file)
        df=pd.read_csv(file,sep="\t",header=None,skiprows=1)

        df.rename(columns={df.columns[0]:'chromosome'},inplace=True)
        df.rename(columns={df.columns[0]:'chromosome'},inplace=True)

        df['position']=df[1]/1000
        print (df.head())
        with PdfPages(str(file)+'.pdf') as pdf:

            for chrm in chrms:

                df_sub=df.loc[df['chromosome']==chrm]
                df_sub.plot.scatter(x='position',y=3,title=chrm,figsize=(25,5))
                sns.scatterplot
                pdf.savefig()

make_plots()
