#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib as plt
import os

def make_figures():
    files=[i for i in os.listdir() if i.endswith("tsv")]
    for file in files:
        df=pd.read_csv(file,sep="\t")
        df['Position']=df['Position']/1e6
        ax=df.plot('Position',y="X1",xlim=(0,4.67),figsize=(10,2),ylim=(0.45,0.55),legend=False,yticks=(0.450,0.50,0.55))
        ax.fill_between(df.Position,0,df.X1)
        x_axis = ax.axes.get_xaxis()
        x_label = x_axis.get_label()
        x_label.set_visible(False)

        ax.figure.savefig(file+"plot.png")
        
    
make_figures()

