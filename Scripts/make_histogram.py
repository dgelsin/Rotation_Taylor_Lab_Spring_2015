#!/usr/bin/env python


import sys
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 
from pylab import show

f = sys.argv[1]
data_table=DataFrame.from_csv(f, header=0, sep="\t")

plt.figure();
muscle_fig = data_table.plot(kind='bar', stacked=True)
get_fig = muscle_fig.get_figure()
get_fig.savefig('Histogram.png') 
#hum_greater_than_zero = data_table[data_table.Human_Muscle_Average_RPKM > 0.1]

#range = data_table.apply(lambda x: x.max() - x.min())
#print range

#maxi = data_table.max()
#print maxi

#mini = data_table.min()
#print mini
#data_table.plot(kind='bar')
#plt.title("Gaussian Histogram")
#plt.xlabel("Genes")
#plt.ylabel("Gene Expression (RPKM)")
#plt.show()

#DataFrame.plot(data_table, x=None, y=None, kind='bar')
#plt.show()