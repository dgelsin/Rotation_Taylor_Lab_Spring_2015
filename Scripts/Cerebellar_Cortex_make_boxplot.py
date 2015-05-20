#!/usr/bin/env python


import sys
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 
import pylab
from pylab import show
import brewer2mpl

f = sys.argv[1]
data_table=DataFrame.from_csv(f, header=0, sep="\t")

plt.figure();
boxprops = dict(linewidth=3, color='red')
medianprops = dict(linestyle='solid', linewidth=2.5, color='yellow')
whiskerprops = dict(linestyle='solid', color='green')


cerebellar_fig = data_table.boxplot(widths=0.15, return_type='axes', vert=False, whis=100, bootstrap=1000, patch_artist=True, fontsize=14, rot=45, figsize=[100,100], positions=[3,2,1], grid=False, whiskerprops=whiskerprops, boxprops=boxprops, medianprops=medianprops)
plt.title("Average Gene Expression Distribution per Species; Tissue Type: Cerebellar Cortex", fontsize = 26, color='blue')

plt.xlabel("Average Expression per Gene (RPKM)", fontsize=18, ha='center') #make a y axis label

#pylab.xticks([4], ['Average Gene Expresion per Gene (RPKM)']) # label each box plot "column"
cerebellar_fig.yaxis.set_ticklabels(['human', 'chimp', 'rhesus monkey', 'mouse'], fontsize=16)# add axis (tick) labels
#pylab.yticks([1, 2, 3, 4], ['Human', 'Chimp', 'Rhesus Monkey', 'Mouse'])# same as above

#plt.show()
get_fig = cerebellar_fig.get_figure()
get_fig.savefig('/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/Tissue_Specific_RNA-seq_data_sets/Plots/Boxplots/1_cerebellar_cortex_boxplot.png',bbox_inches='tight',dpi=300) 