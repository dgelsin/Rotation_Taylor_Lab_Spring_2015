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

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#total RPKM abundance expression
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

muscle_fig = data_table.boxplot(widths=0.15, return_type='axes', vert=False, whis=100, bootstrap=1000, patch_artist=True, fontsize=14, figsize=[100,100], positions=[15,14,13,12,11,10,9,8,7,6,5,4,3,2,1], grid=False, whiskerprops=whiskerprops, boxprops=boxprops, medianprops=medianprops)
#yaxis
plt.xlabel("Expression Abundance per Gene (RPKM)", fontsize=18, ha='center') #make a x axis label
muscle_fig.set_yticklabels(['hM', 'cM', 'rM','hK', 'cK', 'rK','hC', 'cC', 'rC','hP', 'cP', 'rP','hV', 'cV', 'rV'], fontsize=16)# add axis (tick) labels
#colors = ['b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g']
#for ytick, color in zip(muscle_fig.get_yticklabels(), colors):
#	ytick.set_color(color)
#set colors for y axis
get_fig = muscle_fig.get_figure()
get_fig.savefig('/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/Tissue_Specific_RNA-seq_data_sets/Plots/Boxplots/y-axis_all_avg-abundance_expression_distribution_boxplot_novert.png',bbox_inches='tight',dpi=300)


#muscle_fig = data_table.boxplot(widths=0.15, return_type='axes', vert=True, whis=100, bootstrap=1000, patch_artist=True, fontsize=14, rot=45, figsize=[100,100], positions=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], grid=False, whiskerprops=whiskerprops, boxprops=boxprops, medianprops=medianprops)
#xaxis
#plt.ylabel("Expression Abundance per Gene (RPKM)", fontsize=18, ha='center') #make a y axis label
#muscle_fig.set_xticklabels(['hM', 'cM', 'rM','hK', 'cK', 'rK','hC', 'cC', 'rC','hP', 'cP', 'rP','hV', 'cV', 'rV'], fontsize=16)# add axis (tick) labels
#colors = ['b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g']
#for xtick, color in zip(muscle_fig.get_xticklabels(), colors):
#	xtick.set_color(color)
#set colors for x axis

#get_fig = muscle_fig.get_figure()
#get_fig.savefig('/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/Tissue_Specific_RNA-seq_data_sets/Plots/Boxplots/x-axis_all_avg-abundance_expression_distribution_boxplot_vert.png',bbox_inches='tight',dpi=300)



""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#difference expression
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#muscle_fig = data_table.boxplot(widths=0.15, return_type='axes', vert=True, whis=100, bootstrap=1000, patch_artist=True, fontsize=14, rot=45, figsize=[100,100], positions=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], grid=False, whiskerprops=whiskerprops, boxprops=boxprops, medianprops=medianprops)
#xaxis
#plt.ylabel("log2 Expression Difference per Gene (RPKM)", fontsize=18, ha='center') #make a y axis label
#muscle_fig.set_xticklabels(['hum/chimp', 'hum/rhesus','rhesus/chimp', 'hum/chimp', 'hum/rhesus','rhesus/chimp', 'hum/chimp', 'hum/rhesus','rhesus/chimp', 'hum/chimp', 'hum/rhesus','rhesus/chimp', 'hum/chimp', 'hum/rhesus','rhesus/chimp'], fontsize=10)# add axis (tick) labels
#colors = ['b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g']
#for xtick, color in zip(muscle_fig.get_xticklabels(), colors):
#	xtick.set_color(color)
#set colors for x axis
#get_fig = muscle_fig.get_figure()
#get_fig.savefig('/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/Tissue_Specific_RNA-seq_data_sets/Plots/Boxplots/x-axis_all_avg-difference_log2_expression_distribution_boxplot_vert.png',bbox_inches='tight',dpi=300)


#muscle_fig = data_table.boxplot(widths=0.15, return_type='axes', vert=False, whis=100, bootstrap=1000, patch_artist=True, fontsize=14, figsize=[100,100], positions=[15,14,13,12,11,10,9,8,7,6,5,4,3,2,1], grid=False, whiskerprops=whiskerprops, boxprops=boxprops, medianprops=medianprops)
#yaxis
#plt.xlabel("log2 Expression Difference per Gene (RPKM)", fontsize=18, ha='center') #make a x axis label
#muscle_fig.set_yticklabels(['hum/chimp', 'hum/rhesus','rhesus/chimp', 'hum/chimp', 'hum/rhesus','rhesus/chimp', 'hum/chimp', 'hum/rhesus','rhesus/chimp', 'hum/chimp', 'hum/rhesus','rhesus/chimp', 'hum/chimp', 'hum/rhesus','rhesus/chimp'], fontsize=10)# add axis (tick) labels
#colors = ['b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g']
#for ytick, color in zip(muscle_fig.get_yticklabels(), colors):
#	ytick.set_color(color)
#set colors for y axis
#get_fig = muscle_fig.get_figure()
#get_fig.savefig('/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/Tissue_Specific_RNA-seq_data_sets/Plots/Boxplots/y-axis_all_avg-difference_log2_expression_distribution_boxplot_novert.png',bbox_inches='tight',dpi=300)



#pylab.xticks([4], ['Average Gene Expresion per Gene (RPKM)']) # label each box plot "column"
#muscle_fig.xaxis.set_ticklabels(['m_hum/chimp', 'm_hum/rhesus','m_rhesus/chimp', 'k_hum/chimp', 'k_hum/rhesus','k_rhesus/chimp', 'c_hum/chimp', 'c_hum/rhesus','c_rhesus/chimp', 'p_hum/chimp', 'p_hum/rhesus','p_rhesus/chimp', 'v_hum/chimp', 'v_hum/rhesus','v_rhesus/chimp'], fontsize=10)# add axis (tick) labels
#muscle_fig.yaxis.set_ticklabels(['m_hum/chimp', 'm_hum/rhesus','m_rhesus/chimp', 'k_hum/chimp', 'k_hum/rhesus','k_rhesus/chimp', 'c_hum/chimp', 'c_hum/rhesus','c_rhesus/chimp', 'p_hum/chimp', 'p_hum/rhesus','p_rhesus/chimp', 'v_hum/chimp', 'v_hum/rhesus','v_rhesus/chimp'], fontsize=10)# add axis (tick) labels


#plt.title("Average Gene Expression Distribution per Species; Tissue Type: Muscle", fontsize = 26)

#plt.show()
