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
#num_rows = len(data_table)

#print num_rows
#count the number of rows

#empty = data_table.apply(lambda col: pd.isnull(col))

#print empty
for hum_relative_RPKM in ["Relative_Expression_Difference_Hum_Chimp(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: r["Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"] - r["Chimp_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"], axis=1)

for hum_relative_RPKM in ["Relative_Expression_Difference_Hum_Rhesus_Monkey(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: r["Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"] - r["Rhesus_Monkey_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"], axis=1)

for hum_relative_RPKM in ["Relative_Expression_Difference_Chimp_Rhesus_Monkey(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: r["Chimp_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"] - r["Rhesus_Monkey_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"], axis=1)

data_table[["Relative_Expression_Difference_Hum_Chimp(RPKM)", "Relative_Expression_Difference_Hum_Rhesus_Monkey(RPKM)", "Relative_Expression_Difference_Chimp_Rhesus_Monkey(RPKM)"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/Cerebellar_Cortex_Expression_Difference.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Cerebellar Cortex tissue RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["log2_hum-chimp_difference"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2(r["Relative_Expression_Difference_Hum_Chimp(RPKM)"]), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["log2_hum-rhe_monkey_difference"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2(r["Relative_Expression_Difference_Hum_Rhesus_Monkey(RPKM)"]), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["log2_chimp-rhe_monkey_difference"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2(r["Relative_Expression_Difference_Chimp_Rhesus_Monkey(RPKM)"]), axis=1)
#calculate relative expression for rhesus monkey
data_table[["log2_hum-chimp_difference", "log2_hum-rhe_monkey_difference", "log2_chimp-rhe_monkey_difference"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/cerebellar_cortex_log2_relative_expression_with_gene_names_RPKM.tsv", sep="\t", header=True)
