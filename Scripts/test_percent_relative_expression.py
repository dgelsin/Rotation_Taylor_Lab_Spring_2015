#!/usr/bin/env python

import sys
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 

f = sys.argv[1]

data_table=DataFrame.from_csv(f, header=0, sep="\t")

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Muscle tissue RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]*100)), axis=1)
#calculate relative expression for human
#for hum_relative_RPKM in ["Chimp_Muscle_Relative_Gene_Expression_Difference(RPKM)"]:
#	col_name = hum_relative_RPKM
#	data_table[col_name] = data_table.apply(lambda r: (((r["Chimp_Muscle_Average_RPKM"])/((r["Human_Muscle_Average_RPKM"] + r["Chimp_Muscle_Average_RPKM"] + r["Rhesus_Monkey_Muscle_Average_RPKM"])/3))*100), axis=1)
#calculate relative expression for chimp
#for hum_relative_RPKM in ["Rhesus_Monkey_Muscle_Relative_Gene_Expression_Difference(RPKM)"]:
#	col_name = hum_relative_RPKM
#	data_table[col_name] = data_table.apply(lambda r: (((r["Rhesus_Monkey_Muscle_Average_RPKM"])/((r["Human_Muscle_Average_RPKM"] + r["Chimp_Muscle_Average_RPKM"] + r["Rhesus_Monkey_Muscle_Average_RPKM"])/3))*100), axis=1)
#calculate relative expression for rhesus monkey
data_table[["Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/muscle_percent_relative_expression_with_gene_names_RPKM.tsv", sep="\t", header=True)
