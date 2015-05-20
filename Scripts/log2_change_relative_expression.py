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
for hum_relative_RPKM in ["Human_Muscle_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Human_Muscle_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["Chimp_Muscle_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Chimp_Muscle_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["Rhesus_Monkey_Muscle_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Rhesus_Monkey_Muscle_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for rhesus monkey
data_table[["Human_Muscle_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Muscle_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Muscle_Relative_Gene_Expression_Difference(RPKM)"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/muscle_log2_relative_expression_with_gene_names_RPKM.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Kidney tissue RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["Human_Kidney_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Human_Kidney_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["Chimp_Kidney_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Chimp_Kidney_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["Rhesus_Monkey_Kidney_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Rhesus_Monkey_Kidney_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for rhesus monkey
data_table[["Human_Kidney_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Kidney_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Kidney_Relative_Gene_Expression_Difference(RPKM)"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/kidney_log2_relative_expression_with_gene_names_RPKM.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Cerebellar Cortex tissue RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["Chimp_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Chimp_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["Rhesus_Monkey_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Rhesus_Monkey_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for rhesus monkey
data_table[["Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/cerebellar_cortex_log2_relative_expression_with_gene_names_RPKM.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Prefrontal Cortex RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["Human_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Human_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["Chimp_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Chimp_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["Rhesus_Monkey_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Rhesus_Monkey_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for rhesus monkey
data_table[["Human_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/prefrontal_cortex_log2_relative_expression_with_gene_names_RPKM.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Primary Visual Cortex RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["Human_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Human_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["Chimp_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Chimp_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["Rhesus_Monkey_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: np.log2((r["Rhesus_Monkey_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)"]) + 0.1), axis=1)
#calculate relative expression for rhesus monkey
data_table[["Human_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/primary_visual_cortex_log2_relative_expression_with_gene_names_RPKM.tsv", sep="\t", header=True)


data_table[["Human_Muscle_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Muscle_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Muscle_Relative_Gene_Expression_Difference(RPKM)", "Human_Kidney_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Kidney_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Kidney_Relative_Gene_Expression_Difference(RPKM)", "Human_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Cerebellar_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Human_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Prefrontal_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Human_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Chimp_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)", "Rhesus_Monkey_Primary_Visual_Cortex_Relative_Gene_Expression_Difference(RPKM)"]].to_csv("/Users/cmdb/Desktop/Rotation_4_-_Taylor_Lab/Data/Programs/all_species_tissue_log2_relative_expression_with_gene_names_RPKM.tsv", sep="\t", header=True)
