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
for hum_relative_RPKM in ["m_Hum-Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Muscle_Average_RPKM"] +0.1) - np.log2(r["Chimp_Muscle_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["m_Hum-Rhesus_Monkey"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Muscle_Average_RPKM"] +0.1) - np.log2(r["Rhesus_Monkey_Muscle_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["m_Rhesus_Monkey_Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Rhesus_Monkey_Muscle_Average_RPKM"]+0.1) - np.log2(r["Chimp_Muscle_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for rhesus monkey
data_table[["m_Hum-Chimp", "m_Hum-Rhesus_Monkey", "m_Rhesus_Monkey_Chimp"]].to_csv("muscle_relative_expression_with_gene_names_RPKM_difference_LOG2.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Kidney tissue RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["k_Hum-Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Kidney_Average_RPKM"]+0.1) - np.log2(r["Chimp_Kidney_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["k_Hum-Rhesus_Monkey"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Kidney_Average_RPKM"]+0.1) - np.log2(r["Rhesus_Monkey_Kidney_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["k_Rhesus_Monkey_Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Rhesus_Monkey_Kidney_Average_RPKM"]+0.1) - np.log2(r["Chimp_Kidney_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for rhesus monkey
data_table[["k_Hum-Chimp", "k_Hum-Rhesus_Monkey", "k_Rhesus_Monkey_Chimp"]].to_csv("kidney_relative_expression_with_gene_names_RPKM_difference_LOG2.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Cerebellar Cortex tissue RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["c_Hum-Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Cerebellar_Cortex_Average_RPKM"]+0.1) - np.log2(r["Chimp_Cerebellar_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["c_Hum-Rhesus_Monkey"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Cerebellar_Cortex_Average_RPKM"]+0.1) - np.log2(r["Rhesus_Monkey_Cerebellar_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["c_Rhesus_Monkey_Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Rhesus_Monkey_Cerebellar_Cortex_Average_RPKM"]+0.1) - np.log2(r["Chimp_Cerebellar_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for rhesus monkey
data_table[["c_Hum-Chimp", "c_Hum-Rhesus_Monkey", "c_Rhesus_Monkey_Chimp"]].to_csv("cerebellar_cortex_relative_expression_with_gene_names_RPKM_difference_LOG2.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Prefrontal Cortex RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["p_Hum-Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Prefrontal_Cortex_Average_RPKM"]+0.1) - np.log2(r["Chimp_Prefrontal_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["p_Hum-Rhesus_Monkey"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Prefrontal_Cortex_Average_RPKM"]+0.1) - np.log2(r["Rhesus_Monkey_Prefrontal_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["p_Rhesus_Monkey_Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Rhesus_Monkey_Prefrontal_Cortex_Average_RPKM"]+0.1) - np.log2(r["Chimp_Prefrontal_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for rhesus monkey
data_table[["p_Hum-Chimp", "p_Hum-Rhesus_Monkey", "p_Rhesus_Monkey_Chimp"]].to_csv("prefrontal_cortex_relative_expression_with_gene_names_RPKM_difference_LOG2.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Primary Visual Cortex RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_relative_RPKM in ["v_Hum-Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Primary_Visual_Cortex_Average_RPKM"]+0.1) - np.log2(r["Chimp_Primary_Visual_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for human
for hum_relative_RPKM in ["v_Hum-Rhesus_Monkey"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Human_Primary_Visual_Cortex_Average_RPKM"]+0.1) - np.log2(r["Rhesus_Monkey_Primary_Visual_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for chimp
for hum_relative_RPKM in ["v_Rhesus_Monkey_Chimp"]:
	col_name = hum_relative_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((np.log2(r["Rhesus_Monkey_Primary_Visual_Cortex_Average_RPKM"]+0.1) - np.log2(r["Chimp_Primary_Visual_Cortex_Average_RPKM"]+0.1))), axis=1)
#calculate relative expression for rhesus monkey
data_table[["v_Hum-Chimp", "v_Hum-Rhesus_Monkey", "v_Rhesus_Monkey_Chimp"]].to_csv("primary_visual_cortex_relative_expression_with_gene_names_RPKM_difference_LOG2.tsv", sep="\t", header=True)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#All Species & Tissue RNA-seq relative expression data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

data_table[["m_Hum-Chimp", "m_Hum-Rhesus_Monkey", "m_Rhesus_Monkey_Chimp", "k_Hum-Chimp", "k_Hum-Rhesus_Monkey", "k_Rhesus_Monkey_Chimp", "c_Hum-Chimp", "c_Hum-Rhesus_Monkey", "c_Rhesus_Monkey_Chimp", "p_Hum-Chimp", "p_Hum-Rhesus_Monkey", "p_Rhesus_Monkey_Chimp", "v_Hum-Chimp", "v_Hum-Rhesus_Monkey", "v_Rhesus_Monkey_Chimp"]].to_csv("all_species_tissue_relative_expression_Average_RPKM_difference_LOG2.tsv", sep="\t", header=True)
