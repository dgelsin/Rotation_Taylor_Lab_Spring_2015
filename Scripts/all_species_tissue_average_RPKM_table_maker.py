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
#Muscle tissue RNA-seq average data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_mus_average_RPKM in ["Human_Muscle_Average_RPKM"]:
	col_name = hum_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["hM1"] + r["hM2"] + r["hM3"] + r["hM4"] + r["hM5"] + r["hM6"])/6), axis=1) #average each row
#average each row of the human columns and write the values in a new column named "Human_Muscle_Average_RPKM"
for chimp_mus_average_RPKM in ["Chimp_Muscle_Average_RPKM"]:
	col_name = chimp_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["cM1"] + r["cM2"] + r["cM3"] + r["cM4"] + r["cM5"] + r["cM6"])/6), axis=1)
#average each row of the chimp columns and write the values in a new column named "Chimp_Muscle_Average_RPKM"
for Rh_Mon_mus_average_RPKM in ["Rhesus_Monkey_Muscle_Average_RPKM"]:
	col_name = Rh_Mon_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["rM1"] + r["rM2"] + r["rM3"] + r["rM4"] + r["rM5"] + r["rM6"])/6), axis=1)
#average each row of the rhesus monkey columns and write the values in a new column named "Rhesus_Monkey_Muscle_Average_RPKM"
for Mouse_mus_average_RPKM in ["Mouse_Muscle_Average_RPKM"]:
	col_name = Mouse_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["mM1"] + r["mM2"] + r["mM3"] + r["mM4"] + r["mM5"] + r["mM6"])/6), axis=1)
#average each row of the mouse columns and write the values in a new column named "Mouse_Muscle_Average_RPKM"

data_table[["Human_Muscle_Average_RPKM", "Chimp_Muscle_Average_RPKM", "Rhesus_Monkey_Muscle_Average_RPKM", "Mouse_Muscle_Average_RPKM"]].to_csv("muscle_Average_RPKM_with_gene_names.tsv", sep="\t", header=True)
#write the averaged columns to a tsv file

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Kidney tissue average RNA seq data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

for hum_mus_average_RPKM in ["Human_Kidney_Average_RPKM"]:
	col_name = hum_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["hK1"] + r["hK2"] + r["hK3"] + r["hK4"] + r["hK5"] + r["hK6"])/6), axis=1) #average each row
#average each row of the human columns and write the values in a new column named "Human_Muscle_Average_RPKM"
for chimp_mus_average_RPKM in ["Chimp_Kidney_Average_RPKM"]:
	col_name = chimp_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["cK1"] + r["cK2"] + r["cK3"] + r["cK4"] + r["cK5"] + r["cK6"])/6), axis=1)
#average each row of the chimp columns and write the values in a new column named "Chimp_Muscle_Average_RPKM"
for Rh_Mon_mus_average_RPKM in ["Rhesus_Monkey_Kidney_Average_RPKM"]:
	col_name = Rh_Mon_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["rK1"] + r["rK2"] + r["rK3"] + r["rK4"] + r["rK5"] + r["rK6"])/6), axis=1)
#average each row of the rhesus monkey columns and write the values in a new column named "Rhesus_Monkey_Muscle_Average_RPKM"
for Mouse_mus_average_RPKM in ["Mouse_Kidney_Average_RPKM"]:
	col_name = Mouse_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["mK1"] + r["mK2"] + r["mK3"] + r["mK4"] + r["mK5"] + r["mK6"])/6), axis=1)
#average each row of the mouse columns and write the values in a new column named "Mouse_Muscle_Average_RPKM"
data_table[["Human_Kidney_Average_RPKM", "Chimp_Kidney_Average_RPKM", "Rhesus_Monkey_Kidney_Average_RPKM", "Mouse_Kidney_Average_RPKM"]].to_csv("kidney_Average_RPKM_with_gene_names.tsv", sep="\t", header=True)
#write the averaged columns to a tsv file

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Cerebellar cortex tissue average RNA seq data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for hum_mus_average_RPKM in ["Human_Cerebellar_Cortex_Average_RPKM"]:
	col_name = hum_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["hC1"] + r["hC2"] + r["hC3"] + r["hC4"] + r["hC5"] + r["hC6"])/6), axis=1) #average each row
#average each row of the human columns and write the values in a new column named "Human_Muscle_Average_RPKM"
for chimp_mus_average_RPKM in ["Chimp_Cerebellar_Cortex_Average_RPKM"]:
	col_name = chimp_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["cC1"] + r["cC2"] + r["cC3"] + r["cC4"] + r["cC5"] + r["cC6"])/6), axis=1)
#average each row of the chimp columns and write the values in a new column named "Chimp_Muscle_Average_RPKM"
for Rh_Mon_mus_average_RPKM in ["Rhesus_Monkey_Cerebellar_Cortex_Average_RPKM"]:
	col_name = Rh_Mon_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["rC1"] + r["rC2"] + r["rC3"] + r["rC4"] + r["rC5"] + r["rC6"])/6), axis=1)
#average each row of the rhesus monkey columns and write the values in a new column named "Rhesus_Monkey_Muscle_Average_RPKM"
for Mouse_mus_average_RPKM in ["Mouse_Cerebellar_Cortex_Average_RPKM"]:
	col_name = Mouse_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["mC1"] + r["mC2"] + r["mC3"] + r["mC4"] + r["mC5"] + r["mC6"])/6), axis=1)
#average each row of the mouse columns and write the values in a new column named "Mouse_Muscle_Average_RPKM"
data_table[["Human_Cerebellar_Cortex_Average_RPKM", "Chimp_Cerebellar_Cortex_Average_RPKM", "Rhesus_Monkey_Cerebellar_Cortex_Average_RPKM", "Mouse_Cerebellar_Cortex_Average_RPKM"]].to_csv("Cerebellar_Cortex_Average_RPKM_with_gene_names.tsv", sep="\t", header=True)
#write the averaged columns to a tsv file

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Prefrontal cortex tissue average RNA seq data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

for hum_mus_average_RPKM in ["Human_Prefrontal_Cortex_Average_RPKM"]:
	col_name = hum_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["hP1"] + r["hP2"] + r["hP3"] + r["hP4"] + r["hP5"] + r["hP6"])/6), axis=1) #average each row
#average each row of the human columns and write the values in a new column named "Human_Muscle_Average_RPKM"
for chimp_mus_average_RPKM in ["Chimp_Prefrontal_Cortex_Average_RPKM"]:
	col_name = chimp_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["cP1"] + r["cP2"] + r["cP3"] + r["cP4"] + r["cP5"] + r["cP6"])/6), axis=1)
#average each row of the chimp columns and write the values in a new column named "Chimp_Muscle_Average_RPKM"
for Rh_Mon_mus_average_RPKM in ["Rhesus_Monkey_Prefrontal_Cortex_Average_RPKM"]:
	col_name = Rh_Mon_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["rP1"] + r["rP2"] + r["rP3"] + r["rP4"] + r["rP5"] + r["rP6"])/6), axis=1)
#average each row of the rhesus monkey columns and write the values in a new column named "Rhesus_Monkey_Muscle_Average_RPKM"
for Mouse_mus_average_RPKM in ["Mouse_Prefrontal_Cortex_Average_RPKM"]:
	col_name = Mouse_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["mP1"] + r["mP2"] + r["mP3"] + r["mP4"] + r["mP5"] + r["mP6"])/6), axis=1)
#average each row of the mouse columns and write the values in a new column named "Mouse_Muscle_Average_RPKM"
data_table[["Human_Prefrontal_Cortex_Average_RPKM", "Chimp_Prefrontal_Cortex_Average_RPKM", "Rhesus_Monkey_Prefrontal_Cortex_Average_RPKM", "Mouse_Prefrontal_Cortex_Average_RPKM"]].to_csv("Prefrontal_Cortex_Average_RPKM_with_gene_names.tsv", sep="\t", header=True)
#write the averaged columns to a tsv file

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Primary visual cortex tissue average RNA seq data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

for hum_mus_average_RPKM in ["Human_Primary_Visual_Cortex_Average_RPKM"]:
	col_name = hum_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["hV1"] + r["hV2"] + r["hV3"] + r["hV4"] + r["hV5"] + r["hV6"])/6), axis=1) #average each row
#average each row of the human columns and write the values in a new column named "Human_Muscle_Average_RPKM"
for chimp_mus_average_RPKM in ["Chimp_Primary_Visual_Cortex_Average_RPKM"]:
	col_name = chimp_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["cV1"] + r["cV2"] + r["cV3"] + r["cV4"] + r["cV5"] + r["cV6"])/6), axis=1)
#average each row of the chimp columns and write the values in a new column named "Chimp_Muscle_Average_RPKM"
for Rh_Mon_mus_average_RPKM in ["Rhesus_Monkey_Primary_Visual_Cortex_Average_RPKM"]:
	col_name = Rh_Mon_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["rV1"] + r["rV2"] + r["rV3"] + r["rV4"] + r["rV5"] + r["rV6"])/6), axis=1)
#average each row of the rhesus monkey columns and write the values in a new column named "Rhesus_Monkey_Muscle_Average_RPKM"
for Mouse_mus_average_RPKM in ["Mouse_Primary_Visual_Cortex_Average_RPKM"]:
	col_name = Mouse_mus_average_RPKM
	data_table[col_name] = data_table.apply(lambda r: ((r["mV1"] + r["mV2"] + r["mV3"] + r["mV4"] + r["mV5"] + r["mV6"])/6), axis=1)
#average each row of the mouse columns and write the values in a new column named "Mouse_Muscle_Average_RPKM"
data_table[["Human_Primary_Visual_Cortex_Average_RPKM", "Chimp_Primary_Visual_Cortex_Average_RPKM", "Rhesus_Monkey_Primary_Visual_Cortex_Average_RPKM", "Mouse_Primary_Visual_Cortex_Average_RPKM"]].to_csv("Primary_Visual_Cortex_Average_RPKM_with_gene_names.tsv", sep="\t", header=True)
#write the averaged columns to a tsv file

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#All Species & Tissue RNA-seq Avg RPKM data
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

data_table[["Human_Muscle_Average_RPKM", "Chimp_Muscle_Average_RPKM", "Rhesus_Monkey_Muscle_Average_RPKM", "Mouse_Muscle_Average_RPKM", "Human_Kidney_Average_RPKM", "Chimp_Kidney_Average_RPKM", "Rhesus_Monkey_Kidney_Average_RPKM", "Mouse_Kidney_Average_RPKM","Human_Cerebellar_Cortex_Average_RPKM", "Chimp_Cerebellar_Cortex_Average_RPKM", "Rhesus_Monkey_Cerebellar_Cortex_Average_RPKM", "Mouse_Cerebellar_Cortex_Average_RPKM","Human_Prefrontal_Cortex_Average_RPKM", "Chimp_Prefrontal_Cortex_Average_RPKM", "Rhesus_Monkey_Prefrontal_Cortex_Average_RPKM", "Mouse_Prefrontal_Cortex_Average_RPKM","Human_Primary_Visual_Cortex_Average_RPKM", "Chimp_Primary_Visual_Cortex_Average_RPKM", "Rhesus_Monkey_Primary_Visual_Cortex_Average_RPKM", "Mouse_Primary_Visual_Cortex_Average_RPKM"]].to_csv("all_species_tissue_Average_RPKM_with_gene_names.tsv", sep="\t", header=True)
#write the averaged columns to a tsv file

