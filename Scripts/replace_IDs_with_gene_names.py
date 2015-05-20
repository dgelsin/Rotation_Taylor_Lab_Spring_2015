#!/usr/bin/env python

import sys
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 
from pylab import show

f=sys.argv[1]
g=sys.argv[2]
#i=sys.argv[3]

df_ensID=DataFrame.from_csv(f, header=0, sep="\t") #make a data frame that will have its index values replaced with values from the dictionary that will be made below (1st file input)
#print df_ensID

#df_gene_name=DataFrame.from_csv(g, header=0, sep="\t") #make a data frame that will be used to make a dictionary (2nd file input)
#print df_gene_name


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Make a dictionary of gene names to ensemble IDs (key = gene name)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#how to make a dictionary with pandas dataframe
#dictionary2 = df_gene_name.set_index(df_gene_name.index).to_dict()
#print len(dictionary2)

##########
#how to make a dictionary with opening the file directly
##########
with open(g, 'rU') as document:
	gene_name = {} # initiate dictionary
	for line in document:
		key, value = line[:-1].split("\t") #split by tab delimiter
		gene_name[key] = value
		#[:-1] give everything but the last character (newline character in this case)
#print gene_name


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Replace ENSEMBLE IDS with gene name and rename IDs in index that don't match with a new name
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#index = df_ensID.index
#print index

index = np.array(df_ensID.index).ravel().tolist() # make index values into a list
isoform = []
for i in index:
	if i in gene_name.keys():
		genes = gene_name[i]
		isoform.append(genes)
	if i not in gene_name.keys():
		isoform.append("iso")
#make a list of the gene names with a '0' inserted if there is no match into the list isoform
#print isoform
#print len(isoform)

df_ensID.index = isoform
#df_ensID.reset_index(level=0, inplace=True)
#this makes the index into a column, and the new index to be numberically numbered

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Take out ENSEMBL IDS (now named "iso") with no match
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#df_iso = df_ensID.ix["iso"]
#append all index rows that did not have an ENSEMBL ID match into a new data frame

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Write both data frames to file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#df_iso.to_csv("isoform.tsv", sep="\t", header=True)

#df_ensID.drop(["iso"], axis=0, inplace=True)
#drop all index rows with the name "iso" and return all index rows that do not have that
df_ensID.to_csv("genes.tsv", sep="\t", header=True)




