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
#print "Step 1: Making dataframe" sys.stdout.flush()
df_ensID=DataFrame.from_csv(f, header=0, sep="\t") #make a data frame that will have its index values replaced with values from the dictionary that will be made below (1st file input)
#print df_ensID

#df_gene_name=DataFrame.from_csv(g, header=0, sep="\t") #make a data frame that will be used to make a dictionary (2nd file input)
#print df_gene_name
#print "Done with Step 1..." sys.stdout.flush()

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Make a dictionary of gene names to ensemble IDs (key = gene name)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#how to make a dictionary with pandas dataframe
#dictionary2 = df_gene_name.set_index(df_gene_name.index).to_dict()
#print len(dictionary2)

##########
#how to make a dictionary with opening the file directly
##########
#print "Step 2: Making dictionary" sys.stdout.flush()
with open(g, 'rU') as document:
	gene_name = {} # initiate dictionary
	for line in document:
		key, value = line[:-1].split("\t") #split by tab delimiter
		gene_name[key] = value
		#[:-1] give everything but the last character (newline character in this case)
#print gene_name
#print "Done with Step 2..."

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Replace ENSEMBLE IDS with gene name and rename IDs in index that don't match with a new name
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#index = df_ensID.index
#print index
#print "Step 3: Replacing ENSEMBL_IDs with gene names"
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
#print "Done with Step 3..."
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Remove blank (NaN) index names
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#print "Step 4: Removing lines with blank (NaN) gene names"
df_ensID= df_ensID.dropna(axis=0)
#dropna removes NaN lines, axis=0 specifies to look only in the index (ie only droplines where it there is no gene name)
#print "Done with Step 4..."
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Remove duplicate gene names
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#print "Step 5: Removing duplicate gene names"
df_ensID["index"] = df_ensID.index
#cannot specify to remove duplicate index names with drop_duplicates() so have to make the index into a column called "index"
df_ensID.drop_duplicates(cols='index', take_last=True, inplace=True)
#now drop any duplicate genes in the new index column
del df_ensID["index"]
#remove that index column
#print "Done with Step 5..."
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Write both data frames to file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#print "Writing to file"


#drop all index rows with the name "iso" and return all index rows that do not have that
df_ensID.to_csv(sys.stdout, sep="\t", header=False)
#print "All finished :)"



