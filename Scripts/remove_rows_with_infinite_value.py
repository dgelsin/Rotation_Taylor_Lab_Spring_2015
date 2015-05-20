#!/usr/bin/env python

import sys
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 
from pylab import show

f=sys.argv[1]

df=DataFrame.from_csv(f, header=0, sep="\t")

data_table= df.dropna(axis=0, how='any')

#df = df.drop(["inf"], axis=0, inplace=True)

data_table.to_csv("muscle_relative_expression_with_gene_names_noNaN_RPKM.tsv", sep="\t", header=True)


f=sys.argv[2]

df=DataFrame.from_csv(f, header=0, sep="\t")

data_table= df.dropna(axis=0, how='any')

#df = df.drop(["inf"], axis=0, inplace=True)

data_table.to_csv("kidney_relative_expression_with_gene_names_noNaN_RPKM.tsv", sep="\t", header=True)


f=sys.argv[3]

df=DataFrame.from_csv(f, header=0, sep="\t")

data_table= df.dropna(axis=0, how='any')

#df = df.drop(["inf"], axis=0, inplace=True)

data_table.to_csv("cerebellar_cortex_relative_expression_with_gene_names_noNaN_RPKM.tsv", sep="\t", header=True)


f=sys.argv[4]

df=DataFrame.from_csv(f, header=0, sep="\t")

data_table= df.dropna(axis=0, how='any')

#df = df.drop(["inf"], axis=0, inplace=True)

data_table.to_csv("prefrontal_cortex_relative_expression_with_gene_names_noNaN_RPKM.tsv", sep="\t", header=True)


f=sys.argv[5]

df=DataFrame.from_csv(f, header=0, sep="\t")

data_table= df.dropna(axis=0, how='any')

#df = df.drop(["inf"], axis=0, inplace=True)

data_table.to_csv("primary_visual_cortex_relative_expression_with_gene_names_noNaN_RPKM.tsv", sep="\t", header=True)


f=sys.argv[6]

df=DataFrame.from_csv(f, header=0, sep="\t")

data_table= df.dropna(axis=0, how='any')

#df = df.drop(["inf"], axis=0, inplace=True)

data_table.to_csv("all_species_tissue_relative_expression_with_gene_names_noNaN_RPKM.tsv", sep="\t", header=True)

