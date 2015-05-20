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



df["index"] = df.index
df.drop_duplicates(cols='index', take_last=True, inplace=True)
del df["index"]

#df = df.drop(["inf"], axis=0, inplace=True)

df.to_csv("readcounts_matched_no_duplicates", sep="\t", header=True)