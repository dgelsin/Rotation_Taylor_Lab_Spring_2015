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

data_table= df.dropna(axis=0)

#df = df.drop(["inf"], axis=0, inplace=True)

data_table.to_csv("output_noNAN", sep="\t", header=True)