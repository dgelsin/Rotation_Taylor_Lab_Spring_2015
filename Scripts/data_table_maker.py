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

df = df.ix[:64]

df.to_csv(sys.stdout, sep="\t", header=False)
