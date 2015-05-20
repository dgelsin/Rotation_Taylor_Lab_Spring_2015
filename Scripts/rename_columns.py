#!/usr/bin/env python

import sys
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 
from pylab import show

f = sys.argv[1]
data_table=DataFrame.from_csv(f, header=0, sep="\t")
data_table.columns = ["hM1","hM2","hM3","hM4","hM5","hM6","hK1","hK2","hK3","hK4","hK5","hK6","hP1","hP2","hP3","hP4",	"hP5",	"hP6",	"hV1",	"hV2",	"hV3",	"hV4",	"hV5",	"hV6",	"hC1",	"hC2",	"hC3",	"hC4",	"hC5",	"hC6",	"cM1",	"cM2",	"cM3",	"cM4",	"cM5",	"cM6",	"cK1",	"cK2",	"cK3",	"cK4",	"cK5",	"cK6",	"cP1",	"cP2",	"cP3",	"cP4",	"cP5",	"cP6",	"cV1",	"cV2",	"cV3",	"cV4",	"cV5",	"cV6",	"cC1",	"cC2",	"cC3",	"cC4",	"cC5",	"cC6",	"rM1",	"rM2",	"rM3",	"rM4",	"rM5",	"rM6",	"rK1",	"rK2",	"rK3",	"rK4",	"rK5",	"rK6",	"rP1",	"rP2",	"rP3",	"rP4",	"rP5",	"rP6",	"rV1",	"rV2",	"rV3",	"rV4",	"rV5",	"rV6",	"rC1",	"rC2",	"rC3",	"rC4",	"rC5",	"rC6",	"mM1",	"mM2",	"mM3",	"mM4",	"mM5",	"mM6",	"mK1",	"mK2",	"mK3",	"mK4",	"mK5",	"mK6",	"mP1",	"mP2",	"mP3",	"mP4",	"mP5",	"mP6",	"mV1",	"mV2",	"mV3",	"mV4",	"mV5",	"mV6",	"mC1",	"mC2",	"mC3", "mC4", "mC5", "mC6"]

print data_table
data_table.to_csv("RNA-seq_hum_chimp-mod.tsv", sep = "\t")

