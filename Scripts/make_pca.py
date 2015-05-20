#!/usr/bin/env python
import sys
import numpy
import pandas
import matplotlib.pyplot as plot
from sklearn.decomposition import PCA
from pandas import DataFrame
from matplotlib.mlab import PCA as mlabPCA



f = sys.argv[1]
df=DataFrame.from_csv(f, header=0, sep="\t")

df.columns = ["","","","","","","","","","","","","","",""]
#df.columns=["Human Muscle", "Chimp Muscle", "Rhesus Monkey Muscle",	"Mouse Muscle",	"Human Kidney", "Chimp Kidney", "Rhesus Monkey Kidney",	"Mouse Kidney", "Human CC", "Chimp CC", "Rhesus Monkey CC", "Mouse CC",	"Human PFC", "Chimp PFC", "Rhesus Monkey PFC", "Mouse PFC", "Human VC", "Chimp VC", "Rhesus Monkey VC", "Mouse VC"]
#df.columns=["hM1", "hM2", "hM3", "hM4", "hM5", "hM6", "hK1",	"hK2",	"hK3",	"hK4",	"hK5",	"hK6",	"hP1",	"hP2",	"hP3",	"hP4",	"hP5",	"hP6",	"hV1",	"hV2",	"hV3",	"hV4",	"hV5",	"hV6",	"hC1",	"hC2",	"hC3",	"hC4",	"hC5",	"hC6",	"cM1",	"cM2",	"cM3",	"cM4",	"cM5",	"cM6",	"cK1",	"cK2",	"cK3",	"cK4",	"cK5",	"cK6",	"cP1",	"cP2",	"cP3",	"cP4",	"cP5",	"cP6",	"cV1",	"cV2",	"cV3",	"cV4",	"cV5",	"cV6",	"cC1",	"cC2",	"cC3",	"cC4",	"cC5",	"cC6",	"rM1",	"rM2",	"rM3",	"rM4",	"rM5",	"rM6",	"rK1",	"rK2",	"rK3",	"rK4",	"rK5",	"rK6",	"rP1",	"rP2",	"rP3",	"rP4",	"rP5",	"rP6",	"rV1",	"rV2",	"rV3",	"rV4",	"rV5",	"rV6",	"rC1",	"rC2",	"rC3",	"rC4",	"rC5",	"rC6"]
#df.columns=["hM","cM",	"rM", "hK",	"cK",	"rK",	"hC",	"cC",	"rC",	"hP",	"cP", "rP",	"hV",	"cV","rV"]
#df.columns=['m_hum/chimp', 'm_hum/rhesus','m_rhesus/chimp', 'k_hum/chimp', 'k_hum/rhesus','k_rhesus/chimp', 'c_hum/chimp', 'c_hum/rhesus','c_rhesus/chimp', 'p_hum/chimp', 'p_hum/rhesus','p_rhesus/chimp', 'v_hum/chimp', 'v_hum/rhesus','v_rhesus/chimp']


col_names = df.columns.values.tolist()
#print col_names
df = df.T
#transpose the dataframe

n_samples =15
n_features = 3
pca = PCA(n_components=15, whiten=False)

fit = pca.fit(df)

plot.bar(range(15), fit.explained_variance_ratio_ )

x = fit.transform(df)

plot.figure()
fig = plot.scatter(x[:,0], x[:,1], c=(u'g', u'r', u'b'), marker=('o'), s=100)

for i in range(15):
    plot.gca().annotate(col_names[i], (x[i,0], x[i,1]))

#plot.title("Principle Component Plot: Gene Expression between Tissue and Species Type", fontsize = 26)


plot.tick_params(axis='x', which='both', bottom='off', labelbottom='off')
plot.tick_params(axis='y', which='both', left='off', labelleft='off')
#plot.axis('off')
plot.xlabel('component 1')
plot.ylabel('component 2')

#plot.legend(('Human' 'Chimp' 'Rhesus Monkey'), ncol=3, scatterpoints=3, loc='lower left', fontsize=16)

plot.savefig("no_label_average-difference_log2_RPKM_pca1.png",bbox_inches='tight',dpi=300)