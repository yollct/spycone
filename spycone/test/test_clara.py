import matplotlib.pyplot as plt
from matplotlib import cm, colors
import pandas as pd
import numpy as np
import itertools as it
import time
import cProfile 
import networkx as nx
import cProfile, pstats, io
pr = cProfile.Profile()
import pyximport
pyximport.install()
from sklearn.metrics import silhouette_samples, silhouette_score
import seaborn as sns
import sys

from ticone.inputdata.DataSet import DataSet
from ticone.inputdata.BioNetwork import BioNetwork
from ticone.similarity.similarity_matrix import similarity_matrix
from ticone.clustering.hcluster import hierarchical_cluster
from ticone.clustering.pamk_cluster import PAMK_cluster
from ticone.clustering.clara import clara_cluster
#from ticone.visualization.plotting import plotting_clusters, extract_subnetworks
from ticone.visualization.plotting_clusters import plotting_clusters
from ticone.preprocessing._preprocessing import preprocess
#from stat.pq_val import _cal_pqvals
#from enrichment.goterms import go_enrich
from ticone.util_stat.qc import plot_silhouette, _qc_silhouette
from ticone.shuffle import shuffle
from ticone.connectivity.connectivity_results import connectivity
from ticone.feature.get_features import featuresObj
from ticone.clustering import similarities as sim

sys.stdout = open('output.txt','wt')
data1 =  pd.read_csv("/home/chit/Desktop/lrz_ticone/play_sample/src_test_resources_test_files_clustering_qual-patterns3-objects80-seed42-bg120-samples2-conftrue-tp5_mergeTestCase.txt", sep="\t", header=None)


data1.columns = ["id","rep", "obj","tp1","tp2","tp3","tp4","tp5"]
del data1['id']
sample0 = data1[data1.iloc[:,0]=="sample0"]
sample1 = data1[data1.iloc[:,0]=="sample1"]
mergeddata1 = pd.merge(sample0.iloc[:,1:],sample1.iloc[:,1:], on="obj")
gene_list = mergeddata1['obj']
ts = mergeddata1.iloc[:,1:]

gendata1 = DataSet(ts=ts, 
        gene_list = gene_list, 
        reps1 = 2, 
        timepts = 5,
        discretization_steps=[4,4])
# print(sample0.shape)

# clara = PAMK_cluster(gendata1, n_cluster=10, similarityfunction = "pearson")
# clara.find_clusters()
clara = clara_cluster(gendata1, n_cluster=10, n_start=5, n_samples=3, similarityfunction="pearson", prototypefunction="mean", seed=14)
clara.find_clusters()
print(np.array(clara.genelist_clusters))
print(np.array(clara._prototype))



