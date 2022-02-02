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
from array import array

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

alltimeseries = np.array([[[1,2,3],[3,4,5],[4,5,6]], [[1,2,2],[3,4,4],[4,5,5]], [[1,3,3],[3,5,5],[4,6,6]]], dtype="double")
obj_comp = [np.array([1,2,3], dtype="double"),np.array([3,4,5], dtype="double"),np.array([4,5,6], dtype="double")]
simfunc = "euclidean"

result = sim.calculate_center_medoid(alltimeseries, {"1":[0,1,2]}, "pearson")
print(np.array(result[0]))

###before that set cpdef in similarities.pyx
print("Test pearson corr function")
result = sim.pearson_corr(np.array([0,1,2,3,4],dtype="double"), np.array([0,-1,-2,-3,-4],dtype="double"), 5)
assert result == (-1+1)/2

result = sim.pearson_corr(np.array([2,43,3.5],dtype="double"), np.array([10,3,9],dtype="double"), 2)
assert result == round(0.0025074884692379285,1)

print("Test within cluster similarity")
alltimeseries = np.array([[[0,1,2,3,4],[0,1,2,2,4]],[[0,1,2,3,3],[0,3,4,3,5]]], dtype="double")
result = sim.within_cluster_similarity(alltimeseries, {'1':[0,1]}, "pearson")
print(np.array(result))
assert result[0] ==0.9379408266263856


print("Test euclidean")
result= sim.euclidean_dist(np.array([0,1,2,3,4],dtype="double"), np.array([0,1,2,3,4], dtype="double"), 5)
assert result == 0

#generate a dataset
dataset_ts = np.zeros(1,100,5)
for l in dataset_ts.shape[0]:
    for k in dataset_ts.shape[1]:
        dataset_ts[l][k] = [np.random.normal(10,2) for _ in range(dataset_ts.shape[2])]

objlist = ["obj{}".format(x) for x in range(dataset_ts.shape[1])]
