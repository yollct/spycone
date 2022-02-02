import pandas as pd
import numpy as np
import time
import random
import warnings
import math
import itertools as it
from sklearn import metrics

from ..DataSet import DataSet
from ..BioNetwork import BioNetwork
from .._clustering.clusterobj import clusterObj
from .._clustering import similarities as sim

def silhouette_index(clusterObj):
    data = np.array(clusterObj.DataSet.timeserieslist, dtype="double")

    genelist = clusterObj.DataSet.gene_list
    replicates = data.shape[0]

    clustering = clusterObj.index_clusters
    
    pred = clusterObj._labels
    print(pred)

    sil = 0
    for r in range(replicates):
        sil += metrics.silhouette_score(data[r,:,:], pred, metric=clusterObj.metric) 
        #sil += silhouette(data[r,:,:], clustering, simfunc=simfunc)
    
    return sil/replicates

