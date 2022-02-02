import pandas as pd
import numpy as np
import time
import random
import warnings
from numba import jit, njit, prange
import multiprocessing as mp

from inputdata.DataSet import DataSet
from inputdata.BioNetwork import BioNetwork
from similarity.calculate_similarities import _calculate_similarities_timeseriesobj_to_medoid



def _add_obj_to_most_sim_cluster(nonsampledata, nonsampledata_id, mapped_clusters, prototypes, simfunc):
    #calculate similarities within clusters
    if not hasattr(nonsampledata, "shape"):
        nonsampledata = np.array(nonsampledata)

    clusters = dict([(str(d), []) for d in range(len(prototypes))])
    #get clusters number for the closest
    minid = _calculate_similarities_timeseriesobj_to_medoid(nonsampledata, prototypes)


    for c,d in enumerate(minid):
        clusters[str(d)].append(nonsampledata_id[c])

    return clusters

def cluster_map_to_item(DataSet, clusters):
    gene_list = DataSet.gene_list[0]
    symbs = DataSet.symbs[0]
    
    g_clusters = {}
    s_clusters = {}

    for u,v in clusters.items():
        g_tmp = [gene_list[x] for x in v]
        s_tmp = [symbs[x] for x in v]

        g_clusters["Cluster "+str(u)] = g_tmp
        s_clusters = s_tmp

    return g_clusters, s_clusters