import numpy as np
from sklearn.cluster import *
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import pairwise_distances, silhouette_score
from tslearn.metrics import cdist_soft_dtw
from tslearn.clustering import TimeSeriesKMeans
from tslearn.clustering import silhouette_score as tssilhouette_score
from sklearn.metrics import davies_bouldin_score
from sklearn.metrics import pairwise_distances, silhouette_score

PROTOTYPE = {'median':np.median, 'mean':np.mean}

def intra_clusters_similarity(timeserieslist, clusters, prototype, metric):
    ##get array of gene in order of clustering
    a = [y for x in clusters.index_clusters.values() for y in x]
    ##get array of cluster in order of clustering
    b = [u for u,v in clusters.index_clusters.items() for _ in v]
    ## arg sort the gene array
    ## map the cluster array with arg sorted gene array
    cluster_array = [b[i] for i in np.argsort(a)]
    ## calculate similarity the time series object with corresponding prototype
    prototype_arr = [clusters._prototype[x] for x in cluster_array]

    assert len(prototype_arr) == timeserieslist.shape[1]
    if metric == "soft_dtw":
        proto_cluster_sim_arr = [cdist_soft_dtw([prototype_arr[x], np.squeeze(timeserieslist.tiu_timeserieslist[:,x,:],axis=0)]) for x in range(len(prototype_arr))]
    else:
        proto_cluster_sim_arr = [pairwise_distances([prototype_arr[x], np.squeeze(timeserieslist.tiu_timeserieslist[:,x,:],axis=0)], metric="euclidean")[0,1] for x in range(len(prototype_arr))]

    ## calculate average similarities of the clusters 
    all_clu_sim = []
    for clu in clusters.index_clusters.keys():
        this_clu = [x==clu for x in cluster_array]
        avg_sim = 0
        for sim_clu in range(len(this_clu)):
            if this_clu[sim_clu]:
                avg_sim+=proto_cluster_sim_arr[sim_clu]
        avg_sim = avg_sim/np.sum(this_clu)
        all_clu_sim.append(avg_sim)
    
    return all_clu_sim