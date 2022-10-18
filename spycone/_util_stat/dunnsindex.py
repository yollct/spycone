import numpy as np

from .._clustering import similarities as sim

def dunnsindex(clusterObj):
    clusters = clusterObj.index_clusters
    alltimeseries = clusterObj.DataSet.timeserieslist
    simfunc = clusterObj.metric


    #get the min clusters pair similarity (min dist/max corr)
    minmin = sim.clusters_pair_similarity(clusters, simfunc)
    inter_min = np.max(minmin)
    print(np.array(minmin))

    #max intra cluster sim
    #(max dist/min cor)
    intra_max = np.min(sim.intra_clusters_similarity(alltimeseries, clusters, clusterObj._prototype, simfunc))
    print(np.array(intra_max))
    #min/max    
    return inter_min/intra_max