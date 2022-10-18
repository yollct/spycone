import numpy as np

from .._clustering import similarities as sim


def dbindex(clusterObj, simfunc="pearson"):
    alltimeseries = np.array(clusterObj.DataSet.timeserieslist, dtype="double")
    clusters = clusterObj.index_clusters
    prototype= clusterObj._prototype
    simfunc = clusterObj.similarityfunction
    #average dist/corr between each cluster and its medoid
    avg_per_clusters = sim.intra_clusters_similarity(alltimeseries, clusters, prototype, simfunc)

    #compute distance of each clusters pairs
    clusters_pair_sim = sim.clusters_pair_similarity(clusterObj._prototype, simfunc)

    #sum of average dist of each clusters pairs / their distances --> get the max
    values=np.zeros(len(clusters))
    clist = list(clusters.keys())

    dbi=0
    n=0
    for o in range(len(clist)):
        k=clist[o]
        m=0
        for l in clist[0:o]+clist[o+1:]:
            values[o] = (avg_per_clusters[o]+avg_per_clusters[m])/clusters_pair_sim[n][m]

            
            m+=1
        dbi += np.max(values)
        n+=1

    result = dbi/len(clusters)
    
    
    return result