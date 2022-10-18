import pandas as pd
import numpy as np
import itertools as it

#from . import similarities as sim

#from similarities import calculate_similarities_timeseriesobj_to_medoid, calculate_center_medoid

#from similarity.calculate_similarities import calculate_similarities_timeseriesobj_to_medoid, calculate_similarities_clusterobject, calculate_center_medoid

# def _add_obj_to_most_sim_cluster(timeseriesobjs, nonmedoids_ids, cluster_medoids, cluster_keys, simfunc):
#         ##assign objects to the nearest medoids 
#         ##get clustering
#         ##calculate new medoids
#         clusters={}
        
#         #to cython function
#         #get index for the closest medoid for each object

#         clusters, thisSim = sim.fast_calculate_similarities_timeseriesobj_to_medoid(timeseriesobjs, nonmedoids_ids, cluster_medoids, simfunc)
        
#         return clusters, thisSim

# def _new_add_obj_to_most_sim_cluster(nonsampledata, nonsampledata_id, mapped_clusters, prototypes, simfunc):
#     #calculate similarities within clusters
#     #get clusters number for the closest

#     if simfunc == "pearson":
#         mat = [sim.corr_mat(nonsampledata[x], prototypes) for x in range(nonsampledata.shape[0])]
        
#     if simfunc == "euclidean":
#         mat = [sim.dist_mat(nonsampledata[x], prototypes) for x in range(nonsampledata.shape[0])]

#     mat = np.average(mat, axis=0)
#     #maxforobj is the index of medoid which the timeseries is closest
#     maxforobj = np.argmax(mat, axis=1)

#     for id, maxi in enumerate(maxforobj):
#         if str(maxi+1) in mapped_clusters.keys():
#             mapped_clusters[str(maxi+1)].append(nonsampledata_id[id])
#         else:
#             mapped_clusters.setdefault(str(maxi+1), [nonsampledata_id[id]])
    

#     return mapped_clusters

# def within_cluster_similarity(alltimeseries, clusters, similarityfunction):
#     results = np.zeros(len(clusters), dtype="double")
#     c=0
#     for u,v in clusters.items():
#         if similarityfunction == "pearson":
#             mat = [sim.corr_mat(alltimeseries[x,v,:], alltimeseries[x,v,:]) for x in range(alltimeseries.shape[0])]

#         if similarityfunction == "euclidean":
#             mat = [sim.dist_mat(alltimeseries[x,v,:], alltimeseries[x,v,:]) for x in range(alltimeseries.shape[0])]
#         #take average and along axis 1
#         matq = np.average(mat, axis=0)
#         within_sim = np.sum(matq, axis=1)
#         results[c] = within_sim[np.argmax(within_sim)]/len(v)
#         c+=1  

#     return results

def cluster_map_to_item(DataSet, clusters):
    gene_list = DataSet.gene_list
    symbs = DataSet.symbs
    
    g_clusters = {}
    s_clusters = {}


    for u,v in clusters.items():
        g_tmp = [np.array(gene_list)[x] for x in v]
        s_tmp = [np.array(symbs)[x] for x in v]

        g_clusters["Cluster "+str(u)] = g_tmp
        s_clusters["Cluster "+str(u)] = s_tmp

    return g_clusters, s_clusters

def merge_above_similarity_threshold_clusters(clusters, prototypes, threshold=0.9):
    sims = []
    for k in it.combinations(list(prototypes.keys()),2):
        tmp_sims = sim.pearson_corr(np.array(prototypes[k[0]][0]), np.array(prototypes[k[1]][0]), len(prototypes[k[0]]))
        print(k)
        print(tmp_sims)
        # if tmp_sims > threshold:
        #     #make clusters prototypes a dict
        #     clusters[k[0]].append(clusters[k[1]])
        #     clusters.pop(k[1])
        #     prototypes.pop(k[1])

    return clusters, prototypes

    