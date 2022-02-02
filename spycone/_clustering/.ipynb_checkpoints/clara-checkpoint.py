import pandas as pd
import numpy as np
import time
import random
import warnings
from numba import jit, njit, prange
import multiprocessing as mp

from inputdata.DataSet import DataSet
from inputdata.BioNetwork import BioNetwork
from similarity.calculate_similarities import _calculate_similarities_timeseriesobj_to_medoid, _calculate_similarities_clusterobject
from prototype.prototype_builder import _prototype_builder
from util_stat.qc import _qc_silhouette
from clustering.pamk_cluster import PAMK_cluster
from clustering.build_cluster import _add_obj_to_most_sim_cluster, cluster_map_to_item

class clara_cluster:
    def __init__(self,
                DataSet,
                n_cluster,
                n_start,
                similarity = "euclidean",
                prototype = "median", 
                n_samples = 5,
                samplesize=None,
                cluster_prototype=None,
                ):
        self.DataSet = DataSet 
        self.n_samples = n_samples
        self.samplesize = samplesize
        self.prototype = prototype
        self.similarity = similarity
        self.n_cluster = n_cluster
        self.n_start = n_start
        self.cluster_prototype = []
        

        def recommended_samplesize(self):
            n = self.DataSet.ts[0].shape[0]
            k = self.n_cluster

            if n <= 250:
                return n
            
            return 100+2*k
        
        if samplesize == None:
            self.samplesize = recommended_samplesize(self)

    def find_sample_data(self, alltimeseriesobject, samplesize):
        sampledata = np.array(alltimeseriesobject.copy())
        indexlist=[]
        
        if len(sampledata) > samplesize:
            chosen_index = np.random.randint(0,len(sampledata), size=(samplesize,))
            sampledata = [sampledata[x] for x in chosen_index]
            
        return sampledata, chosen_index


    def mapping_clusters(self, id_list, clusters):
        mapped_clusters = {}
        for u,v in clusters.items():
            tmp = 0
            tmp = [id_list[x] for x in v]
            
            mapped_clusters[str(u)] = tmp
        
        return mapped_clusters


    def find_clusters(self):
        print("CLARA clustering initiated.")
        clusters = {}
        alltimeseriesobjects = self.DataSet.ts[0]
        n_cluster = self.n_cluster
        n_start = self.n_start
        samplesize = self.samplesize
        n_samples = self.n_samples
        simfunc = self.similarity

        bestSim = -1
        bestclusters = {}
        for i in range(n_samples):
            sampledata, chosen_index = self.find_sample_data(alltimeseriesobjects, samplesize)
            nonsampledata  = np.delete(alltimeseriesobjects, chosen_index, axis=0)
            nonsampledata_id = np.delete([n for n in range(alltimeseriesobjects.shape[0]+1)], chosen_index)
            
            #run pamk
            pam = PAMK_cluster(sampledata, n_cluster=self.n_cluster, prototype=self.prototype, for_clara=True)
            clusters = pam.find_clusters()

            #map index from all objects to subset
            mapped_clusters = self.mapping_clusters(chosen_index, clusters) 

            #add objects to the nearest cluster
            clusters = _add_obj_to_most_sim_cluster(nonsampledata, nonsampledata_id, mapped_clusters,  pam.cluster_prototype[0], simfunc)
            
            
            #get similarity scores
            nouse, s_score = _qc_silhouette(alltimeseriesobjects, clusters, l=1)

            

            if bestSim >= s_score:
                bestclusters = clusters
                cluster_prototypes = _prototype_builder(alltimeseriesobjects, n_cluster, clusters, simfunc)
                try:
                    self.cluster_prototype[0]
                    self.cluster_prototype[0] = cluster_prototypes
                except:
                    self.cluster_prototype.append(cluster_prototypes)

                print("Clara similarities overall {}".format(bestSim))
                print("Clara reaching convergence. Stop clustering.")
                break
            else:
                bestSim = float(s_score)      
                print("sample number : {}".format(i))
                #print("Cluster similarities {}".format(s_score))
                i += 1 

        
        final_clusters = cluster_map_to_item(self.DataSet, bestclusters)
        
        return bestclusters, final_clusters

    