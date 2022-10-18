import numpy as np
import random

from collections import defaultdict
from sklearn.cluster import *
from sklearn_extra.cluster import KMedoids
#from similarity.calculate_similarities import calculate_similarities_timeseriesobj_to_medoid, calculate_similarities_clusterobject, calculate_center_medoid
from .._prototype.prototype_builder import _prototype_builder
from .build_cluster import cluster_map_to_item
from .clusterobj import clusterObj

class clara(clusterObj):
    def __init__(self, DataSet, n_cluster, n_start, metric="euclidean",prototypefunction = "median", n_samples=5,
                samplesize=None, genelist_clusters=None, index_clusters=None, symbs_clusters=None, _prototype=None, _final_n_cluster=None):
        self.DataSet= DataSet
        self.n_cluster = n_cluster
        self.n_start = n_start
        self.metric = metric
        self.prototypefunction=prototypefunction
        self.n_samples = n_samples
        self.samplesize = samplesize
        super().__init__(genelist_clusters, index_clusters, symbs_clusters, _prototype, _final_n_cluster)

        def _recommended_samplesize(self):
                n = self.DataSet.ts[0].shape[0]
                k = self.n_cluster

                if n <= 250:
                    return n
                
                return 40+2*k # original ticone is 40+2*k
            
        if samplesize == None:
            self.samplesize = _recommended_samplesize(self)

    def _add_clusters(self, genelistclusters, symbsclusters, indexclusters):
        self.genelist_clusters=genelistclusters
        self.symbs_clusters=symbsclusters
        self.index_clusters=indexclusters

        return None

    def fit(self):
        #return clusterings 

        #start iteration
        bestSim = -np.inf
        bestclusters = defaultdict(list)
        i=0

        while i < self.n_samples:
            #select random samples
            sampledata = self.DataSet.timeserieslist.copy() #which is 3D
            print(sampledata.shape)

            timeseriesobj_ids = np.arange(sampledata.shape[1])
            chosen_index = random.sample(set(timeseriesobj_ids), self.samplesize)

            #fit to k medoids
            sampledata_median = np.median(sampledata[:,chosen_index,:], axis=0)

            print(sampledata_median.shape)

            thisclustering = KMedoids(self.n_cluster)
            thisclustering_fit = thisclustering.fit(sampledata_median)

            #predict non samples
            non_sample_id = np.delete(timeseriesobj_ids, chosen_index)
            nonsampledata = np.median(sampledata[:,non_sample_id,:], axis=0)

            predict_sample = thisclustering.predict(nonsampledata)

            
            #calculate scores
            allindex = np.hstack([chosen_index, non_sample_id])
            alllabels = np.hstack([thisclustering_fit.labels_,predict_sample])
            print(len(chosen_index), non_sample_id.shape)
            print(thisclustering_fit, predict_sample.shape)

            assert allindex.shape[0] == alllabels.shape[0]

            x = np.transpose(np.vstack([alllabels, allindex]))
            
            index_clusters = dict()
            for k,v in x:
                if k in index_clusters.keys():
                    index_clusters[k].append(v)
                else:
                    index_clusters.setdefault(k,[v])
            
            #medoids dict
            thismedoids = dict()
            for u in thisclustering.medoid_indices_:
                thismedoids.setdefault(u,self.DataSet.timeserieslist[:,u,:])


            totalsim = thisclustering.inertia_
            print(totalsim)

            if bestSim <= totalsim:
                bestSim = totalsim
                bestclusters = index_clusters

            i+=1

        genelist_clusters, symbs_clusters = cluster_map_to_item(self.DataSet, bestclusters)
        cluster_prototypes = _prototype_builder(self.DataSet.timeserieslist, bestclusters, self.prototypefunction)
        self._prototype = cluster_prototypes
        self._add_clusters(genelist_clusters, symbs_clusters, bestclusters)
        self._final_n_cluster = len(bestclusters)

        return bestclusters