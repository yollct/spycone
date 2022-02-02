import pandas as pd
import numpy as np
import time
import random
import warnings
from numba import jit, njit, prange
import multiprocessing as mp

from inputdata.DataSet import DataSet
from inputdata.BioNetwork import BioNetwork
from similarity.calculate_similarities import _calculate_similarities_timeseriesobj_to_medoid, _calculate_similarities_clusterobject, _calculate_closest_medoid
from prototype.prototype_builder import _prototype_builder
from util_stat.qc import _qc_silhouette
from clustering.build_cluster import cluster_map_to_item

class PAMK_cluster():
    '''
    PAMK clustering 
    '''
    def __init__(self,
                DataSet,
                n_cluster,
                similarity = "euclidean",
                prototype = "median",
                n_start = 100,
                cluster_medoids=None,
                sim_medoids=None,
                cluster_prototype = None,
                until_convergence = True,
                swap_probability = 0.2,
                swap_iteration = 10,
                for_clara = False):
        self.DataSet = DataSet
        self.n_cluster = n_cluster
        self.n_start = n_start
        self.similarity = similarity
        self.prototype = prototype
        self.cluster_medoids= []
        self.sim_medoids = []
        self.cluster_prototype = []
        self.until_convergence = until_convergence
        self.swap_probability = swap_probability
        self.swap_iteration = swap_iteration
        self.for_clara = for_clara

    def _set_initial_medoids(self, timeseriesobjs, n_cluster):
        # set initial medoids 
        # assigning random objects to each medoid
        n_objects = timeseriesobjs.shape[0]
        cluster_medoids = []
        nonmedoidarray = np.array(timeseriesobjs.copy(), dtype="float32")

        chosen_medoid = 0

        while (chosen_medoid < n_cluster):
            random_choice = random.randint(1, nonmedoidarray.shape[0]-1)
            cluster_medoids.append(np.array(nonmedoidarray[random_choice],dtype="float32"))
            
            nonmedoidarray = np.delete(nonmedoidarray,(random_choice), axis=0)
            chosen_medoid+=1
        
        
        #save in object
        try:
            self.cluster_medoids[0]
            self.cluster_medoids[0] = cluster_medoids
        except:
            self.cluster_medoids.append(cluster_medoids)

        return nonmedoidarray
    
    def _assign_to_nearest_medoids(self, timeseriesobjs):
        ##assign objects to the nearest medoids 
        ##get clustering
        ##calculate new medoids
        cluster_medoids = self.cluster_medoids[0]
        clusters = dict([(str(d+1), []) for d in range(len(cluster_medoids))])
        minid = _calculate_similarities_timeseriesobj_to_medoid(timeseriesobjs, cluster_medoids)

        #get index for the closest medoid for each object
        #lenght = object; value = medoid number;
        #max_id = np.argmin(assigned, axis=1)

        for c,d in enumerate(minid):
            clusters[str(d+1)].append(c)

        print('Get {} clusters.'.format(len(clusters)))

        return clusters

    def _create_mapping():
        pass


    
    

    def find_clusters(self, prototype_func="median"):
        try:
            alltimeseriesobjects = self.DataSet.ts[0]
        except:
            alltimeseriesobjects = self.DataSet
        
        if not hasattr(alltimeseriesobjects, "shape"):
            alltimeseriesobjects = np.array(alltimeseriesobjects, dtype="object")

        print("The time series obj has ({})".format(alltimeseriesobjects.shape))
        n_cluster = self.n_cluster
        n_start = self.n_start
        swap_probability = self.swap_probability
        
        print("PAMK clustering initiated.")
        #set initital medoids randomly
        nonmedoids = self._set_initial_medoids(alltimeseriesobjects, n_cluster)
        
        
        
        #start iterate
        n = 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            
            if self.until_convergence:
                n_start=10000

            bestSim = -1
            while n < n_start:

                #assign objs to nearest medoids
                clusters = self._assign_to_nearest_medoids(alltimeseriesobjects)
                
                #update prototype for new clusters
                new_medoids, nonmedoids = _calculate_closest_medoid(alltimeseriesobjects, clusters)
                
                
                try: 
                    self.cluster_medoids[0]
                    self.cluster_medoids[0] = new_medoids
                except:
                    self.cluster_medoids = new_medoids
                
                #similarity between objects
               
                nouse, s_score = _qc_silhouette(alltimeseriesobjects, clusters, l=1)
                

                if bestSim >= s_score:  
                    cluster_prototypes = _prototype_builder(alltimeseriesobjects, n_cluster, clusters, prototype_func)
                    try:
                        self.cluster_prototype[0]
                        self.cluster_prototype[0] = cluster_prototypes
                    except:
                        self.cluster_prototype.append(cluster_prototypes)

                    print("PAMK similarities overall {}".format(bestSim))
                    print("PAMK reaching convergence. Stop clustering.")
                    break
                else:
                    bestSim = float(s_score)
                    cluster_medoids = new_medoids.copy()
                    n += 1 
                    print("iteration: {}".format(n))
                    #print("With similarities overall {}".format(s_score))
                    
                #SWAP
                iteration = 0
                nonmedoidsarray = nonmedoids.copy()
                randomtoswap = np.random.uniform(0,1,size=(len(cluster_medoids)*len(nonmedoidsarray),))
                swapcheck = False

                if (swap_probability<1 & iteration < swap_iteration):
                    for m in range(len(cluster_medoids)):
                        for o in range(len(nonmedoidsarray)):
                            ids = m * len(nonmedoidsarray) + o
                            if randomtoswap[ids] < swap_probability:
                                ##do swap
                                tmp = cluster_medoids[m]
                                cluster_medoids[m] = nonmedoidsarray[o]
                                nonmedoidsarray[o] = tmp
                                swapcheck = True
                                print("Swapped.")
                        
                    swap_iteration += 1

                if swapcheck==False:
                    print("not swapped. continue cluster")
                    continue
            
            if not self.for_clara:
                final_clusters = cluster_map_to_item(self.DataSet, clusters)

                return final_clusters
            
            else:
                return clusters

        