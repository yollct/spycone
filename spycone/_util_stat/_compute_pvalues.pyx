import pandas as pd
import numpy as np
import time
import random
import warnings

cimport cython
from libc.math cimport fmin
import sys
from joblib import Parallel, delayed

from ..DataSet import DataSet
from ..BioNetwork import BioNetwork
from .._feature import get_features as fea
from .._fitness.fitness_score import fitness_score
from .._shuffle import shuffle
#from ..clustering.clusterobj import clusterObj
from ..iso import iso_function

cdef class basic_pvalue():
    cdef str object_type
    cdef str shuffle
    cdef double[:,:] fitness_scores
    cdef bint fitness_scores_two_sided
    #cdef double[:] combine_pvalues
    cdef int n_permutations
    cpdef public object testobject

    def __cinit__(self, testobject, object_type, shuffle=None, fitness_scores = None, fitness_scores_two_sided = True, n_permutations=1000):
        self.testobject=testobject
        self.object_type = object_type #connectivity, clusters, iso_pairs
        self.shuffle = shuffle
        self.fitness_scores = fitness_scores
        self.fitness_scores_two_sided = fitness_scores_two_sided
        #self.combine_pvalues = combine_pvalues
        self.n_permutations = n_permutations

    cdef initialize_features_for_permutation(self, features, shuffled):
        pass

    cdef calculate_original_object_fitness(self):
        cdef fitness_score
        pass

    # Print iterations progress
    def progressbar(self, it, prefix="", size=60, file=sys.stdout):
        count = len(it)
        def show(j):
            x = int(size*j/count)
            file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
            file.flush()        
        show(0)
        for i, item in enumerate(it):
            yield item
            show(i+1)
        file.write("\n")
        file.flush()


    cdef double empirical_pvalue(self, double total_observation, double better_observation):
        cdef double result = 0
        if total_observation > 0:
            result = (1+better_observation)/ (1+total_observation)
            return result
        
    cdef double combine_pvalues(self, double[:] p_values):
        cdef double result= 1
        cdef double i
        cdef double pseudocount = 1
        for i in p_values:
            if i == 0:
                result *= pseudocount
            else:
                result *= i
        return result
        
    cpdef double[:] do_calculate(self):
        cdef:
            int n_object
            double[:,:,:] timeserieslist
            features_original
            double[:,:] original_fitness_score
            double[:,:] permuted_fitness

        if self.object_type=="clusters":
            n_object = self.testobject._final_n_cluster
            timeserieslist = self.testobject.DataSet.timeserieslist
            #get features of clusters
            print("Calculating features.")
            features_original = fea.featuresObj(self.testobject, timeserieslist, feature_type="clusters")
            #get fitness of clusters
            print("Calculating fitness score.")
            original_fitness_score = features_original.feature_store

            #data structure for storing permuted features #[features][permutation]
            permuted_fitness = np.zeros((original_fitness_score.shape[0], n_object))

        # else:
        #     if not isinstance(testobject, clusterObj):
        #         raise ValueError("ClusterObj object is needed.")
            
        if self.object_type == "iso_pairs":
            n_object = len(self.testobject.isopairs['major'])
            timeserieslist = self.testobject.timeserieslist
            features_original = fea.featuresObj(self.testobject, timeserieslist, feature_type="iso_pairs")
            original_fitness_score = features_original.feature_store
            permuted_fitness = np.zeros((original_fitness_score.shape[0], n_object))
 

        cdef double[:,:] batch_smaller=np.zeros((original_fitness_score.shape[0], n_object))
        cdef double[:,:] batch_larger=np.zeros((original_fitness_score.shape[0], n_object))
        cdef double[:,:] batch_equal=np.zeros((original_fitness_score.shape[0], n_object))
        cdef double[:,:] batch_comp=np.zeros((original_fitness_score.shape[0], n_object))

       
        #shuffle feature and fitness 
        cdef int p=0 #permutation loop
        cdef int fs, co, o
        cdef shuffleObj 
        cdef double[:,:,:] shuffle1   

        def _cal_permutations(p):
        #for p in self.progressbar(it=range(self.n_permutations), prefix="Permutation:", size=self.n_permutations/10):
            shuffleObj = shuffle.shuffling(timeserieslist)
            shuffle1 = np.array(shuffleObj.shuffle_globally(), dtype="double")
            #calculate features for shuffle obj

            if self.object_type=="clusters":
                permuted_features = fea.featuresObj(testobject=self.testobject, timeserieslist=shuffle1, feature_type = "shuffle")
            #calculate fitness for shuffle obj 
                permuted_fitness = permuted_features.feature_store

            elif self.object_type=="iso_pairs":
                permuted_features = fea.featuresObj(testobject=self.testobject, timeserieslist = shuffle1, feature_type="shuffle_iso_pairs")
                permuted_fitness = permuted_features.feature_store


            for fs in range(original_fitness_score.shape[0]):
                fitness = original_fitness_score[fs]
                permuted_fitness = permuted_features.feature_store
                shuffledfitness = permuted_fitness[fs]


                #comparison
                for o in range(len(fitness)):
                    compare = shuffledfitness[o] - fitness[o]
                    if compare < 0:
                        batch_smaller[fs][o] += 1
                    elif compare > 0:
                        batch_larger[fs][o] += 1
                    else:
                        batch_equal[fs][o] += 1
                    batch_comp[fs][o] += 1
            p+=1
            time.sleep(0.1)
            return bs, bl, be, bc

        res = Parallel(backend="multiprocessing") (delayed(_cal_permutations)(p) for p in self.progressbar(it=range(self.n_permutations), prefix="Permutation:", size=self.n_permutations/10))
            
        #print(np.array(batch_smaller), np.array(batch_larger), np.array(batch_equal), np.array(fitness), np.array(shuffledfitness))
        #end permutation
        print("Done permutation")
        cdef double[:,:] pvals=np.zeros((original_fitness_score.shape[0], n_object))
        cdef double total_observation = 0
        cdef double better_observation = 0
        
        def _count_permutations(fs,o):
        # for fs in range(original_fitness_score.shape[0]):
        #     for o in range(n_object):
                total_observation = batch_comp[fs][o]
                

                if self.fitness_scores_two_sided is True:
                    better_observation = fmin(batch_larger[fs][o], batch_smaller[fs][o])*2 + batch_equal[fs][o]
                else:
                    better_observation = batch_equal[fs][o] + batch_larger[fs][o]

                pv = self.empirical_pvalue(total_observation, better_observation)
                return pv
        
        p_vals = Parallel(backend="multiprocessing") (delayed(_count_permutations)(fs, o) for o in range(n_object) for fs in range(original_fitness_score.shape[0])) 

        p_vals = np.array(p_vals).reshape(n_object, original_fitness_score.shape[0])


        ##combine pvalues
        cdef double[:] result_pvalue = np.zeros(n_object)
        cdef double[:] p_values_list

        for o in range(n_object):
            p_values_list = np.zeros(original_fitness_score.shape[0])
            p_values_list = pvals[:,o]
            result_pvalue[o] = self.combine_pvalues(p_values_list)
            
        return result_pvalue

        
        

