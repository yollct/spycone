import pandas as pd
import numpy as np
import time
import random
import warnings
import sys
from scipy.stats import chi2

sys.path.insert(0, "../")
from .._feature import get_features as fea
from .._shuffle import shuffle

from .. import clustering

class basic_pvalue():
    def __init__(self, testobject, object_type, iso_switch = True, shuffle=None, fitness_scores = None, fitness_scores_two_sided = True, n_permutations=1000):
        self.testobject=testobject ##DataSet
        self.object_type = object_type #connectivity, clusters, iso_pairs
        self.shuffle = shuffle
        self.fitness_scores = fitness_scores
        self.fitness_scores_two_sided = fitness_scores_two_sided
        #self.combine_pvalues = combine_pvalues
        self.n_permutations = n_permutations
        self.iso_switch = iso_switch

    def initialize_features_for_permutation(self, features, shuffled):
        pass

    def calculate_original_object_fitness(self):
        pass

    # Print iterations progress
    def progressbar(self, it, prefix="", size=60, file=sys.stdout):
        count = len(it)
        def show(j):
            x = int(size*j/count)
            file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*int((size-x)), j, count))
            file.flush()        
        show(0)
        for i, item in enumerate(it):
            yield item
            show(i+1)
        file.write("\n")
        file.flush()


    def empirical_pvalue(self, total_observation,better_observation):
        if total_observation > 0:
            result = (1+better_observation)/ (1+total_observation)
            return result
        
    def combine_pvals(self, pvals):
        chisqprob = lambda chisq, df: chi2.sf(chisq, df)
        s = -2 * np.sum(np.log(pvals))
        
        return chisqprob(s, 2 * len(pvals))

        
    def do_calculate(self):
        if self.object_type=="clusters":
            n_object = self.testobject.clusterobj._final_n_cluster
            timeserieslist = self.testobject.timeserieslist
            #get features of clusters
            print("Calculating features.")
            features_original = fea.featuresObj(self.testobject, timeserieslist, feature_type="clusters")
            #get fitness of clusters
            print("Calculating fitness score.")
            original_fitness_score = features_original.feature_store
            print(original_fitness_score)

            #data structure for storing permuted features #[features][permutation]
            permuted_fitness = np.zeros((original_fitness_score.shape[0], n_object))

        # else:
        #     if not isinstance(testobject, clusterObj):
        #         raise ValueError("ClusterObj object is needed.")
            
        if self.object_type == "iso_pairs":
            n_object = len(self.testobject.isoobj.isopairs['major_transcript'])
            timeserieslist = self.testobject.timeserieslist
            features_original = fea.featuresObj(self.testobject, timeserieslist, feature_type="iso_pairs")
            original_fitness_score = features_original.feature_store
            permuted_fitness = np.zeros((original_fitness_score.shape[0], n_object))
 

        batch_smaller=np.zeros((original_fitness_score.shape[0], n_object))
        batch_larger=np.zeros((original_fitness_score.shape[0], n_object))
        batch_equal=np.zeros((original_fitness_score.shape[0], n_object))
        batch_comp=np.zeros((original_fitness_score.shape[0], n_object))

        #shuffle feature and fitness  

        #def _cal_permutations(p):
        for p in self.progressbar(it=range(self.n_permutations), prefix="Permutation:", size=self.n_permutations/10):
            shuffleObj = shuffle.shuffling(timeserieslist)
            shuffle1 = np.array(shuffleObj.shuffle_dataset_rowwise(), dtype="double")

            self.testobject.timeserieslist = shuffle1 ##reassign testobject time series to the shuffled dataset
            #calculate features for shuffle obj
            iso = iso_function(self.testobject)
            corrdf = iso.detect_isoform_switch(filtering=False, min_diff=self.testobject.isoobj.min_diff, corr_cutoff=self.testobject.isoobj.corr_cutoff, event_im_cutoff=self.testobject.isoobj.event_im_cutoff)
            self.testobject = iso.total_isoform_usage(corrdf)
            shuffle_clu = clustering.clustering(self.testobject, input_type="isoformusage", algorithm=self.testobject.clusterobj.algorithm, metric=self.testobject.clusterobj.metric, linkage=self.testobject.clusterobj.linkage, n_clusters=self.testobject.clusterobj.n_clusters)


            if self.object_type=="clusters":
                permuted_features = fea.featuresObj(testobject=self.testobject, timeserieslist=shuffle1, feature_type = "clusters")
            #calculate fitness for shuffle obj 
                permuted_fitness = permuted_features.feature_store


            elif self.object_type=="iso_pairs":
                permuted_features = fea.featuresObj(testobject=self.testobject, timeserieslist = shuffle1, feature_type="shuffle_iso_pairs")
                permuted_fitness = permuted_features.feature_store

            
            for fs in range(original_fitness_score.shape[0]):
                fitness = original_fitness_score[fs]
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

        #res = Parallel(backend="multiprocessing") (delayed(_cal_permutations)(int(p)) for p in self.progressbar(it=range(self.n_permutations), prefix="Permutation:", size=self.n_permutations/10))
        #print(res)
       
        #end permutation
        print("Done permutation")
        pvals=np.zeros((original_fitness_score.shape[0], n_object))
        total_observation = 0
        better_observation = 0
        

        for fs in range(original_fitness_score.shape[0]):
            for o in range(n_object):
                total_observation = batch_comp[fs][o]
                

                if self.fitness_scores_two_sided is True:
                    better_observation = min(batch_larger[fs][o], batch_smaller[fs][o])*2 + batch_equal[fs][o]
                else:
                    better_observation = batch_equal[fs][o] + batch_larger[fs][o]

                pvals[fs][o] = self.empirical_pvalue(total_observation, better_observation)
                
        # p_vals = Parallel(backend="multiprocessing") (delayed(_count_permutations)(fs, o) for o in range(n_object) for fs in range(original_fitness_score.shape[0])) 
       
        #p_vals = np.array(p_vals).reshape(n_object, original_fitness_score.shape[0])

        ##combine pvalues
        result_pvalue = np.zeros(n_object)

        for o in range(n_object):
            p_values_list = np.zeros(original_fitness_score.shape[0])
            p_values_list = pvals[:,o]
            result_pvalue[o] = self.combine_pvals(p_values_list)

        return result_pvalue

        
        

