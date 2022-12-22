import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import mannwhitneyu
from scipy.stats import chi2
import warnings
from functools import reduce
from itertools import combinations
import pickle, os, sys
warnings.simplefilter("ignore")
from collections import defaultdict

from ._util_stat.multipletesting import pvalue_correction
from ._feature import featuresObj 
from ._shuffle import shuffling
# Disable
def _blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def _enablePrint():
    sys.stdout = sys.__stdout__

dir_path = os.path.dirname(os.path.realpath(__file__))

# def _iso_switch_between_arrays(arr1, arr2, orgarr1, orgarr2):
#         ab = np.sign(arr1-arr2)

#         ##look for switch points
#         points=[[x for x in range(ab.shape[1]) if ab[r,x] != np.roll(ab[r], 1)[x] and x!=0] for r in range(ab.shape[0])]

#         #points: switch points of all replicates
#         ##check if all r have switch points
#         #then take the points where all appeared
#         def multiple_intersect(points):
#             return np.unique(list(reduce(lambda x, y: set(x).intersection(set(y)), points)))

#         #calculate differences
#         #[(abs(a[x-1] - a[x]) + abs(b[x-1] -b[x]))/2 for x in points], a,b
#         if any(points):
#             final_sp = multiple_intersect(points) 
#             if not any(final_sp):
#                 final_sp = [list(set(x[0]).intersection(set(x[1]))) for x in combinations(points,2)]
#                 final_sp = np.unique(reduce(lambda x,y: x+y, final_sp))
#                 if not any(final_sp):
#                     final_sp = reduce(lambda x,y : x+y, points)
        

#             final_sp = np.sort(final_sp)
#             # mean of max difference at the time points between 2 arrays for all replicates
#             allsp_diff = [np.mean([np.max([abs(arr1[r,x-1] - arr2[r,x-1]), abs(arr2[r,x] -arr1[r,x])]) for r in range(arr1.shape[0])]) for x in final_sp]
            
#         else:
#             final_sp = []
#             allsp_diff=[]

#         #calculate corr
#         cors = [(1-np.corrcoef(orgarr1[r], orgarr2[r])[0,1])/2 for r in range(arr1.shape[0])]
        
#         #for each switch point, the max switch diff before and after 
#         #allsp : all the before and after switch differences for all switch points
#         return len(final_sp)/arr1.shape[1], allsp_diff, final_sp, np.mean(cors)
class basic_pvalue():
    def __init__(self, testobject, object_type, iso_switch = True, shuffle=None, fitness_scores = None, fitness_scores_two_sided = False, n_permutations=1000):
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
            features_original = featuresObj(self.testobject, timeserieslist, feature_type="clusters")
            #get fitness of clusters
            print("Calculating fitness score.")
            original_fitness_score = features_original.feature_store
        

            #data structure for storing permuted features #[features][permutation]
            permuted_fitness = np.zeros((original_fitness_score.shape[0], n_object))

        # else:
        #     if not isinstance(testobject, clusterObj):
        #         raise ValueError("ClusterObj object is needed.")
            
        if self.object_type == "iso_pairs":
            n_object = len(self.testobject.isoobj.isopairs['major_transcript'])
            timeserieslist = self.testobject.timeserieslist
            features_original = featuresObj(self.testobject, timeserieslist, feature_type="iso_pairs")
            original_fitness_score = features_original.feature_store
            permuted_fitness = np.zeros((original_fitness_score.shape[0], n_object))
 

        batch_smaller=np.zeros((original_fitness_score.shape[0], n_object))
        batch_larger=np.zeros((original_fitness_score.shape[0], n_object))
        batch_equal=np.zeros((original_fitness_score.shape[0], n_object))
        batch_comp=np.zeros((original_fitness_score.shape[0], n_object))

        #shuffle feature and fitness  

        #def _cal_permutations(p):
        for p in self.progressbar(it=range(self.n_permutations), prefix="Permutation:", size=self.n_permutations/10):
            _blockPrint()
  
            shuffleObj = shuffling(timeserieslist=timeserieslist)
            #shuffle1 = np.array(shuffleObj.shuffle_dataset_rowwise(), dtype="double")
            #self.testobject.timeserieslist = shuffle1 ##reassign testobject time series to the shuffled dataset
            
            #calculate features for shuffle obj
            iso = iso_function(self.testobject)
            #corrdf = iso.detect_isoform_switch(filtering=False, min_diff=self.testobject.isoobj.min_diff, corr_cutoff=self.testobject.isoobj.corr_cutoff, event_im_cutoff=self.testobject.isoobj.event_im_cutoff)
            #self.testobject = iso.total_isoform_usage(corrdf)
            #shuffle_clu = clustering.clustering(self.testobject, input_type="isoformusage", algorithm=self.testobject.clusterobj.algorithm, metric=self.testobject.clusterobj.metric, linkage=self.testobject.clusterobj.linkage, n_clusters=self.testobject.clusterobj.n_clusters)


            if self.object_type=="clusters":
                permuted_features = featuresObj(testobject=self.testobject, timeserieslist=timeserieslist, feature_type = "clusters", seed=p)
            #calculate fitness for shuffle obj 
                permuted_fitness = permuted_features.feature_store


            elif self.object_type=="iso_pairs":
                permuted_features = featuresObj(testobject=self.testobject, timeserieslist = timeserieslist, feature_type="shuffle_iso_pairs", seed=p)
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
            _enablePrint()
           

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

class iso_function():
    """
    Isoform level analysis

    Parameters 
    ----------
    DataSet : 
        DataSet object
    Species : 
        Species ID

    Methods
    --------
    detect_isoform_switch :

        Parameters
        ----------
        combine : {'median', 'mean'}
            aggregation methods for replicates. Default="median"

        filtering : (boolean, default=True) 
            if True, low expression genes will be filtered out.

        filter_cutoff : default=2
            expression mean cutoff to filter.

        corr_cutoff : default = 0.7
            minimum correlation of isoform pairs to be included in the output.

        p_val_cutoff : default = 0.05   
            significant p-value cutoff to be included in the output.

        min_diff : default = 0.1
            minimum differences of relative abundance to be included in the output.

        event_im_cutoff : default = 0.1
            minimum event importance to be included in the output.

        adjustp : str {'fdr_bh' (default), 'holm_bonf', 'bonf'}
            Method for multiple testing
            bonf: Bonferroni method
            holm_bonf: holm-bonferroni method
            fdr_bh: Benjamin-hochberg false discovery rate

        n_permutations :
            Number of permutations if permutation test is used.

    total_isoform_usage : 
    
        Parameters
        ----------
        ids_result 
            the result dataframe of isoform switch detection

        norm : `boolean`, default=True
            if True, it normalizes time series matrix to relative abundance to gene expression.

        gene_level : `boolean`, default=True
            if True, it calculates total isoform usage for each gene, otherwise individual isoform usage for each isoform


    Return 
    --------
    Create an instance for isoform switch analysis. 
    
    """
    def __init__(self, dataset):
        self.dataset = dataset
        self.species = dataset.species
    
        if "." in dataset.transcript_id[0]:
            infile = open(os.path.join(dir_path,f"data/network/{dataset.species}_ver_trans_pfam.pkl"),'rb')
            exfile = open(os.path.join(dir_path,f"data/network/exons_{dataset.species}_ensembl.pkl"),'rb')
        else:
            infile = open(os.path.join(dir_path,f"data/network/{dataset.species}_ver_trans_pfam.pkl"),'rb')
            exfile = open(os.path.join(dir_path,f"data/network/exons_{dataset.species}_ensembl.pkl"),'rb')
        self.trans_pfam_dict = pickle.load(infile)
        self.trans_exon_dict = pickle.load(exfile)

        groupdict = defaultdict(lambda: defaultdict(list))
        for idx, sym in enumerate(self.dataset.gene_id):
            if sym != 'n':
                groupdict[sym]['trans_id'].append(self.dataset.transcript_id[idx])
                groupdict[sym]['indices'].append(idx)
                groupdict[sym]['array'].append(self.dataset.timeserieslist[:,idx,:])
        
        #normlize work
        for sym, df in groupdict.items():
            org = df['array'][0]
            for x in range(1,len(df['array'])):
                org = np.hstack([org, df['array'][x]])
            
            df['array'] = org.reshape(self.dataset.reps1, len(df['array']), self.dataset.timepts)
            
            tmp = np.array([df['array'][x]/np.sum(df['array'][x], axis=0) for x in range(df['array'].shape[0])])
            tmp[np.isnan(tmp)] = 1
            groupdict[sym]['normarr'] = tmp

        self.normdf = groupdict

    def _take_interval(self, sp, timepts):
        newsp = np.insert(sp,0,0)
        newsp = np.append(newsp, timepts)
        intervals = []
        for i in range(len(sp)):
            intervals.append(np.min([newsp[i+1]-newsp[i], newsp[i+2]-newsp[i+1]]))
        return intervals

    def _iso_switch_between_arrays(self, arr1, arr2, orgarr1, orgarr2):
        ab = np.sign(arr1-arr2)

        if self.dataset.reps1 ==1:
            points = [x for x in range(1, ab.shape[1]) if (ab[0][x] != np.roll(ab[0], 1)[x])]
            final_sp=points
            allsp_diff =[np.mean([abs(arr1[:,x-1] - arr2[:,x-1]), abs(arr2[:,x] -arr1[:,x])]) for x in final_sp]
            cors = (1-np.corrcoef(orgarr1, orgarr2)[0,1])/2
            iso_ratio=len(final_sp)/arr1.shape[0]
        else:
        ##look for switch points
            points=[[x for x in range(ab.shape[1]) if ab[r,x] != np.roll(ab[r], 1)[x] and x!=0] for r in range(ab.shape[0])]
        #points: switch points of all replicates
        ##check if all r have switch points
        #then take the points where all appeared
            def multiple_intersect(points):
                return np.unique(list(reduce(lambda x, y: set(x).intersection(set(y)), points)))

            #calculate differences
            #[(abs(a[x-1] - a[x]) + abs(b[x-1] -b[x]))/2 for x in points], a,b
            rep_agree = round(self.dataset.reps1*0.6)
            if any(points):
                final_sp = multiple_intersect(points) 
                if not any(final_sp):
                    if rep_agree==1:
                        final_sp=[]
                    else:
                        final_sp = [list(set(x[0]).intersection(set(x[1]))) for x in combinations(points,rep_agree)]
                        final_sp = np.unique(reduce(lambda x,y: x+y, final_sp))
                        if not any(final_sp):
                            final_sp = []
            

                final_sp = np.sort(final_sp)
                # mean of max difference at the time points between 2 arrays for all replicates
                allsp_diff = [np.mean([np.max([abs(arr1[r,x-1] - arr2[r,x-1]), abs(arr2[r,x] -arr1[r,x])]) for r in range(arr1.shape[0])]) for x in final_sp]
                
            else:
                final_sp = []
                allsp_diff=[]

        
        #calculate corr
            cors = [(1-np.corrcoef(orgarr1[r], orgarr2[r])[0,1])/2 for r in range(arr1.shape[0])]
            iso_ratio = len(final_sp)/arr1.shape[1]
        #for each switch point, the max switch diff before and after 
        #allsp : all the before and after switch differences for all switch points
        return iso_ratio, allsp_diff, final_sp, np.mean(cors)


    def _event_enrichness(self, normdf, sp, arr1, arr2):
        ##penalize low expressed genes
        penal1 = np.mean([arr1[r]/np.max(normdf[r],axis=0) for r in range(normdf.shape[0])],axis=0)
        penal2 = np.mean([arr2[r]/np.max(normdf[r],axis=0) for r in range(normdf.shape[0])],axis=0)

        return np.mean(np.append(penal1[sp-1:sp+1],penal2[sp-1:sp+1]))

        
    def _remove_low_expressed(self):
        #return isoform indices
        timedf = pd.concat([pd.DataFrame(self.timeserieslist[x]) for x in range(self.timeserieslist.shape[0])], axis=1)
        geneexp = pd.concat([pd.Series(self.symbs, name='symbol'), timedf],axis=1)

        #filter out those higher than 2 by default (can change)
        lowgenes_idx = np.where(geneexp.groupby('symbol').agg('sum').mean(axis=1) < 2)[0]
        #get the gene names lower than cutoff
        lowgenes = geneexp.iloc[lowgenes_idx, 0] 


        return lowgenes

    def _remove_low_abundance(self, pre_normdf):
        #return isoform indices
    

        #calculate fold change of major isoform with respect to the minor isoforms
        lowexp = np.mean(np.max(pre_normdf, axis=2),axis=0)
        
        #filter out those higher than 2 by default (can change)
        lowgenes = np.where(lowexp<=0.1)[0]
        
        filtered_normdf = np.delete(pre_normdf, lowgenes, axis=1)
        return filtered_normdf, lowgenes

    def _combine_pvals(self, pvals):
        ##no used
        chisqprob = lambda chisq, df: chi2.sf(chisq, df)
        s = -2 * np.sum(np.log(pvals))
        
        return chisqprob(s, 2 * len(pvals))

        
    def _take_interval(self, sp):
        #check interval between switch points
        newsp = np.insert(sp,0,0)
        newsp = np.append(newsp, self.dataset.timepts)
        intervals = []
        for i in range(len(sp)):
            intervals.append(np.min([newsp[i+1]-newsp[i], newsp[i+2]-newsp[i+1]]))
        return intervals

    def _diff_before_after_switch(self, normdf, arr1, arr2, allsp_diff, switch_points):
        tmp_sw =  np.insert(switch_points,0,0)
        tmp_sw = np.append(tmp_sw, self.dataset.timepts)
        if self.dataset.reps1>1:
            #get the intervals between time points
            if any(switch_points):
                #the best sp has the highest mean differences before and after switch points between the two time series 
                bestpval1=0
                bestpval2=0
                allbestp=1
                finaldiff = 0
                finalenrich = 0
                best_switch_point=0
                this_prob = 0
                
                thispt=0 

                for bp in range(1, len(tmp_sw)-1):
                    #if len(switch_points) != bp+1:
                    try:
                        ##if the exp before switch is zero expressions
                        iso1_pval = mannwhitneyu(arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0], arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0])
                        iso2_pval = mannwhitneyu(arr2[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0], arr2[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0])
                    except:
                        iso1_pval=[1,1]
                        iso2_pval=[1,1]
                        
                    thisdiff = np.mean([abs(np.mean(arr1[:,tmp_sw[bp-1]:tmp_sw[bp]]) - np.mean(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]])), abs(np.mean(arr2[:,tmp_sw[bp-1]:tmp_sw[bp]]) - np.mean(arr2[:,tmp_sw[bp]:tmp_sw[bp+1]]))])
                    thisenrich = self._event_enrichness(normdf, tmp_sw[bp], arr1, arr2)
                    iso_prob = ((np.sum(arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0]>arr2[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0]))/arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0].shape[0]) + ((np.sum(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0]<arr2[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0]))/(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0].shape[0]))

                    thispt=tmp_sw[bp]

                    if (iso1_pval[1] < allbestp) & (thisenrich > finalenrich):
                        bestpval1, bestpval2=iso1_pval[1], iso2_pval[1]
                        allbestp=iso1_pval[1]
                        finaldiff = thisdiff
                        finalenrich = thisenrich
                        best_switch_point=tmp_sw[bp]
                        this_prob = iso_prob/2
                    elif (iso2_pval[1] < allbestp) & (thisenrich > finalenrich):
                        bestpval1, bestpval2=iso1_pval[1], iso2_pval[1]
                        allbestp=iso2_pval[1]
                        finaldiff = thisdiff
                        finalenrich = thisenrich
                        best_switch_point=tmp_sw[bp]
                        this_prob = iso_prob/2
                    else:
                        continue

                return finaldiff, bestpval1, bestpval2, best_switch_point, finalenrich, this_prob
            else:
                return 0, 1, 1, 0, 0, 0

        else:             
            #the best sp has the highest mean differences before and after switch points between the two time series 
            if any(allsp_diff):
                bp = np.argmax(allsp_diff)
                finalenrich = self._event_enrichness(normdf, bp, arr1, arr2)
                iso_prob = ((np.sum(arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0]>arr2[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0]))/arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0].shape[0]) + ((np.sum(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0]<arr2[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0]))/(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0].shape[0]))
                return allsp_diff[bp], 1, 1, switch_points[bp], finalenrich, iso_prob
            else:
                return 0,1,1,0,0,0

    def _isoform_differential_usage(self, thisnormdf, thisexp):
        #define major transcript
        #highest mean across time in averaged replicates
    
        #vis
        # for n in range(normdf.shape[1]):
        #     plt.plot(normdf[0,n,:])
        # plt.show()
        
        iso_ratio = np.zeros((thisnormdf.shape[1],thisnormdf.shape[1]))
        iso_diff_value = np.zeros((thisnormdf.shape[1],thisnormdf.shape[1]))
        maj_pval = np.zeros((thisnormdf.shape[1],thisnormdf.shape[1]))
        min_pval = np.zeros((thisnormdf.shape[1],thisnormdf.shape[1]))
        best_switch_point = [[] for _ in range(thisnormdf.shape[1])]
        corr_list = np.zeros((thisnormdf.shape[1],thisnormdf.shape[1]))
        enrichness = np.zeros((thisnormdf.shape[1],thisnormdf.shape[1]))
        iso_prob = np.zeros((thisnormdf.shape[1],thisnormdf.shape[1]))
        all_switch_point = [[] for _ in range(thisnormdf.shape[1])]

        #check for all pairs
        for maj in range(thisnormdf.shape[1]):
            orgarr1 = thisexp[:,maj,:] 
            arr1 = thisnormdf[:,maj,:]
            for x in range(thisnormdf.shape[1]): #each minor transcript
                if x == maj:
                    iso_ratio[maj, x], iso_diff_value[maj, x] = 0, 0
                    best_switch_point[maj].append(0)
                    all_switch_point[maj].append(0)
                    corr_list[maj, x] = 0
                    continue
                arr2 = thisnormdf[:,x,:]
                orgarr2 = thisexp[:,x,:]

                ##finding switch points, correlation
               
                iso_ratio[maj, x], allsp, final_sp, corr_list[maj, x] = self._iso_switch_between_arrays(arr1, arr2, orgarr1, orgarr2)
                ##calculate diff. value, p-values
                iso_diff_value[maj, x], maj_pval[maj, x], min_pval[maj, x], bs, enrichness[maj,x], iso_prob[maj, x] = self._diff_before_after_switch(thisnormdf, arr1, arr2, allsp, final_sp)
                
                best_switch_point[maj].append(bs)
                all_switch_point[maj].append(final_sp)


        return iso_ratio, maj_pval, min_pval, iso_diff_value, best_switch_point, corr_list, enrichness, all_switch_point, iso_prob


    def _correlation_coef(self, thisnormdf, maj):
        m_normdf = np.median(thisnormdf, axis=0)
        
        corr_mat = (1-np.corrcoef(m_normdf))/2


        return corr_mat
    

    def _exclusive_domain(self, iso1, iso2):
        if iso1 in self.trans_pfam_dict.keys():
            pfams1 = self.trans_pfam_dict[iso1]
        else: 
            pfams1 = []

        if iso2 in self.trans_pfam_dict.keys():
            pfams2 = self.trans_pfam_dict[iso2]
        else: 
            pfams2 = []
        
        diff = set(pfams1).symmetric_difference(set(pfams2))
        diff = [x for x in list(diff) if not type(x)==float]
        return diff

    def _diff_exons(self, iso1, iso2, final_sp, isoid):
        if iso1 in self.trans_exon_dict.keys():
            exons1 = self.trans_exon_dict[iso1]
        else: 
            exons1 = []

        if iso2 in self.trans_exon_dict.keys():
            exons2 = self.trans_exon_dict[iso2]
        else: 
            exons2 = []

        ##major or minor isoforms
        ##define isoform_from

        iso1_id = np.where(np.isin(self.normdf[isoid]['trans_id'], iso1))[0][0]
        iso2_id = np.where(np.isin(self.normdf[isoid]['trans_id'], iso2))[0][0]
        iso1_values = np.median(self.normdf[isoid]['array'][:,iso1_id, final_sp-1])
        iso2_values = np.median(self.normdf[isoid]['array'][:,iso2_id, final_sp-1])
        
        if iso1_values > iso2_values:
            majoriso, minoriso = iso1, iso2
            majoriso_exons, minoriso_exons = exons1, exons2
        else:
            majoriso, minoriso = iso2, iso1
            majoriso_exons, minoriso_exons = exons2, exons1

        gain_exons = set(majoriso_exons).difference(set(minoriso_exons))
        loss_exons = set(minoriso_exons).difference(set(majoriso_exons))
        nondiff = set(exons1).intersection(set(exons2))
        gain_exons = [x for x in list(gain_exons) if not type(x)==float]
        loss_exons = [x for x in list(loss_exons) if not type(x)==float]
        nondiff = [x for x in list(nondiff) if not type(x)==float]

        return majoriso, minoriso, gain_exons, loss_exons, nondiff
        
                
    def detect_isoform_switch(self, filtering = True, prob_cutoff=0.5, corr_cutoff=0.5, event_im_cutoff = 0.1, min_diff = 0.05, p_val_cutoff = 0.05, n_permutations=100, adjustp='fdr_bh', diff_domain=True):
        """
        Detect isoform switching events along time series.
        
        

        Return
        ------

        DataFrame 

        """
        #calculate probability and correlation for switching
        #
        #isodf = pd.DataFrame({'symbol':DataSet.symbs, 'isoform_id':DataSet.isoform_id})
        self.min_diff = min_diff
        self.corr_cutoff = corr_cutoff
        self.event_im_cutoff = event_im_cutoff
        self.timeserieslist = self.dataset.timeserieslist
        self.transcript_id = self.dataset.transcript_id
        self.symbs = self.dataset.symbs
        self.gene_id = self.dataset.gene_id
        if filtering==True:
            lowgenes = self._remove_low_expressed()
        else:
            lowgenes=set()
        print(str(len(set(lowgenes))) + " of genes with low expression are filtered out!!.")

        print(f"calculating for {len(self.normdf.keys())} genes...")
              
        iso_pairs_id = defaultdict(list)
        iso_dict_info = defaultdict(tuple)
        for gene, values in self.normdf.items():
            idxs = values['indices']
            thisnormdf = values['normarr']
            thisexp = values['array']

            #calculate isoform ratio
            res_maj = self._isoform_differential_usage(thisnormdf, thisexp)

            #return major transcript and transcript that switch
            if self.dataset.reps1 > 1:
                for ids1, iso_ratio0 in enumerate(res_maj[0]):
                    for ids2, iso_ratio in enumerate(iso_ratio0):
                        if iso_ratio >0 :
                            ##
                            majoriso, minoriso, gain_exons, loss_exons, nondiff = self._diff_exons(self.transcript_id[idxs[ids1]], self.transcript_id[idxs[ids2]], res_maj[4][ids1][ids2], str(gene))
                            iso_pairs_id['gene'].append(str(gene))
                            iso_pairs_id['gene_symb'].append(self.symbs[idxs[ids1]])
                            iso_pairs_id['major_transcript'].append(majoriso)
                            iso_pairs_id['minor_transcript'].append(minoriso)
                            iso_pairs_id['switch_prob'].append(res_maj[8][ids1][ids2])
                            iso_pairs_id['diff'].append(res_maj[3][ids1, ids2])
                            iso_pairs_id['maj_pval'].append(res_maj[1][ids1, ids2])
                            iso_pairs_id['min_pval'].append(res_maj[2][ids1, ids2])
                            iso_pairs_id['p_value'].append(self._combine_pvals([res_maj[1][ids1,ids2], res_maj[2][ids1, ids2]]))
                            iso_pairs_id['corr'].append(res_maj[5][ids1, ids2])
                            iso_pairs_id['best_switch_point'].append(res_maj[4][ids1][ids2])
                            iso_pairs_id['event_importance'].append(res_maj[6][ids1][ids2])
                            iso_pairs_id['gain_exons'].append(gain_exons)
                            iso_pairs_id['loss_exons'].append(loss_exons)
                            iso_pairs_id['constant_exons'].append(nondiff)
                            iso_pairs_id['all_switch_points'].append(res_maj[7][ids1][ids2])

                            iso_dict_info[self.transcript_id[idxs[ids1]]] = (gene, ids1)
                            iso_dict_info[self.transcript_id[idxs[ids2]]] = (gene, ids2)

                            if diff_domain == True:
                                iso_pairs_id['exclusive_domains'].append(self._exclusive_domain(self.transcript_id[idxs[ids1]],self.transcript_id[idxs[ids2]]))


            else:

                for ids1, iso_ratio0 in enumerate(res_maj[0]):
                    for ids2, iso_ratio in enumerate(iso_ratio0):
                        if iso_ratio >0 :
                            ##
                            majoriso, minoriso, gain_exons, loss_exons, nondiff = self._diff_exons(self.transcript_id[idxs[ids1]], self.transcript_id[idxs[ids2]], res_maj[4][ids1][ids2], str(gene))
                            iso_pairs_id['gene'].append(str(gene))
                            iso_pairs_id['gene_symb'].append(self.symbs[idxs[ids1]])
                            iso_pairs_id['major_transcript'].append(majoriso)
                            iso_pairs_id['minor_transcript'].append(minoriso)
                            iso_pairs_id['switch_prob'].append(res_maj[8][ids1][ids2])
                            iso_pairs_id['diff'].append(res_maj[3][ids1, ids2])
                            iso_pairs_id['corr'].append(res_maj[5][ids1, ids2])
                            iso_pairs_id['best_switch_point'].append(res_maj[4][ids1][ids2])
                            iso_pairs_id['event_importance'].append(res_maj[6][ids1][ids2])
                            iso_pairs_id['gain_exons'].append(gain_exons)
                            iso_pairs_id['loss_exons'].append(loss_exons)
                            iso_pairs_id['constant_exons'].append(nondiff)
                            iso_pairs_id['all_switch_points'].append(res_maj[7][ids1][ids2])

                            iso_dict_info[self.transcript_id[idxs[ids1]]] = (gene, ids1)
                            iso_dict_info[self.transcript_id[idxs[ids2]]] = (gene, ids2)

                            if diff_domain == True:
                                iso_pairs_id['exclusive_domains'].append(self._exclusive_domain(self.transcript_id[idxs[ids1]],self.transcript_id[idxs[ids2]]))

        self.isopairs = iso_pairs_id
        self.iso_dict_info = iso_dict_info
        self.dataset.isoobj = self

        if self.dataset.reps1 == 1:
            cal = basic_pvalue(testobject=self.dataset, object_type="iso_pairs", n_permutations = n_permutations)
            pval = cal.do_calculate() #self iso_function object
            iso_pairs_id['p_value'] = pval

            self.isopairs = iso_pairs_id
        
        if pd.DataFrame(iso_pairs_id).shape[0]==0:
            print("We didn't found any significant switch.")
            return None 
        
        res = pd.DataFrame(iso_pairs_id).sort_values('p_value')
        mct = pvalue_correction(res['p_value'], method=adjustp)
        res['adj_pval'] = mct.corrected_pvals

        res = res[(res['switch_prob']>prob_cutoff) & (res['diff']>min_diff) & (res['adj_pval']<p_val_cutoff) & (res['corr']>corr_cutoff) & (res['event_importance']>event_im_cutoff)]
        res = res.sort_values("switch_prob", ascending=False).drop_duplicates(['major_transcript', 'minor_transcript'])
        print("----Result statistics----")
        print(f"Total genes with IS genes: {res['gene'].unique().shape[0]}")
        print(f"Events found: {res.shape[0]}")
        print(f"Events with affecting domains: {np.sum([1 for x in res[  'exclusive_domains'] if len(x)>0])}")
        print(f"Mean event importance: {res['event_importance'].mean()}")
        print(f"Mean difference before and after switch: {res['diff'].mean()}")
        print("--------------------------")
        print("DONE")
        
        self.is_result = res
        self.dataset.isoobj = self

        return res.reset_index(drop=True)[['gene', 'gene_symb', 'major_transcript', 'minor_transcript', 'switch_prob', 'corr', 'diff', 'event_importance', 'exclusive_domains', 'p_value', 'adj_pval']]


    def _isoform_usage(self, eachgene, norm):
        #eachgene 3D array [replicates, isoforms, timepoints]
        result = np.zeros((eachgene.shape[0], eachgene.shape[1], eachgene.shape[2]-1))
        
        if norm:
            #calculate the size factors 
            sizefactor = np.sum(eachgene, axis=1)
            sizefactor[sizefactor==0] = 1

            #normalize to transcript abundance for each replicates
            normdf = np.array([eachgene[x,:,:] / sizefactor[x] for x in range(eachgene.shape[0])])
        
        else:
            normdf = eachgene
        
        #subtract usage from next time point
        for i in range(eachgene.shape[2]-1):
            result[:,:,i] = np.abs(normdf[:,:,i]-normdf[:,:,i+1])
    
        #sum up and take median for one gene across all replicates
        result = np.median(np.sum(result, axis=1),axis=0)
        
        return result
            

    def total_isoform_usage(self, isd_result=None, norm=True, gene_level=True):
        '''
        only genes with more than 1 isoform are computed.

        Parameters
        ----------
        ids_result 
            the result dataframe of isoform switch detection

        norm : `boolean`, default=True
            if True, it normalizes time series matrix to relative abundance to gene expression.

        gene_level : `boolean`, default=True
            if True, it calculates total isoform usage for each gene, otherwise individual isoform usage for each isoform

        Return
        -------
        New DataSet Object
        
        '''
        isodf = pd.DataFrame({'symbol':self.dataset.gene_id, 'isoform_id':self.dataset.transcript_id})
        groupdict = isodf.groupby(['symbol']).indices
        symbolmap = dict(zip(self.dataset.gene_id, self.dataset.symbs))
        tp = self.dataset.timepts
        
        result_df = pd.DataFrame(columns=['entrez','gene'] + ['tp{}'.format(x) for x in range(self.dataset.timepts-1)])
                
        if isd_result is not None:
            keepgenes = isd_result['gene'].unique().tolist()
        else:
            keepgenes = np.unique(self.dataset.gene_id)  

        if gene_level is True:        
            for gene, iso in groupdict.items():
                if len(iso) > 1 and gene in set(keepgenes):
                # iso_switch = isd_result[isd_result['gene']==gene]
                # all_isos = pd.concat([iso_switch.major_isoform, iso_switch.minor_isoform]).unique().tolist()
                # keepiso = np.where(pd.Series(DataSet.isoform_id).isin(all_isos))[0]

                    eachgene = np.array(self.dataset.timeserieslist[:,iso,:], dtype="double")
                    
                    #calculate isoform usage for each gene
                    res = self._isoform_usage(eachgene, norm)
                
                    result_df.loc[result_df.shape[0]] = [gene, symbolmap[gene]] + res.tolist()

            #make changes in DataSet and return new dataset
            self.dataset.tiu_ts = np.array(result_df.iloc[:,2:])
            self.dataset.tiu_timeserieslist = np.expand_dims(np.array(result_df.iloc[:,2:]), axis=0)
            self.dataset.tiu_gene_id = result_df['entrez']
            self.dataset.tiu_symbs = result_df['gene']
            self.dataset.tiu_timepts = tp-1

        
        elif isd_result is not None and gene_level is False:
            keepgenes = isd_result['gene'].unique().tolist()
            keep_genes_arrays = np.hstack([self.normdf[x]['normarr'] for x in keepgenes])

            result = np.zeros((keep_genes_arrays.shape[0], keep_genes_arrays.shape[1], keep_genes_arrays.shape[2]-1))

            keep_trans_indices = reduce(lambda x,y: x+y, [self.normdf[x]['indices'] for x in keepgenes])
            result_trans = self.dataset.transcript_id[keep_trans_indices]
            result_entrez = self.dataset.gene_id[keep_trans_indices]
            result_sym = self.dataset.symbs[keep_trans_indices]

            for i in range(keep_genes_arrays.shape[2]-1):
                result[:,:,i] = np.abs(keep_genes_arrays[:,:,i]-keep_genes_arrays[:,:,i+1])


            self.dataset.tiu_ts = result
            self.dataset.tiu_timeserieslist = result
            self.dataset.tiu_gene_id = result_entrez
            self.dataset.tiu_symbs = result_sym
            self.dataset.tiu_transcript_id = result_trans
            self.dataset.tiu_timepts = tp-1

        return self.dataset


    def isoform_abundance(self, is_result=None):
        """
        Get isoform relative abundance based on isoform switched genes
        """
        if is_result is None:
            keepgenes = is_result['gene'].unique().tolist()
        else:
            keepgenes = self.dataset.gene_id

        result = np.hstack([self.normdf[x]['normarr'] for x in keepgenes])

        keep_trans_indices = reduce(lambda x,y: x+y, [self.normdf[x]['indices'] for x in keepgenes])
        result_trans = self.dataset.transcript_id[keep_trans_indices]
        result_entrez = self.dataset.gene_id[keep_trans_indices]

        tiuds = self.dataset.__copy__()
        tiuds.ts[0] = result
        tiuds.timeserieslist = result
        tiuds.gene_id = result_entrez
        tiuds.transcript_id = result_trans
        tiuds.timepts = self.dataset.timepts
        tiuds.reps1 = self.dataset.reps1

        return tiuds
