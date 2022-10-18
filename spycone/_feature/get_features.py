import numpy as np
import warnings
from joblib import Parallel, delayed

from .._clustering.similarities import intra_clusters_similarity
from .._prototype.prototype_builder import _prototype_builder


##features "number of objects", "similarity value"
class featuresObj():
    def __init__(self, testobject, timeserieslist, feature_type=None, feature_store = None, features=["prototype_standard_variance","average_similarity","number_of_object"], seed=1234):
        self.features = features #prototypesv, averagesim, noofobj
        self.testobject = testobject
        self.timeserieslist = timeserieslist #np.array
        self.feature_store = feature_store ##np.array
        self.feature_type = feature_type #connectivity /clusters / shuffle / iso_pairs / shuffle_iso_pairs
        self.seed = seed
        
        if feature_store is None:
            if feature_type == "clusters" or feature_type == "shuffle":
                try:
                    self.get_clusters_features()
                except:
                    raise ValueError("Clustering object is needed for clusters feature.")

            if feature_type == "iso_pairs":
                #test object data frame of isopairs id
                self.get_switching_features()
            
            elif feature_type == "shuffle_iso_pairs":
                self.get_shuffle_switching_features()

    
    def _cluster_prototype_standard_variance(self):
        if self.feature_type == "clusters":
            prototypes = self.testobject.clusterobj._prototype
            values = []
            for index,i in prototypes.items():
                values.append(np.var(i))
            return values

        if self.feature_type == "shuffle":
            
            clusters = self.testobject.clusterobj.index_clusters
            timeseries_object = self.timeserieslist
            prototype_func = self.testobject.clusterobj.prototypefunction

            prototypes = _prototype_builder(timeseries_object, clusters, prototype_func=prototype_func)
            values = []
            for index,i in prototypes.items():
                values.append(np.var(i))
            return 1/values

    def _cluster_average_similarity(self) -> np.ndarray:
        if self.feature_type == "clusters":
            values = intra_clusters_similarity(self.testobject.timeserieslist, self.testobject.clusterobj.index_clusters, self.testobject.clusterobj._prototype, self.testobject.clusterobj.metric)

            values = np.array(values)
            return values

        if self.feature_type == "shuffle":
            timeseries_object = self.timeserieslist
            clusters = self.testobject.clusterobj.index_clusters

            all_score = intra_clusters_similarity(timeseries_object, clusters, self.testobject.clusterobj._prototype, self.testobject.clusterobj.metric)

            values = all_score
            return values

    def _cluster_information_content(self):
        pass

    def _cluster_number_of_objects(self) -> np.ndarray:
        if self.feature_type == "clusters":
            clusterobj = self.testobject.clusterobj.genelist_clusters
            values = []
            
            for u, v in clusterobj.items():
                values.append(len(v))

            values = np.array(values)
            return values

        if self.feature_type == "shuffle":
            clusterobj = self.clusterobj.genelist_clusters
            values = []
            
            for u, v in clusterobj.items():
                values.append(len(v))

            values = np.array(values)
            return values

    def _clustering_total_obj_cluster_similarity(self):
        pass

    def _clustering_number_of_clusters(self):
        pass

    def _take_interval(self, sp):
            #check interval between switch points
        newsp = np.insert(sp,0,0)
        newsp = np.append(newsp, self.testobject.timepts)
        intervals = []
        for i in range(len(sp)):
            intervals.append(newsp[i+2]-newsp[i+1])
        return intervals

    def _iso_pairs_switch_ratio(self):
        #return iso switching ratio for each pairs
        isopairs = self.testobject.isoobj.isopairs
        nisopairs = len(self.testobject.isoobj.isopairs['major_transcript'])
        iso_ratio = np.zeros(nisopairs)
        iso_dict_info = self.testobject.isoobj.iso_dict_info

        def _cal_iso_ratio(n):
        # for n in range(nisopairs):
            iso1 = iso_dict_info[isopairs['major_transcript'][n]]
            iso2 = iso_dict_info[isopairs['minor_transcript'][n]]
            

            try:
                arr1 = self.testobject.isoobj.normdf[iso1[0]]['normarr'][:,iso1[1],:]
                arr2 = self.testobject.isoobj.normdf[iso2[0]]['normarr'][:,iso2[1],:]

                orgarr1 = self.testobject.isoobj.normdf[iso1[0]]['array'][:,iso1[1],:]
                orgarr2 = self.testobject.isoobj.normdf[iso2[0]]['array'][:,iso2[1],:]
            except:
                raise ValueError("{},{} doesn't seem valid.".format(iso1,iso2))
           
            ir, alldiff, final_sp, cor = self.testobject.isoobj._iso_switch_between_arrays(arr1,arr2,orgarr1,orgarr2)

            if any(alldiff):
                bp = np.argmax(alldiff)
                tmp_sw=np.insert(final_sp, 0,0)
                tmp_sw=np.append(tmp_sw, self.testobject.timepts)
                iso_inter = np.max(self._take_interval(alldiff))
                #iso_prob = ((np.sum(arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0]>arr2[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0]))/arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0].shape[0]) + ((np.sum(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0]<arr2[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0]))/(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0].shape[0]))
                dv = np.max(alldiff)
            else:
                dv = 0
                iso_prob=0
                iso_inter=0
                
            return dv, cor, iso_inter
        #
        #parallel run 
        res = Parallel(backend="loky") (delayed(_cal_iso_ratio)(n) for n in range(nisopairs))
        diff, corr, iso_inter = map(list,zip(*res))
        return corr, iso_inter

    def _shuffle_iso_pairs_switch_ratio(self):
        #return iso switching ratio for each pairs
        isopairs = self.testobject.isoobj.isopairs #id of iso pairs
        nisopairs = len(isopairs['major_transcript'])
        iso_dict_info = self.testobject.isoobj.iso_dict_info

        #make normdf for shuffle
        warnings.simplefilter("ignore")
        from collections import defaultdict
        from functools import partial
        
        groupdict = defaultdict(lambda: defaultdict(list))
        for idx, sym in enumerate(self.testobject.gene_id):
            groupdict[sym]['indices'].append(idx)
            groupdict[sym]['array'].append(self.timeserieslist[:,idx,:])
    
        
        #normlize work
        for sym, df in groupdict.items():
            df['array'] = np.concatenate(df['array']).reshape(self.testobject.reps1,len(df['array']),self.testobject.timepts)
            tmp = df['array']/np.sum(df['array'], axis=1)
            tmp[np.isnan(tmp)] = 0
            for i in range(tmp.shape[0]):
                #np.random.shuffle(tmp[i,:,:])
                for j in range(tmp.shape[1]):
                    np.random.seed(self.seed)
                    np.random.shuffle(tmp[i,j,:])

            groupdict[sym]['normarr'] = tmp

        def _cal_iso_ratio(n):
        # for n in range(nisopairs):
            iso1 = iso_dict_info[isopairs['major_transcript'][n]]
            iso2 = iso_dict_info[isopairs['minor_transcript'][n]]

            arr1 = groupdict[iso1[0]]['normarr'][:,iso1[1],:]
            arr2 = groupdict[iso2[0]]['normarr'][:,iso2[1],:]

            # orgarr1 = groupdict[iso1[0]]['array'][:,iso1[1],:]
            # orgarr2 = groupdict[iso2[0]]['array'][:,iso2[1],:]

            ir, alldiff, final_sp, cor = self.testobject.isoobj._iso_switch_between_arrays(arr1, arr2, arr1, arr2)
            
            if any(alldiff):
                bp = np.argmax(alldiff)
                tmp_sw=np.insert(final_sp, 0,0)
                tmp_sw=np.append(tmp_sw, self.testobject.timepts)
                iso_inter = self._take_interval(alldiff)
                #iso_prob = (((np.sum(arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0]>arr2[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0]))/arr1[:,tmp_sw[bp-1]:tmp_sw[bp]].reshape(1,-1)[0].shape[0]) + ((np.sum(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0]<arr2[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0]))/(arr1[:,tmp_sw[bp]:tmp_sw[bp+1]].reshape(1,-1)[0].shape[0])))/2
                return np.max(alldiff), cor, np.max(iso_inter)
            else:
                return 0, cor, 0

        res = Parallel(backend="loky") (delayed(_cal_iso_ratio)(n) for n in range(nisopairs))

        diff, corr, iso_inter = map(list,zip(*res))
        return corr, iso_inter



    def get_clusters_features(self):
        psd = self._cluster_prototype_standard_variance()
        cas = self._cluster_average_similarity()
        #noo = self.cluster_number_of_objects()

        self.feature_store = np.array([psd, cas])

        return self.feature_store


    def get_switching_features(self):
        #features list correspond to gene id
        ##iso switch ratio
        ##
        fea1 = self._iso_pairs_switch_ratio()

        self.feature_store = np.array(list(fea1))
        return self.feature_store

    def get_shuffle_switching_features(self):
        fea1 = self._shuffle_iso_pairs_switch_ratio()
        
        self.feature_store = np.array(list(fea1))
        return self.feature_store




    