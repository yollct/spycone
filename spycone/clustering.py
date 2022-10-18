import pandas as pd
import numpy as np
import time
from functools import reduce
from collections import defaultdict
from sklearn.cluster import *
from sklearn.metrics import pairwise_distances, silhouette_score
from tslearn.metrics import cdist_soft_dtw
from tslearn.clustering import TimeSeriesKMeans
from sklearn.metrics import davies_bouldin_score

from ._clustering.kmedoids import KMedoids
from ._clustering.clusterobj import clusterObj
# from ._util_stat import compute_pvalues 

class clustering(clusterObj):
    """
    Clustering object

    Parameters
    -----------
    DataSet
        DataSet object

    BioNetwork (needed if composite is True)
        BioNetwork object

    input_type: str
        clustering expression data put "expression", and clustering total isoform usage put "isoformusage"
    
    algorithm : {'kmeans', 'kmedoids', 'dbscan', 'hierarchical', 'optics'}
        clustering algorithms from sklearn

    composite : `boolean`, default=True
        if True, distance metrics is composited with inverse shortest path
    
    metric : {'euclidean', 'correlation'}
        metrics from sklearn

    linkage : only for 'hierarchical' clustering {default='average', 'complete', 'ward'}

    prototypefunction : {default='median', 'mean'}
        aggregation function for cluster prototypes

    searchspace : (default=20)
        range to search for optimal number of clusters
    
    Attributes
    -----------
    _prototype : dictionary of prototypes for each cluster (keys)
    
    _lables : array with length of object
        The cluster label of each object

    genelist_clusters : dictionary of clusters
        Key and values pair of clustering with entrez ID (gene_list ID)

    index_clusters : dictionary of clusters
        Key and values pair of clustering with indices 

    symbs_clusters : dictionary of clusters
        Key and values pair of clustering with gene symbols 

    _final_n_cluster : 
        number of clusters

    _silhouette_index :
        silhouette index of this clustering
    
    Examples 
    ---------
    

    """
    def __init__(self,
                DataSet,
                algorithm,
                input_type,
                n_clusters=10,
                composite = False,
                BioNetwork=None,
                metric = "euclidean",
                prototypefunction = "median",
                linkage="average",  #scikit learn recommended
                searchspace = 20,
                seed = 1234321,
                transform = None,
                **kwargs
                ):
        self.DataSet = DataSet 
        self.algorithm = algorithm
        self.n_clusters= n_clusters
        self.input_type = input_type ###expression or isoformusage
        self.BioNetwork = BioNetwork
        self.composite = composite
        self.linkage = linkage
        self.metric = metric
        self.prototypefunction = prototypefunction
        self.seed = seed
        self.transform = transform

        self.kwargs = kwargs
        self.params_init = {}
        self.params_init.update((key, value) for key, value in self.kwargs.items())

        if self.algorithm == "kmeans":
            assert self.metric == "euclidean" or self.metric == "soft_dtw", "Kmeans can only use euclidean or soft dtw distances as metric."


    CLUSTERING_ALGORITHMS = {'dbscan':DBSCAN, 'optics':OPTICS, 'kmeans': TimeSeriesKMeans, 'kmedoids':KMedoids, 'hierarchical': AgglomerativeClustering}
    PROTOTYPE = {'median':np.median, 'mean':np.mean}
    TRANSFORMATION = {'log':np.log10}
    
    def find_clusters(self):
        """
        Perform clustering 
        """
        t1 = time.time()
        
        alltimeseriesobjects, gene_id = self.__get_timeserieslist_and_geneids()
        func = self.CLUSTERING_ALGORITHMS[self.algorithm]       
        simfunc = self.PROTOTYPE[self.prototypefunction]

        ###todo###
        #combine replicates
        allrep = simfunc(alltimeseriesobjects, axis=0)
        
        if self.transform is not None:
            transfunc = self.TRANSFORMATION[self.transform] 
            allrep = transfunc(allrep)
            #if value returns infinite , set it as 0
            allrep[~np.isfinite(allrep)]=0

        self.allrep=allrep 

        dist = self.__cal_dist_mat(allrep)

        #compute inverse shortest distance 
        if self.composite:
            sp, gene_index = inverse_shortest_path(gene_id, self.BioNetwork)

            ixgrid = np.ix_(gene_index, gene_index)
            sp_final = sp[ixgrid]

            #composite metrics
            dist = (dist+sp_final)/2

        #distance matrix diagonal should be 0
        np.fill_diagonal(dist, 0) #return none


        bestsim = -np.inf
        if self.algorithm == "hierarchical" and self.linkage != "ward":
            cluster = func(affinity="precomputed", linkage=self.linkage, n_clusters=self.n_clusters).fit(dist)
            #sil = silhouette_score(dist, cluster.labels_, metric="precomputed")

        elif self.algorithm == "hierarchical" and self.linkage=="ward":
            cluster = func(affinity="euclidean", linkage=self.linkage, n_clusters=self.n_clusters).fit(allrep)
            #sil = silhouette_score(allrep, cluster.labels_, metric="euclidean")

        elif self.algorithm == "kmedoids":
            cluster = func(n_clusters=self.n_clusters).fit(dist)
            #sil =  silhouette_score(dist, cluster.labels_, metric="precomputed")
        
        elif self.algorithm == "kmeans" and self.metric == "soft_dtw":
            func = TimeSeriesKMeans
            cluster = func(n_clusters=self.n_clusters, metric="softdtw").fit(allrep)
            #sil = tssilhouette_score(allrep, cluster.labels_, metric="softdtw", n_jobs=5)

        else:
            cluster = func(metric=self.metric, **self.kwargs).fit(allrep)
            #sil = silhouette_score(dist, cluster.labels_, metric="precomputed")

        sil = silhouette_score(allrep, cluster.labels_, metric=self.metric)
        db = davies_bouldin_score(allrep, cluster.labels_)
        
        self._labels = cluster.labels_+1
        self.silhouette_index = sil
        self.davies_bouldin = db
        self.__map_clusters(self._labels)
        self._prototype = {}
        for u,v in self.index_clusters.items():
            self._prototype.update({u: simfunc(allrep[v], axis=0)})

        
        self.DataSet.clusterobj = self
        print(f"clustering took {time.time()-t1}s. ")
        return self.genelist_clusters

    def __cal_dist_mat(self, allrep):
        if self.metric == "soft_dtw":
            dist = cdist_soft_dtw(allrep)
            dist[np.isnan(dist)] = 0
        else:
            #compute euclidean/correlation 
            dist = pairwise_distances(allrep, metric=self.metric)
            dist[np.isnan(dist)] = 0 
        return dist

    def __get_timeserieslist_and_geneids(self):
        if self.input_type == "expression":
            alltimeseriesobjects = np.array(self.DataSet.timeserieslist, dtype="double")
            gene_id = self.DataSet.gene_id
        else:
            try:
                alltimeseriesobjects = np.array(self.DataSet.tiu_timeserieslist, dtype="double")
                gene_id = self.DataSet.tiu_gene_id
            except ValueError:
                print("Please run iso.total_isoform_usage() to calculate total isoform usage before clustering.")
        return(alltimeseriesobjects, gene_id)

    def __map_clusters(self, clusterlabel):
        if self.input_type == "expression":
            genelist = self.DataSet.gene_id
            symbol = self.DataSet.symbs
            geneindex = range(0, len(genelist))
        else:
            genelist = self.DataSet.tiu_gene_id
            symbol = self.DataSet.tiu_symbs
            geneindex = range(0, len(genelist))

        
        #map clusters to gene list dict
        self.genelist_clusters = defaultdict(list)
        self.index_clusters = defaultdict(list)
        self.symbs_clusters = defaultdict(list)
        for idx, cluster in enumerate(clusterlabel):
            self.genelist_clusters[int(cluster)].append(genelist[idx])
            self.index_clusters[int(cluster)].append(geneindex[idx])
            self.symbs_clusters[int(cluster)].append(symbol[idx])

        self._final_n_cluster = len(self.genelist_clusters)

    def combine_clusters(self, clusters_to_combine):
        update_genelist = defaultdict(list)
        update_index = defaultdict(list)
        update_symbs = defaultdict(list)
        for u in self.genelist_clusters.keys():
            if u in clusters_to_combine:
                update_genelist[clusters_to_combine[0]].append(self.genelist_clusters[u])
                update_index[clusters_to_combine[0]].append(self.index_clusters[u])
                update_symbs[clusters_to_combine[0]].append(self.symbs_clusters[u])
            else:
                update_genelist[u].append(self.genelist_clusters[u])
                update_index[u].append(self.index_clusters[u])
                update_symbs[u].append(self.symbs_clusters[u])

        for u in update_genelist.keys():
            update_genelist[u] = reduce(lambda a,y:a+y, update_genelist[u])
            update_index[u] = reduce(lambda a,y: a+y, update_index[u])
            update_symbs[u] = reduce(lambda a,y:a+y, update_symbs[u])

        self.genelist_clusters = update_genelist
        self.index_clusters = update_index
        self.symbs_clusters = update_symbs
        self._prototype = {}
        for u,v in self.index_clusters.items():
            self._prototype.update({u: self.PROTOTYPE[self.prototypefunction](self.allrep[v], axis=0)})

        ## update silhouette score
        new_label=[]
        new_index=[]
        for u,v in self.index_clusters.items():
            for x in v:
                new_label.append(u)
                new_index.append(x)

        new_df = pd.DataFrame({"label": new_label, "ind": new_index})
        new_df.sort_values('ind')
        self.silhouette_index = silhouette_score(self.allrep, new_df['label'].to_list(), metric=self.metric)
        #self.silhouette_index = silhouette_score(self.allrep, cluster.labels_, metric="euclidean")

    def split_clusters(self, clusters_to_split, n_sub=5):
        """
        clusters_to_split: input the number of cluster (int)
        """
        #get the clusters to split
        allTimeSeriesObject, allGeneId = self.__get_timeserieslist_and_geneids()
        thisTimeSeriesObject = allTimeSeriesObject[:,self.index_clusters[clusters_to_split],:]
        
        func = self.CLUSTERING_ALGORITHMS[self.algorithm]       
        simfunc = self.PROTOTYPE[self.prototypefunction]
        
        thisReplicates = simfunc(thisTimeSeriesObject, axis=0)
        dist = self.__cal_dist_mat(thisReplicates)
        np.fill_diagonal(dist, 0) #return none

        ###TODO functionalize this bunch of code (?)
        if self.algorithm == "hierarchical" and self.linkage != "ward":
            cluster = func(affinity="precomputed", linkage=self.linkage, n_clusters=n_sub).fit(dist)
            #sil = silhouette_score(dist, cluster.labels_, metric="precomputed")

        elif self.algorithm == "hierarchical" and self.linkage=="ward":
            cluster = func(affinity="euclidean", linkage=self.linkage, n_clusters=n_sub).fit(thisReplicates)
            #sil = silhouette_score(allrep, cluster.labels_, metric="euclidean")

        elif self.algorithm == "kmedoids":
            cluster = func(n_clusters=n_sub).fit(dist)
            #sil =  silhouette_score(dist, cluster.labels_, metric="precomputed")
        
        elif self.algorithm == "kmeans" and self.metric == "soft_dtw":
            func = TimeSeriesKMeans
            cluster = func(n_clusters=n_sub, metric="softdtw").fit(thisReplicates)
            #sil = tssilhouette_score(allrep, cluster.labels_, metric="softdtw", n_jobs=5)

        else:
            cluster = func(metric=self.metric, **self.kwargs).fit(thisReplicates)
            #sil = silhouette_score(dist, cluster.labels_, metric="precomputed")
        
        maxclu = np.max(list(self.genelist_clusters))
        sub_genelist_clusters = defaultdict(list)
        sub_index_clusters = defaultdict(list)
        sub_symbs_clusters = defaultdict(list)
        for idx, cluster in enumerate(cluster.labels_+1):
            sub_genelist_clusters[int(cluster)].append(self.genelist_clusters[clusters_to_split][idx])
            sub_index_clusters[int(cluster)].append(self.index_clusters[clusters_to_split][idx])
            sub_symbs_clusters[int(cluster)].append(self.symbs_clusters[clusters_to_split][idx])

        ##remove the original cluster
        del self.genelist_clusters[clusters_to_split]
        del self.index_clusters[clusters_to_split]
        del self.symbs_clusters[clusters_to_split]

        ##combine the clusters
        for u,v in sub_genelist_clusters.items():
            for x in range(len(v)):
                self.genelist_clusters[u+maxclu].append(v[x])
                self.index_clusters[u+maxclu].append(sub_index_clusters[u][x])
                self.symbs_clusters[u+maxclu].append(sub_symbs_clusters[u][x])

        #update prototype
        self._prototype = {}
        for u,v in self.index_clusters.items():
            self._prototype.update({u: self.PROTOTYPE[self.prototypefunction](self.allrep[v], axis=0)})

        ## update silhouette score
        new_label=[]
        new_index=[]
        for u,v in self.index_clusters.items():
            for x in v:
                new_label.append(u)
                new_index.append(x)

        new_df = pd.DataFrame({"label": new_label, "ind": new_index})
        new_df.sort_values('ind')
        self.silhouette_index = silhouette_score(self.allrep, new_df['label'].to_list(), metric="euclidean")

    # def calculate_pvalue(self, object_type = "clusters", n_permutations=1000, fitness_scores_two_sided = True):
    #     cal = compute_pvalues.basic_pvalue(testobject=self.DataSet, object_type = object_type, n_permutations = n_permutations, fitness_scores_two_sided = fitness_scores_two_sided)
    #     pval = cal.do_calculate()

    #     self.cluster_pvalues = np.array(pval)
        
    #     return np.array(pval)
        #fitness_scores.append()