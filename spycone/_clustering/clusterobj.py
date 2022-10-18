import numpy as np

class clusterObj():
    def __init__(self,
                genelist_clusters=None,
                index_clusters = None,
                symbs_clusters = None,
                clusters_sim = None,
                _prototype = None,
                _optimal_cluster = None,
                _final_n_cluster = None):
        self.genelist_clusters = genelist_clusters
        self.symbs_clusters = symbs_clusters
        self.index_clusters = index_clusters
        self.clusters_sim = clusters_sim
        self._prototype = _prototype
        self._optimal_cluster = _optimal_cluster
        self._final_n_cluster = _final_n_cluster

    def _get_indexobjects_in_cluster(self, cl):
        indexclu = self.index_clusters
        objs = indexclu[cl]
        objs = np.array(objs, dtype="object")
        return objs

    def _get_idobjects_in_cluster(self, cl):
        idclu = self.genelist_clusters
        objs = idclu[cl]
        objs = np.array(objs, dtype="object")
        return objs

    def _get_symobjects_in_cluster(self, cl):
        symbsclu = self.symbs_clusters
        objs = symbsclu[cl]
        objs = np.array(objs, dtype="object")
        return objs

    



        


    