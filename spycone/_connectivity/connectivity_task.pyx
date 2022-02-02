import pandas as pd
import numpy as np
import time
import random
import warnings
import itertools as it
import math
cimport cython
cimport numpy as np
from cython.parallel import prange

from ..DataSet import DataSet
from ..BioNetwork import BioNetwork
from .._clustering.clusterobj import clusterObj

DTYPE = np.int32
ctypedef np.int_t DTYPE_t

@cython.boundscheck(False)
cdef double _get_theoreticalpossible_undirected_edges_btn_nodes_with_degrees(tuple deg, int maxdeg, np.ndarray[double, ndim=1] total_nodedeg_count):
    #print("Calculate theoretical possible edges between node degrees {}".format(deg))
    '''
    possible number of  edges with two nodes with specific node degrees.
    '''
    #cdef int maxdeg = BioNetwork._get_maxdegree()
    # cdef np.ndarray[double, ndim=2] theor = np.empty((maxdeg+1, maxdeg+1))
    # theor.fill(-1)
    
    #cdef np.ndarray[double, ndim=1] total_nodedeg_count = BioNetwork._get_total_node_degree_counts()
    
    cdef double result = 0 
    cdef double prob = 0
    if deg[0] != deg[1]:
        result = total_nodedeg_count[deg[0]] * total_nodedeg_count[deg[1]]
    else:
        result = (total_nodedeg_count[deg[0]] * total_nodedeg_count[deg[0]] -1 ) / 2 + total_nodedeg_count[deg[0]]
    
    # theor[deg[0]][deg[1]] = result
    
    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[double, ndim=2] get_total_undirected_edge_prob_nodedegrees(np.ndarray[double, ndim=2] edgecounts, int maxdeg, np.ndarray[double, ndim=1] total_nodedeg_count):
    cdef np.ndarray[double, ndim=2] edgeprobs
    #cdef np.ndarray[double, ndim=2] edgecounts
    cdef tuple edgescounts_combi 
    
        
    print("calculating total edge counts for node degrees")
    #edgecounts = BioNetwork._get_total_undirected_edge_count_nodedegrees()   
    print("Node degrees edge counts : {}".format(edgecounts))

    edgeprobs = np.empty((len(edgecounts), len(edgecounts)), dtype="double")
    edgeprobs.fill(0)
    edgecounts_combi = it.product(range(len(edgecounts)), range(len(edgecounts)))

    for deg in edgecounts_combi:
        #this thread pool calculate the prob. of edges between node degrees
        #returns theo: the probability
        if edgecounts[deg[0]][deg[1]] == 0:
            continue
        theo = _get_theoreticalpossible_undirected_edges_btn_nodes_with_degrees(deg, maxdeg, total_nodedeg_count)

        if theo > 0:
            prob = edgecounts[deg[0]][deg[1]]/theo
        else:
            prob = 0      
        
        edgeprobs[deg[0]][deg[1]] = prob
    
    return edgeprobs


def get_total_directed_edge_prob_nodedegrees(self):
    pass

def get_directed_edges_between_clusters(self):
    #connected_nodes = self.BioNetwork.adj() #get adjacency matrix

    pass

@cython.boundscheck(False)
@cython.wraparound(False)
cdef dict _get_undirected_edges_between_clusters(dict clusteringobj,  list connected_nodes):
    #cdef list connected_nodes 
    cdef dict counts
    cdef int sourceCluster 
    cdef int targetCluster 
    cdef tuple countdata
    
    #connected_nodes = BioNetwork.get_connected_nodearray()
    counts = {}        
   

    #dict stores count data
    for k in it.product(clusteringobj.keys(),repeat=2):
        counts["{},{}".format(k[0],k[1])] = 0
    

    #look for edges
    ##loop through adj matrix
    for edge in connected_nodes:
        sourceNode = edge[0]
        targetNode = edge[1]
    
    ##check which cluster (source and target)
        
        for u,v in clusteringobj.items():
            
            if sourceNode != targetNode:
                if sourceNode in v:
                    sourceCluster = u

                if targetNode in v:
                    targetCluster = u
                
                

                    
            else:
                continue

        if (sourceCluster not in clusteringobj.keys()) or (targetCluster not in clusteringobj.keys()):
            continue

        try:

            p1 = "{},{}".format(sourceCluster, targetCluster)

            counts[p1] += 1
        except:

            p2 = "{},{}".format(targetCluster, sourceCluster)
            counts[p2] += 1
        

    return counts    


def get_expected_directed_edges_between_clusters(self):
    pass

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double _get_expected_undirected_edges_between_clusters(total_undirected_edge_prob_nodedegrees, dict clusteringobj, dict nodedegrees, int c1,int c2):
    cdef np.ndarray[double, ndim=2] edgeprob

    edgeprob = total_undirected_edge_prob_nodedegrees
    cdef np.ndarray[np.int_t, ndim=2] deg
    network = BioNetwork
    cdef double prob = 0
    cdef list nodes, nodes1, nodes2 
    #cdef dict nodedegrees = network._get_undirected_node_degree()
    cdef list degree1, degree2
    cdef int mindeg, maxdeg

    degs = np.empty((edgeprob.shape[0]), dtype="int32")
    degs.fill(0)

    for i in range(edgeprob.shape[0]):
        for j in range(edgeprob.shape[1]):
            if edgeprob[i][j]>0:
                np.append(degs, [i,j])

    prob = 0.0

    if c1 == c2:
        nodes = clusteringobj[c1]
        for i in it.product(range(len(nodes)),range(len(nodes))):
            node = nodes[i[0]]
            degree1 = [nodedegrees[node] if node in nodedegrees.keys() else 0]
            node2 = nodes[i[1]]
            degree2 = [nodedegrees[node] if node in nodedegrees.keys() else 0]
        
            mindeg = min(degree1[0], degree2[0])
            maxdeg = max(degree1[0], degree2[0])
            #print("caluculate expected edges between {} and {}".format(c1,c2))
            prob += edgeprob[mindeg][maxdeg]

    else:
        nodes1 = clusteringobj[c1]
        nodes2 = clusteringobj[c2]
    
        for i in it.product(range(len(nodes1)),range(len(nodes2))):
            node = nodes1[i[0]]
            degree1 = [nodedegrees[node] if node in nodedegrees.keys() else 0]
            node2 = nodes2[i[1]]
            degree2 = [nodedegrees[node2] if node2 in nodedegrees.keys() else 0]
            
            mindeg = min(degree1[0], degree2[0])
            maxdeg = max(degree1[0], degree2[0])
        #print("caluculate expected edges between {} and {}".format(c1,c2))
        
            prob += edgeprob[mindeg][maxdeg]
    return prob

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef get_enrichment_directed_edges_between_clusters_nodedegree(char * c1, char * c2):
    pass


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef dict get_enrichment_undirected_edges_between_clusters_nodedegrees(np.ndarray total_undirected_edge_prob_nodedegrees, dict clusteringobj, int maxdegree, np.ndarray[double, ndim=2] edgecounts, list connected_nodes, dict nodedegrees):
    cdef dict counts, undirected_enrich, expectation
    cdef int c1, c2
    cdef double observed, expected
    cdef tuple k

    counts = _get_undirected_edges_between_clusters(clusteringobj, connected_nodes)
    #print("Edge count between clusters: {}".format(counts))
    undirected_enrich = {}
    expectation = {}

    #dict stores count data
    for k in it.combinations(clusteringobj.keys(), 2):
        #print("calculate enrichment for {} and {}.".format(k[0],k[1]))

        undirected_enrich["{},{}".format(k[0],k[1])] = 0
        c1 = k[0] #cluster1
        c2 = k[1] #cluster2
        
        expected = _get_expected_undirected_edges_between_clusters(total_undirected_edge_prob_nodedegrees, clusteringobj, nodedegrees, c1, c2)
        
        observed = counts["{},{}".format(k[0],k[1])]

        expectation["{},{}".format(k[0],k[1])] = expected
        if (expected != 0) & (observed !=0):
            enrichment = math.log(observed / expected) / math.log(2) 
        else:
            enrichment = 0
        
        undirected_enrich["{},{}".format(k[0],k[1])] = enrichment

    return undirected_enrich


##############for two clusterings ##################
cdef dict _get_undirected_edges_between_clustering(dict clusteringobj1, dict clusteringobj2,  list connected_nodes):
    #cdef list connected_nodes 
    cdef dict counts
    cdef int sourceCluster 
    cdef int targetCluster 
    cdef tuple countdata
    
    #connected_nodes = BioNetwork.get_connected_nodearray()
    counts = {}        
   

    #dict stores count data
    for k0 in clusteringobj1.keys():
        for k1 in clusteringobj2.keys():
            counts[f"{k0},{k1}"] = 0
    

    #look for edges
    ##loop through adj matrix
    for edge in connected_nodes:
        sourceNode = edge[0]
        targetNode = edge[1]

        sourceCluster = 0
        targetCluster = 0
    
    ##check which cluster (source and target)
        for u,v in clusteringobj1.items():
            if sourceNode != targetNode:
                if sourceNode in v:
                    sourceCluster = u
            else:
                continue

        for u,v in clusteringobj2.items():
            if sourceNode != targetNode:
                if targetNode in v:
                    targetCluster = u           
            else:
                continue

        if sourceCluster == 0 or targetCluster == 0:
            continue

        try:

            p1 = "{},{}".format(sourceCluster, targetCluster)
            counts[p1] += 1
        except:

            p2 = "{},{}".format(targetCluster, sourceCluster)
            counts[p2] += 1

    return counts   

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double _get_expected_undirected_edges_between_clustering(total_undirected_edge_prob_nodedegrees, dict clusteringobj1, dict clusteringobj2, dict nodedegrees, int c1,int c2):
    cdef np.ndarray[double, ndim=2] edgeprob

    edgeprob = total_undirected_edge_prob_nodedegrees
    cdef np.ndarray[np.int_t, ndim=2] deg
    network = BioNetwork
    cdef double prob = 0
    cdef list nodes, nodes1, nodes2 
    #cdef dict nodedegrees = network._get_undirected_node_degree()
    cdef list degree1, degree2
    cdef int mindeg, maxdeg

    degs = np.empty((edgeprob.shape[0]), dtype="int32")
    degs.fill(0)

    for i in range(edgeprob.shape[0]):
        for j in range(edgeprob.shape[1]):
            if edgeprob[i][j]>0:
                np.append(degs, [i,j])

    prob = 0.0

    if c1 == c2:
        nodes = clusteringobj1[c1]
        for i in it.product(range(len(nodes)),range(len(nodes))):
            node = nodes[i[0]]
            degree1 = [nodedegrees[node] if node in nodedegrees.keys() else 0]
            node2 = nodes[i[1]]
            degree2 = [nodedegrees[node] if node in nodedegrees.keys() else 0]
        
            mindeg = min(degree1[0], degree2[0])
            maxdeg = max(degree1[0], degree2[0])
            #print("caluculate expected edges between {} and {}".format(c1,c2))
            prob += edgeprob[mindeg][maxdeg]

    else:
        nodes1 = clusteringobj1[c1]
        nodes2 = clusteringobj2[c2]
    
        for i in it.product(range(len(nodes1)),range(len(nodes2))):
            node = nodes1[i[0]]
            degree1 = [nodedegrees[node] if node in nodedegrees.keys() else 0]
            node2 = nodes2[i[1]]
            degree2 = [nodedegrees[node2] if node2 in nodedegrees.keys() else 0]
            
            mindeg = min(degree1[0], degree2[0])
            maxdeg = max(degree1[0], degree2[0])
        #print("caluculate expected edges between {} and {}".format(c1,c2))
        
            prob += edgeprob[mindeg][maxdeg]
    return prob

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef dict get_enrichment_undirected_edges_between_clustering_nodedegrees(np.ndarray total_undirected_edge_prob_nodedegrees, dict clusteringobj1, dict clusteringobj2, int maxdegree, np.ndarray[double, ndim=2] edgecounts, list connected_nodes, dict nodedegrees):
    cdef dict counts, undirected_enrich, expectation
    cdef int c1, c2
    cdef double observed, expected
    cdef tuple k

    counts = _get_undirected_edges_between_clustering(clusteringobj1, clusteringobj2, connected_nodes)
    #print("Edge count between clusters: {}".format(counts))
    undirected_enrich = {}
    expectation = {}

    #dict stores count data
    for k0 in clusteringobj1.keys():
        for k1 in clusteringobj2.keys():
        #print("calculate enrichment for {} and {}.".format(k[0],k[1]))
            undirected_enrich[f"{k0},{k1}"] = 0
            c1 = k0 #cluster1
            c2 = k1 #cluster2
            
            expected = _get_expected_undirected_edges_between_clustering(total_undirected_edge_prob_nodedegrees, clusteringobj1, clusteringobj2, nodedegrees, c1, c2)
            
            observed = counts[f"{k0},{k1}"]

            expectation[f"{k0},{k1}"] = expected
            if (expected != 0) & (observed !=0):
                enrichment = math.log(observed / expected) / math.log(2) 
            else:
                enrichment = 0
            
            undirected_enrich[f"{k0},{k1}"] = enrichment

    return undirected_enrich

    
        