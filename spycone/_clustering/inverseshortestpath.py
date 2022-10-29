import networkx as nx
import pandas as pd
import numpy as np
import itertools as it

#compute pairwise network location similarity
#input as gene node names, corresponding to time series index and network node names


# def inverse_shortest_path(gene_names, BioNetwork):
#     #total node count - 1
#     G = BioNetwork.g()
#     maxdist = len(BioNetwork.lst_g())-1
#     #dist_map = G.new_vertex_property("double", np.inf)
#     mapper = {}
#     for i in range(0, G.num_vertices()):
#         mapper.setdefault(G.vertex_properties['name'][i],i)


#     #initiate
#     #get unique gene list 
#     gene_index = [mapper[x] for x in gene_names]
#     sp_distmat = np.zeros((G.num_vertices(), G.num_vertices()))
    
#     #
    

#     tmp_distmat = shortest_distance(G)


#     #get result from property map
#     for i in range(G.num_vertices()):
#         sp_distmat[i,:] = tmp_distmat[i].a
    

#     #graph tool return a large number when no path is found...
#     #correct those large number so no negative return
#     sp_distmat[sp_distmat>maxdist] = maxdist
#     return (maxdist-sp_distmat)/maxdist, gene_index


