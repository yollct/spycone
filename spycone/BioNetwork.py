
##network class
#import networkx as nx
import numpy as np
import networkx as nx
import os
dir_path = os.path.dirname(os.path.realpath(__file__))

class BioNetwork:
    """Storage of biological network.

    Parameters
    ----------
    path : 
        dir path to the network file or "human" as default human biogrid network


    Attributes
    -----------
    connected_nodes 

    total_node_degree_counts

    total_undirected_edgecount_nodedegree

    max_degree

    all_degree

    Methods
    ----------
    g : 
        generate network from file path if type is `file` or return itself if type is `network`

    lst_g : 
        return list of node names

    adj : 
        return adjcency matrix

    removing_nodes : 
        remove nodes given the index


    """
    def __init__(
        self, 
        path = "human",
        keytype="entrezid",
        **kwargs
        ):
        self.path = path
        self.keytype=keytype

        self.kwargs = {}
        self.kwargs.update((key, value) for key, value in kwargs.items())
        # self.connected_nodes = connected_nodes
        # self.total_node_degree_counts = total_node_degree_counts
        # self.total_undirected_edgecount_nodedegree = total_undirected_edgecount_nodedegree
        # self.max_degree = max_degree
        # self.all_degree = all_degree
        
   

        
        """
        Read the network path and generate networkx object.
        """
        if self.path == "human":
            #print("Reading edge list and removing loops")
            path = os.path.join(dir_path,f'data/network/9606_ddi_weighted_biogrid_entrez.tab')
            g = nx.read_edgelist(path=path, **self.kwargs)
            tmp = 0
            for gene in g:
                if gene in g[gene]:
                    tmp += 1
                    g.remove_edge(gene, gene)
        elif self.path == "mouse":
            path = os.path.join(dir_path,f'data/network/10090_biogrid_entrez.tab')
            g = nx.read_edgelist(path=path, **self.kwargs)
            tmp = 0
            for gene in g:
                if gene in g[gene]:
                    tmp += 1
                    g.remove_edge(gene, gene)
        else:
            g = nx.read_edgelist(path=path, **self.kwargs)
            tmp = 0
            for gene in g:
                if gene in g[gene]:
                    tmp += 1
                    g.remove_edge(gene, gene)
            # #print("Read edge list, network contains {} nodes and {} edges".format(len(g.nodes()), len(g.edges())))
            # #print('Removed', tmp, 'loops')

            #use graph-tool
            # edgelist = []
            # with open(self.path, "r") as f:
            #     for a in f.readlines():
            #         edgelist.append((a[:-1].split("\t")[0],a[:-1].split("\t")[1]))

            # g = Graph(directed=False)

            # props = g.add_edge_list(edgelist, hashed=True)
            # v_prop = g.new_vertex_property("string")
            # for i in range(g.num_vertices()):
            #     v_prop[i] = props[i]

            # g.vertex_properties['name'] = v_prop
            # remove_parallel_edges(g) #remove duplicates edges
        self.networkz = g

    def g(self):
        """
        Return the networkx object.
        """
        return self.networkz
        
    def lst_g(self):
        """
        Return list of gene name in the network.
        """
        g = self.networkz
        lst_g = list(g)
        lst_g.sort()

        # for i in range(g.num_vertices()):
        #     lst_g.append(g.vertex_properties['name'][i])
        #lst_g = [str(x) for x in lst_g]
        return np.array(lst_g)

    def adj(self):
        """
        Return adjacency matrix of the network.
        """
        #print("computing adjacency matrix...")
        g = self.networkz
        lst_g = self.lst_g()
        adj = nx.to_numpy_matrix(g, nodelist=lst_g)
        
        #graph tool
        # g = self.g()
        # lst_g = self.lst_g()
        # a = adjacency(g)
        # adj = a.A #csr matrix to array
        # adj[adj>1] = 1 
        # return np.array(adj)
        return np.array(adj)


    # def read_graph(self):
    #     self.g = self._getgraph(self.path)
    #     # self.lst_g = self._getgene(self.g)
    #     # self.adj = self._getadj(self.g, self.lst_g)

    def _removing_nodes(self, nodes_to_remove):
        #nodes_to_remove = nodes name
        network = self.networkz

        #network.remove_vertex(nodes_to_remove) #graphtool
        
        network.remove_nodes_from(nodes_to_remove)
        self.networkz = network


    def _get_connected_nodearray(self):
        if "self.connected_nodes" not in locals():
            #basically edge list
            network = self.networkz
            connected_nodes = list(network.edges)

            self.connected_nodes = connected_nodes
            return connected_nodes
        else:
            return self.connected_nodes


    def _get_total_node_degree_counts(self):
        # if self.total_node_degree_counts is None:
        #     #return counts of node degree for all nodes
        #     network = self.g()
        #     #return array
        #     nodedegrees = network.degree_property_map("total").get_array()
        #     result = np.empty((max(nodedegrees)+1), dtype="double")
        #     result.fill(0)
        #     for i in nodedegrees:
        #         result[i] += 1
            
        #     self.total_node_degree_counts = result
        #     return result
        # else:
        #     return self.total_node_degree_counts
        if "self.total_node_degree_counts" not in locals():
            #return counts of node degree for all nodes
            network = self.networkz
            nodedegrees = dict(network.degree())
            result = np.empty((max(nodedegrees.values())+1), dtype="double")
            result.fill(0)
            for i in nodedegrees.values():
                result[i] += 1
            
            self.total_node_degree_counts = result
            return result
        else:
            return self.total_node_degree_counts

        

    def _get_undirected_node_degree(self) -> dict:
        network = self.networkz
        #nodedegrees = network.degree_property_map("total").get_array()
        nodedegrees = dict(network.degree())

        self.all_degree = nodedegrees
        return self.all_degree
        
        

    def _get_maxdegree(self):
        #return the max node degree
        if "self.max_degree" not in locals():
            network = self.networkz
            # nodedegrees = network.degree_property_map("total").get_array()
            # maxdegree = max(nodedegrees)
            nodedegrees = dict(network.degree())
            maxdegree = max(nodedegrees.values())

            self.max_degree = maxdegree

            return self.max_degree
        else:
            return self.max_degree

    def _get_total_undirected_edge_count_nodedegrees(self):
        if "self.total_undirected_edgecount_nodedegree" not in locals():
            #return total node degree for edges for undirected graph
            nodedegrees = dict(self.networkz.degree())
            maxdegree = self._get_maxdegree() 
            total_undirected_edgecount_nodedegree = np.empty((maxdegree+1, maxdegree+1), dtype="double")
            total_undirected_edgecount_nodedegree.fill(0)
            connectednodes = self._get_connected_nodearray()
            
            for i in connectednodes:
                node1 = i[0]
                node2 = i[1]
                degree1 = nodedegrees[node1]
                degree2 = nodedegrees[node2]

                mindeg = min(degree1, degree2)
                maxdeg = max(degree1, degree2)
                
                total_undirected_edgecount_nodedegree[mindeg][maxdeg] += 1

            self.total_undirected_edgecount_nodedegree = total_undirected_edgecount_nodedegree
            return total_undirected_edgecount_nodedegree
        else: 
            return self.total_undirected_edgecount_nodedegree



    def _get_total_directed_edge_count_nodedegrees(self):
        try:
            total_directed_edgecount_nodedegree = self.total_directed_edgecount_nodedegree
        
        except:
            pass
        
        return total_directed_edgecount_nodedegree