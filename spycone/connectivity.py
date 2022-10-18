import numpy as np
import seaborn as sns
from collections import defaultdict

from ._connectivity.connectivity_task import get_enrichment_undirected_edges_between_clusters_nodedegrees, get_total_undirected_edge_prob_nodedegrees, get_enrichment_undirected_edges_between_clustering_nodedegrees

class connectivity():
    def __init__(self,
                clusterObj,
                DataSet,
                BioNetwork,
                clusterObj2 = None):
        self.clusterObj = clusterObj
        self.DataSet = DataSet
        self.BioNetwork = BioNetwork
        self.directed_enrichment = None
        self.undirected_enrichment = None
        self.total_undirected_edge_prob_nodedegrees = None
        self.clusterObj2 = clusterObj2


    def calculate_connectivity(self):
        clusteringresults = self.clusterObj.genelist_clusters
        dataset = self.DataSet

        ##precalculate network attributes
        maxdegree = self.BioNetwork._get_maxdegree()
        edgecounts = self.BioNetwork._get_total_undirected_edge_count_nodedegrees()
        connected_nodes = self.BioNetwork._get_connected_nodearray()
        nodedegrees = self.BioNetwork._get_undirected_node_degree()
        total_nodedeg_count = self.BioNetwork._get_total_node_degree_counts()

        if self.clusterObj2 is None:
            print("Calculating edges probabilities between clusters.")
            if self.total_undirected_edge_prob_nodedegrees is None:
                c = get_total_undirected_edge_prob_nodedegrees(edgecounts, maxdegree, total_nodedeg_count)
                self.total_undirected_edge_prob_nodedegrees = c
            else:
                c = self.total_undirected_edge_prob_nodedegrees

            print("Calculating enriched edges between clusters.")
            d = get_enrichment_undirected_edges_between_clusters_nodedegrees(c, clusteringresults, maxdegree, edgecounts, connected_nodes, nodedegrees)
            
            ###normalize the enrichment odd ratio with number of object in both clusters
            newv = defaultdict(float)
            for u,v in d.items():
                cluster1 = int(u.split(",")[0])
                cluster2 = int(u.split(",")[1])
                totaln = len(clusteringresults[cluster1])+len(clusteringresults[cluster2])
                newv[u] = v/totaln

            self.undirected_enrichment = newv
            print("Calculated enrichment of edges between clusters")
        
        else:
            clusteringresults2 = self.clusterObj2.genelist_clusters
            print("Calculating edges probabilities between clusters.")
            if self.total_undirected_edge_prob_nodedegrees is None:
                c = get_total_undirected_edge_prob_nodedegrees(edgecounts, maxdegree, total_nodedeg_count)
                self.total_undirected_edge_prob_nodedegrees = c
            else:
                c = self.total_undirected_edge_prob_nodedegrees

            print("Calculating enriched edges between clusters.")
            d = get_enrichment_undirected_edges_between_clustering_nodedegrees(c, clusteringresults, clusteringresults2, maxdegree, edgecounts, connected_nodes, nodedegrees)
            
            ###normalize the enrichment odd ratio with number of object in both clusters
            # newv = defaultdict(float)
            # for u,v in d.items():
            #     cluster1 = int(u.split(",")[0])
            #     cluster2 = int(u.split(",")[1])
            #     totaln = len(clusteringresults[cluster1])+len(clusteringresults2[cluster2])
            #     newv[u] = v/totaln

            self.undirected_enrichment = d
            print("Calculated enrichment of edges between clusters")

        return d

    def vis_cluster_network(self):
        logodd = self.undirected_enrichment
        prototypes = self.clusterObj._prototype

        #create edge list
        logodd_mat = np.zeros((len(prototypes), len(prototypes)))
        prototype_corr = np.zeros((len(prototypes), len(prototypes)))
        for cs, lodd in logodd.items():
            css = cs.split(",")
            logodd_mat[int(css[0])-1,int(css[1])-1] = lodd
            proto1 = prototypes[int(css[0])]
            proto2 = prototypes[int(css[1])]
            prototype_corr[int(css[0])-1,int(css[1])-1] = np.corrcoef(proto1, proto2)[0,1]
            
        x= ['Cluster {}'.format(x) for x in list(self.clusterObj.genelist_clusters.keys())]
        y= ['Cluster {}'.format(x) for x in list(self.clusterObj.genelist_clusters.keys())]

        #code from matplotlib to create heatmap
        #
        # def heatmap(data, row_labels, col_labels, ax=None,
        #         cbar_kw={}, cbarlabel="", **kwargs):
        #     """
        #     Create a heatmap from a numpy array and two lists of labels.

        #     Parameters
        #     ----------
        #     data
        #         A 2D numpy array of shape (N, M).
        #     row_labels
        #         A list or array of length N with the labels for the rows.
        #     col_labels
        #         A list or array of length M with the labels for the columns.
        #     ax
        #         A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        #         not provided, use current axes or create a new one.  Optional.
        #     cbar_kw
        #         A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
        #     cbarlabel
        #         The label for the colorbar.  Optional.
        #     **kwargs
        #         All other arguments are forwarded to `imshow`.
        #     """

        #     if not ax:
        #         ax = plt.gca()

        #     # Plot the heatmap
        #     im = ax.imshow(data, **kwargs)

        #     # Create colorbar
        #     cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        #     cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

        #     # We want to show all ticks...
        #     ax.set_xticks(np.arange(data.shape[1]))
        #     ax.set_yticks(np.arange(data.shape[0]))
        #     # ... and label them with the respective list entries.
        #     ax.set_xticklabels(col_labels)
        #     ax.set_yticklabels(row_labels)

        #     # Let the horizontal axes labeling appear on top.
        #     ax.tick_params(top=True, bottom=False,
        #                 labeltop=True, labelbottom=False)

        #     # Rotate the tick labels and set their alignment.
        #     plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
        #             rotation_mode="anchor")

        #     # Turn spines off and create white grid.
        #     for edge, spine in ax.spines.items():
        #         spine.set_visible(False)

        #     ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
        #     ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
        #     ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
        #     ax.tick_params(which="minor", bottom=False, left=False)

        #     return im, cbar

        # def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
        #              textcolors=("black", "white"),
        #              threshold=None, **textkw):
        #     """
        #     A function to annotate a heatmap.

        #     Parameters
        #     ----------
        #     im
        #         The AxesImage to be labeled.
        #     data
        #         Data used to annotate.  If None, the image's data is used.  Optional.
        #     valfmt
        #         The format of the annotations inside the heatmap.  This should either
        #         use the string format method, e.g. "$ {x:.2f}", or be a
        #         `matplotlib.ticker.Formatter`.  Optional.
        #     textcolors
        #         A pair of colors.  The first is used for values below a threshold,
        #         the second for those above.  Optional.
        #     threshold
        #         Value in data units according to which the colors from textcolors are
        #         applied.  If None (the default) uses the middle of the colormap as
        #         separation.  Optional.
        #     **kwargs
        #         All other arguments are forwarded to each call to `text` used to create
        #         the text labels.
        #     """

        #     if not isinstance(data, (list, np.ndarray)):
        #         data = im.get_array()

        #     # Normalize the threshold to the images color range.
        #     if threshold is not None:
        #         threshold = im.norm(threshold)
        #     else:
        #         threshold = im.norm(data.max())/2.

        #     # Set default alignment to center, but allow it to be
        #     # overwritten by textkw.
        #     kw = dict(horizontalalignment="center",
        #             verticalalignment="center")
        #     kw.update(textkw)

        #     # Get the formatter in case a string is supplied
        #     if isinstance(valfmt, str):
        #         valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

        #     # Loop over the data and create a `Text` for each "pixel".
        #     # Change the text's color depending on the data.
        #     texts = []
        #     for i in range(data.shape[0]):
        #         for j in range(data.shape[1]):
        #             kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
        #             text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
        #             texts.append(text)

        #     return texts

        # fig, (ax,ax1) = plt.subplots(2,1, figsize=(10,12))

        # im, _ = heatmap(logodd_mat, y, x, ax=ax, cmap="PiYG")
        # texts = annotate_heatmap(im)

        # im, _ = heatmap(prototype_corr, y, x, ax=ax1, cmap="PiYG")
        # texts = annotate_heatmap(im)

        # fig.tight_layout()
        # plt.show()
        colmap = sns.light_palette("seagreen", as_cmap=True)
        sns.clustermap(logodd_mat, xticklabels=x, yticklabels=y, cmap=colmap)
        sns.clustermap(prototype_corr, xticklabels=x, yticklabels=y, cmap='coolwarm')

        return None

    
    


