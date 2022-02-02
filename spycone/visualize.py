import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from functools import reduce

from .DataSet import DataSet
from .BioNetwork import BioNetwork
from ._clustering.clusterobj import clusterObj

def vis_all_clusters(clusterObj, x_label="time points", y_label="expression", Titles = "Cluster {col_name}", xtickslabels=None, **kwargs):
    """
    Visualize all the clusters with cluster prototype

    Parameters 
    -----------
    clusterObj
        input clustering object with results

    x_label
        x-axis label of the plot
    
    y_label
        y-axis label of the plot
    
    Titles :"Cluster {col_name}"
        titles for each cluster

    """
    if clusterObj.input_type =="transcript_expression":
        alltimeseries = np.array(clusterObj.DataSet.timeserieslist, dtype="double")
        tp = clusterObj.DataSet.timepts
    elif clusterObj.input_type == "gene_expression":
        alltimeseries = np.array(clusterObj.DataSet.genelevel_expression, dtype="double")
        tp = clusterObj.DataSet.timepts
    else:
        try:
            alltimeseries = np.array(clusterObj.DataSet.tiu_timeserieslist, dtype="double")
            tp = clusterObj.DataSet.tiu_timepts
        except ValueError:
            print("Please run iso.total_isoform_usage() to calculate total isoform usage before clustering.")

    median_allts = np.median(alltimeseries, axis=0)

    if clusterObj.transform:
        transfunc = clusterObj.TRANSFORMATION[clusterObj.transform] 
        median_allts = transfunc(median_allts)
        #if value returns infinite , set it as 0
        median_allts[~np.isfinite(median_allts)]=0
        
    

    ##faster probably
    reorder = []
    clusters = []
    for cluster, cluster_obj in clusterObj.index_clusters.items():
        for v in cluster_obj:
            reorder.append(v)
            clusters.append(cluster)

    precluster_tsarray = pd.DataFrame(median_allts)
    precluster_tsarray.columns = ["tp{}".format(x) for x in range(tp)]
    cluster_tsarray = precluster_tsarray.iloc[reorder,:].reset_index(drop=True)
    
    clusters_sns = pd.concat([pd.Series(clusters,name="clusters"),cluster_tsarray], axis=1)
    clusters_sns = pd.melt(clusters_sns, id_vars="clusters",value_vars=cluster_tsarray.columns, value_name="expression",var_name="timepoints", ignore_index=False)
    clusters_sns['gene'] = clusters_sns.index

    #for prototype
    cluster_pro = pd.DataFrame(np.array(list(clusterObj._prototype.values())))
    cluster_pro.columns = ["tp{}".format(x) for x in range(tp)]
    cluster_pro['clusters'] = list(clusterObj.index_clusters.keys())
    cluster_pro = pd.melt(cluster_pro, id_vars='clusters', value_vars=cluster_pro.columns.to_list()[:-1], value_name="expression", var_name="timepoints")
 
    # allsns = pd.concat(cluster_sns_list)
    # sns.relplot(x="timepoints", y="expression",kind="line", hue="cluster", col="cluster",col_wrap=5,  height=3, aspect=.75, linewidth=2.5, palette= "Set2", data=allsns)

    grey = ['lightgray'] * len(clusterObj.genelist_clusters.keys())
    g = sns.relplot(data=clusters_sns, x="timepoints", y="expression", col="clusters", hue="clusters",
    units="gene", kind="line", height=3, aspect=.75,estimator=None, linewidth=1, palette=grey, legend=False, **kwargs)

    g.set_axis_labels(x_label, y_label, fontsize=13)
    g.set_titles(Titles)
    if xtickslabels:
        g.set_xticklabels(labels=xtickslabels,rotation=90, ha="right", fontsize=8)

    ##plot prototypes
    #g = sns.relplot(data=cluster_pro, x="timepoints", y="expression", kind="line", col="clusters",hue="clusters",col_wrap=3, height=3, aspect=.75, linewidth=4, palette="Set2")
    pal = sns.color_palette("dark", len(clusterObj.genelist_clusters.keys()))
    i=0
    for x,ax in g.axes_dict.items():
        subdata = cluster_pro[cluster_pro['clusters']==x]
        
        sns.lineplot(data=subdata, x="timepoints", y="expression", linewidth=4, color=pal[i], ax=ax, legend=False)
        i+=1

    plt.show()



def vis_cluster(clusterObj, cluster):
    alltimeseries = np.array(clusterObj.DataSet.timeserieslist, dtype="double")
    mean_allts = np.mean(alltimeseries, axis=0)

    cluster_obj = clusterObj.index_clusters[cluster]

    cluster_tsarray = pd.DataFrame(mean_allts[cluster_obj,:])
    cluster_tsarray.columns = ["tp{}".format(x) for x in range(clusterObj.DataSet.timepts)]

    cluster_sns = pd.melt(cluster_tsarray, value_vars=cluster_tsarray.columns, value_name="expression", var_name="timepoints")
    cluster_sns['cluster'] = "Cluster {}".format(cluster)

    sns.relplot(x="timepoints", y="expression",kind="line", col_wrap=5,  height=3, aspect=.75, linewidth=2.5, palette= "Set2", data=cluster_sns)

def cluster_heatmap(clusterObj, scalelab="z-score", xlab="Total isoform usage across time", xticklabels =["E17-PN1", "PN1-PN10", "PN10-28","PN28-PN90"],  z_score=None, standard_scale=None, cmap="coolwarm", save_path=None, cluster_label_space = 0.0001):
    """
    Visualize all clusters in a heatmap.

    Parameters
    -----------
    clusterObj
        clusterobj
    
    scalelab
        label for the scale
    
    xlab
        label for x-axis

    xticklabels
        label for x-ticks
    
    z_score : (default:None, 0, 1)
        visualize z_score value rowwise(0) or columnwise(1), cannot be used together with standard_scale

    standard_scale : (default:None, 0, 1)
        visualize standard_scale value rowwise(0) or columnwise(1), cannot be used together with z_score

    cmap : "coolwarm"
        colormap
    
    save_path
        file path for saving the figure.
    
    Return
    -------
    Figure
    
    """
    reorder=[]
    clusters=[]
    for u,v in clusterObj.index_clusters.items():
        for i in v:
            clusters.append(u)
            reorder.append(i)
    pal = sns.color_palette("Set2", 100)
    rowc = dict(zip(list(clusterObj.index_clusters.keys()), pal))
    rowcol = pd.Series(clusters, name='Clusters').map(rowc)
    rowcol.index = pd.DataFrame(clusterObj.DataSet.timeserieslist[0]).iloc[reorder,:].index

    if cmap != "coolwarm":
        colmap = sns.light_palette(cmap, as_cmap=True)

    else:
        colmap = cmap
    cl = sns.clustermap(pd.DataFrame(clusterObj.DataSet.timeserieslist[0]).iloc[reorder,:], standard_scale=standard_scale, z_score=z_score, cmap=colmap, 
                row_cluster=False, 
                col_cluster=False,
                yticklabels=False, row_colors=rowcol, xticklabels=xticklabels)


    cl.ax_heatmap.set(xlabel=xlab)

    ratios=0
    for u,v in clusterObj.index_clusters.items():

        ratios += len(v)*cluster_label_space
        cl.ax_row_dendrogram.annotate("Cluster "+str(u),(0.5,1-ratios))
    plt.xlabel(scalelab)
    plt.show()
    
    if save_path:
        plt.savefig(save_path)

def vis_heatmap(list_genes, dataset, scalelab="z-score", xlab="Total isoform usage across time", xticklabels =None,  z_score=None, standard_scale=None, cmap="coolwarm", save_path=None, cluster_label_space = 0.0001):
    """
    Visualize all clusters in a heatmap.

    Parameters
    -----------
    list_genes
        list of genes ID you want to visualize the expression value
    
    dataset
        a dataset object

    scalelab
        label for the scale
    
    xlab
        label for x-axis

    xticklabels
        label for x-ticks
    
    z_score : (default:None, 0, 1)
        visualize z_score value rowwise(0) or columnwise(1), cannot be used together with standard_scale

    standard_scale : (default:None, 0, 1)
        visualize standard_scale value rowwise(0) or columnwise(1), cannot be used together with z_score

    cmap : "coolwarm"
        colormap
    
    save_path
        file path for saving the figure.
    
    Return
    -------
    Figure
    
    """
    reorder=[]
    clusters=[]
    for u,v in clusterObj.index_clusters.items():
        for i in v:
            clusters.append(u)
            reorder.append(i)
    pal = sns.color_palette("Set2", 100)
    rowc = dict(zip(list(clusterObj.index_clusters.keys()), pal))
    rowcol = pd.Series(clusters, name='Clusters').map(rowc)
    rowcol.index = pd.DataFrame(clusterObj.DataSet.timeserieslist[0]).iloc[reorder,:].index

    targetset = dataset.timeserieslist[:,np.where(np.isin(dataset.gene_id, list_genes))[0],:]
    if cmap != "coolwarm":
        colmap = sns.light_palette(cmap, as_cmap=True)
    else:
        colmap = cmap
        
    cl = sns.clustermap(pd.DataFrame(dataset.timeserieslist[0]).iloc[reorder,:], standard_scale=standard_scale, z_score=z_score, cmap=colmap, 
                row_cluster=False, 
                col_cluster=False,
                yticklabels=False, row_colors=rowcol, xticklabels=xticklabels)


    cl.ax_heatmap.set(xlabel=xlab)

    ratios=0
    for u,v in clusterObj.index_clusters.items():

        ratios += len(v)*cluster_label_space
        cl.ax_row_dendrogram.annotate("Cluster "+str(u),(0.5,1-ratios))
    plt.xlabel(scalelab)
    plt.show()
    
    if save_path:
        plt.savefig(save_path)

def compare_mod(mods , groups, group_type="clusters", figw=3, top=None, figh=None):
        def heatmap(x, y, size, color, figw, figh):
            fig,ax= plt.subplots(figsize=(figw,figh))

            # Mapping from column names to integer coordinates
            x_labels = [v for v in x.unique()]
            y_labels = [v for v in y.unique()]
            x_to_num = {p[1]:p[0]+0.5 for p in enumerate(x_labels)} 
            y_to_num = {p[1]:p[0]+0.5 for p in enumerate(y_labels)} 

            size_scale = 50
            ax.grid()
            #ax.hlines(y=y.map(y_to_num), xmin=0, xmax=x.map(x_to_num), color='black',linewidth=3, alpha =0.8)
            tpc=ax.scatter(
                x=x.map(x_to_num), # Use mapping for x
                y=y.map(y_to_num), # Use mapping for y
                s=size, # Vector of square sizes, proportional to size parameter
                marker='o',# Use square as scatterplot marker
                c=color,
                cmap='cividis'
            )
            cbar = fig.colorbar(tpc)

            # Show column labels on the axes
            ax.set_xticks([x_to_num[v] for v in x_labels])
            ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='right')
            ax.set_yticks([y_to_num[v] for v in y_labels])
            ax.set_yticklabels(y_labels)
            cbar.set_label("-log10(adjusted p-values)")
            plt.xlim([0,len(x.unique())])


        df = pd.DataFrame()

        for e,i in enumerate(mods):
            if isinstance(i, dict):
                for uu,vv in i.items():
                    vv[0]=pd.DataFrame(vv[0])
                    vv[0]['group'] = f"cluster{groups[e]}_module{uu}"
                    vv[0] = vv[0].sort_values(['Adjusted P-value'])
                    if top is not None:
                        df =pd.concat([df, vv[0].iloc[:top,:]])
                    else:
                        df =pd.concat([df, vv[0]])
            else:
                
                if len(i)>0:
                    i[0]['group'] = f"{group_type}{groups[e]}"
                    i[0]=pd.DataFrame(i[0])
                    i[0] = i[0].sort_values(['Adjusted P-value'])
                    if top is not None:
                            df =pd.concat([df, i[0].iloc[:top,:]])
                    else:
                        df =pd.concat([df, i[0]])
            
            #mod_d['group'] = 'TSIS iso'
        if df.shape[0] == 0:
            return None
        else:
            df = df[df['Adjusted P-value']<0.05].sort_values(['group','Adjusted P-value'], ascending=False)
            # if df.shape[0]>20:
            #     df = df.iloc[:21,:]

            heatmap(x=df['group'], y=df['Term'], color=-np.log10(df['Adjusted P-value']), size=-np.log10(df['Adjusted P-value'])*50, figh=df.shape[0]/3, figw=figw)

                
            return df


def switch_plot(gene, DataSet, ascov, xaxis_label=None, all_isoforms=False, relative_abundance=False):
    """
    Switching plot for isoforms / Expression plot for non-switched genes
    if the input gene is not isoform switched, the expression plot will be plotted.
    
    Parameters 
    ----------
    gene : str
        Input gene ID / symbs you would like to plot
    DataSet : DataSet obj
    ascov : DataFrame
        the result dataframe of your isoform switch detection
    xaxis_label : list
        x axis label for the plots
    """
    try:
        entrezid = DataSet.gene_id[np.where(DataSet.symbs==gene)[0]]
    except:
        entrezid = DataSet.gene_id[np.where(DataSet.gene_id==gene)[0]]

    if str(entrezid[0]) in set(ascov['gene']) and all_isoforms==False:
        ncol = 4
        nrow =  2 if len(entrezid) > 4 else len(entrezid) % 4
    
        #fig, axes = plt.subplots(nrows=nrow, ncols=ncol)
        
        for axx, idx in enumerate(np.where(ascov['gene']==str(entrezid[0]))[0]):
            
            iso1 = ascov['major_transcript'].to_list()[idx]
            iso2 = ascov['minor_transcript'].to_list()[idx]

            
            iso1_id = np.where(DataSet.transcript_id == iso1)[0]
            iso2_id = np.where(DataSet.transcript_id == iso2)[0]
            
            
            if relative_abundance:
                val_name="relative abundance"
                iso_dict = DataSet.isoobj.normdf
                try:
                    iso1_ts = pd.DataFrame(iso_dict[entrezid[0]]['normarr'][:,np.where(iso_dict[entrezid[0]]['indices']==iso1_id)[0][0],:])
                    iso2_ts = pd.DataFrame(iso_dict[entrezid[0]]['normarr'][:,np.where(iso_dict[entrezid[0]]['indices']==iso2_id)[0][0],:])
                except:
                    iso1_ts = pd.DataFrame(np.squeeze(iso_dict[entrezid[0]]['normarr'][:,np.where(iso_dict[entrezid[0]]['indices']==iso1_id)[0][0],:], axis=1))
                    iso2_ts = pd.DataFrame(np.squeeze(iso_dict[entrezid[0]]['normarr'][:,np.where(iso_dict[entrezid[0]]['indices']==iso2_id)[0][0],:], axis=1))
            else:
                val_name="expression"
                try:
                    iso1_ts = pd.DataFrame(np.squeeze(DataSet.timeserieslist[:,iso1_id,:], axis=1))
                    iso2_ts = pd.DataFrame(np.squeeze(DataSet.timeserieslist[:,iso2_id,:], axis=1))
                except:
                    iso1_ts = pd.DataFrame(DataSet.timeserieslist[:,iso1_id,:])
                    iso2_ts = pd.DataFrame(DataSet.timeserieslist[:,iso2_id,:])

            iso1_ts = pd.melt(iso1_ts, value_vars=iso1_ts.columns, value_name=val_name, var_name="time point")
            iso2_ts = pd.melt(iso2_ts, value_vars=iso2_ts.columns, value_name=val_name, var_name="time point")
            iso1_ts['isoform'] = iso1
            iso2_ts['isoform'] = iso2

            isos = pd.concat([iso1_ts, iso2_ts], axis=0)

            g=sns.catplot(x="time point", y=val_name, hue="isoform", capsize=.2, kind="point", data=isos)
            plt.title(gene, fontsize=15)
            g.set_xticklabels(labels=xaxis_label, ha="right", fontsize=12)
            sns.catplot(x="time point", y=val_name, hue="isoform", kind="strip", data=isos)
            g.set_xticklabels(labels=xaxis_label, ha="right", fontsize=12)
            plt.title(gene, fontsize=15)
    
    else:
        
        entrezids = np.where(DataSet.gene_id==gene)[0]
        if len(entrezids) == 0:
            entrezids = np.where(DataSet.symbs==gene)[0]

        allts = []
        if relative_abundance:
            iso_dict = DataSet.isoobj.normdf
            val_name = "relative abundance"
            try:
                iso_ts = pd.DataFrame(np.squeeze(iso_dict[entrezid[0]]['normarr'], axis=1))
            except:
                iso_ts = iso_dict[entrezid[0]]['normarr']

            for eachts in range(iso_ts.shape[1]):
                tmp_ts = pd.DataFrame(iso_ts[:,eachts,:])
                tmp_ts = pd.melt(tmp_ts, value_vars=tmp_ts.columns, value_name=val_name, var_name="time point")
                tmp_ts['isoform'] = DataSet.transcript_id[iso_dict[entrezid[0]]['indices'][eachts]]
                allts.append(tmp_ts)
        else: 
            val_name="expression"
            for idx in entrezids:
                try:
                    iso_ts = pd.DataFrame(np.squeeze(DataSet.timeserieslist[:,idx,:], axis=1))
                except:
                    iso_ts = pd.DataFrame(DataSet.timeserieslist[:,idx,:])
                    
                iso_ts = pd.melt(iso_ts, value_vars=iso_ts.columns, value_name=val_name, var_name="time point")
                iso_ts['isoform'] = DataSet.transcript_id[idx]
                allts.append(iso_ts)

        isos = pd.concat(allts, axis=0)
        g=sns.catplot(x="time point", y=val_name, hue="isoform", capsize=.2, kind="point", data=isos)
  
        g.set_xticklabels(labels=xaxis_label, ha="right", fontsize=12)
        plt.title(gene, fontsize=15)

def expression_plot(gene, DataSet):
    entrezids = np.where(DataSet.gene_id==gene)[0]
    
    if len(entrezids) == 0:
        entrezids = np.where(DataSet.symbs==gene)[0]

    allts = []
    for idx in entrezids:
        try:
            iso_ts = pd.DataFrame(np.squeeze(DataSet.timeserieslist[:,idx,:], axis=1))
        except:
            iso_ts = pd.DataFrame(DataSet.timeserieslist[:,idx,:])
            
        iso_ts = pd.melt(iso_ts, value_vars=iso_ts.columns, value_name="expression", var_name="time point")
        allts.append(iso_ts)

    isos = pd.concat(allts, axis=0)
    # sumisos = isos.groupby("time point").agg(sum)
    # sumisos = sumisos.reset_index()

    g= sns.catplot(x="time point", y= "expression", kind="point",data=isos)
    plt.title(gene, fontsize=15)



def gsea_plot(gsea_result, cluster, modules=None, nterms=None):
    """
    Visualizing the functional enrichment

    Parameters
    -----------
    gsea_result : dict 
        the results of gsea
    cluster : str
        the cluster number you would like to visualize
    nterms : int 
        (optional) if you would like to visualize only subset of terms e.g. the top 10 terms 
    """
    cluster = int(cluster)

    if modules is None:
        mod_cluster = gsea_result[cluster]
    else:
        mod_cluster = gsea_result[cluster][modules]

    if isinstance(mod_cluster, dict):
        subset = pd.DataFrame(mod_cluster[0])
    else: 
        if mod_cluster[0].shape[0]>0:
            subset = pd.DataFrame(mod_cluster[0])
            subset = subset.sort_values(['Adjusted P-value'])

        else:
            print("This gsea has no result.")

    if nterms:
        subset = subset.head(nterms)

    plt.figure(figsize=(subset.shape[0],12))
    subset['-log10(adjusted p-value)'] = -np.log10(subset['Adjusted P-value'])
    sns.barplot(y="Term", x="-log10(adjusted p-value)", data=subset, color="salmon")
    if modules is not None:
        plt.title(f"Cluster {cluster} Module {modules}")
    else:
        plt.title(f"Cluster {cluster}")


def vis_modules(mods, dataset=None, size=5, outputpng=None):
    ascov = dataset.isoobj.is_result

    mapping = dict(zip(dataset.gene_id,dataset.symbs))

    def _map_color(x):
        xx = str(x).split("/")[0]
        if xx in set(ascov['gene'].tolist()):
            try:
                domains = np.unique(reduce(lambda x,y: x+y, ascov[ascov['gene_symb']==xx]['exclusive_domains'].tolist()))
            except:
                domains = []

            if len(domains)>0:
                return "#e06666" #pink IS with functional switch
            else:
                return "#fff2cc" ##yellow without functional switch
        else:
            return "#b4f7fe" ##green #interactors

    def _map_entrez2symb(x):
        xx=x.split("/")
        if xx[0] in mapping.keys():
            symb = mapping[xx[0]]
        else: 
            symb = x
        if len(xx)>1:
            return symb+"/"+xx[1]
        else:
            return symb
    
    
    
    for cluster,v in mods.items():
        nmod = np.sum([True if len(eachmod.nodes) > size else False for eachmod in v[0]])#check the how many modules bigger than the threshold size
        nrow = (nmod//3)+1 if nmod/3 - nmod//3 > 0 else nmod//3
        h=nrow*5
        if nrow==0:
            nrow=1
            h=10
        elif nrow==1:
            nrow=2
            h=nrow*5
        
        fig, ax=plt.subplots(nrow, 3, figsize=(h,h))
        if len(ax.shape)==2:
            for e,G in enumerate(zip(v[0], [x for vv in ax for x in vv])):
                if len(G[0].nodes)>size:
                    colors = list(map(_map_color, G[0].nodes))
                    nx.draw(G[0], with_labels=True, labels = dict(zip(G[0].nodes, list(map(_map_entrez2symb, G[0].nodes)))), node_color=colors, node_size=300, font_size=10, ax=G[1])
                # nt.from_nx(G)
                # nt.show('nx.html')
                    G[1].set_title(f"cluster {cluster} module {e} : q-val: {v[1][e]}")
                
            if outputpng:
                plt.savefig(f"{outputpng}/modules_clu{cluster}.png", dpi=200, bbox_inches="tight")
                plt.close()
            else:
                plt.axis('off')
                plt.show()
            

        else:
            for e, G in enumerate(zip(v[0], [x for x in ax])):
                colors = list(map(_map_color, G[0].nodes))
                nx.draw(G[0], with_labels=True, labels = dict(zip(G[0].nodes, list(map(_map_entrez2symb, G[0].nodes)))), node_color=colors, node_size=300, font_size=10, ax=G[1])
                G[1].set_title(f"cluster {cluster} module {e}")

            if outputpng:
                plt.savefig(f"{outputpng}/modules_clu{cluster}.png", dpi=200, bbox_inches="tight")
                plt.close()
            else:
                plt.axis('off')
                plt.show()
            

