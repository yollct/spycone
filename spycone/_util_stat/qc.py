import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.pyplot as plt
from matplotlib import cm, colors

from ..DataSet import DataSet
from ..BioNetwork import BioNetwork



def _qc_silhouette(timeseriesobj, clusters, l=1):
    
    alltimeseries = timeseriesobj.copy()

    # if not hasattr(t1[0], "shape"):
    #     t1 = np.array(t1)
    allrep_s_score = []
    allrep_df_s = []

    for t1 in alltimeseries:
        clu = clusters.copy() ##clusters with index of objects
        
        #get time course expression array for each genes in each clusters
        time_clu = pd.DataFrame(data=None, columns=["TimeCourse", "clusters"])
        for i,j in clu.items():
            if len(j) >= l:
                for jj in j:
                    time_clu = time_clu.append({"TimeCourse":t1[jj], "clusters":str(i)}, ignore_index=True)

        
        time_clu_d = time_clu.to_dict('list')
    
        #calculate sil. for each samples
        with ProcessPoolExecutor() as exe:
            silhouette = exe.submit(silhouette_samples, list(time_clu_d['TimeCourse']), list(time_clu_d['clusters']))

        allrep_s_score.append(np.mean(silhouette.result()))
        time_clu_d['silhouette'] = silhouette.result()
        allrep_df_s.append(time_clu_d)

    return allrep_df_s, np.mean(allrep_s_score)

def _qc_silhouette_within_cluster(alltimeseriesobject, cluster, clusteritem):
    t1 = alltimeseriesobject
    
    if not hasattr(t1, "shape"):
        t1 = np.array(t1)

    #get time course expression array for each genes in each clusters
    time_clu = pd.DataFrame(data=None, columns=["TimeCourse", "clusters"])

    for jj in clusteritem:
        time_clu = time_clu.append({"TimeCourse":t1[jj], "clusters":str(cluster)}, ignore_index=True)

    time_clu_d = time_clu.to_dict('list') #convert to dict
    
    with ProcessPoolExecutor() as exe:
        silhouette = exe.submit(silhouette_samples, list(time_clu_d['TimeCourse']), list(time_clu_d['clusters']))

    time_clu_d['silhouette'] = silhouette
    time_clu = pd.DataFrame.from_dict(time_clu_d) #convert back to df

    #get max index according to sil. value
    maxids = time_clu.loc[time_clu['silhouette'].idxmax()]

    return maxids


def plot_silhouette(DataSet, clusters):
    try:
        t1 = DataSet.ts[0]
    except:
        t1 = DataSet.copy()
        
    time_clu_d, s_score = _qc_silhouette(t1, clusters)
    sil = pd.DataFrame.from_dict(time_clu_d)

    cluster = sil['clusters'].unique()
    y_lower = 10
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(20, 10)
    txt = 'The silhouette score of the cluster is ' + str(s_score) + ' and the number of clusters is ' + str(
        len(cluster))

    for i in range(len(cluster)):
        
        clusteri = sil['clusters'] == cluster[i]
        tmp = sil[clusteri]
        dic = tmp.to_dict('list')

        y_upper = y_lower + tmp.shape[0]

        sortt = sorted(dic['silhouette'])
        color = cm.nipy_spectral(float(i) / len(cluster))
        plt.fill_betweenx(np.arange(y_lower, y_upper), sortt, facecolor=color, edgecolor=color, label=cluster[i])

        y_lower = y_upper + 1000


    fig.text(.5, .01, txt, ha='center')
    ax.set_title('The Silhouette plot')
    ax.set_xlabel('Silhouette coefficient')
    ax.set_ylabel('Clusters')
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)    
    #plt.show()
    
