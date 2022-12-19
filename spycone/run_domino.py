import networkx as nx
import os
from collections import defaultdict
from scipy.stats import chi2
import numpy as np
import pandas as pd
from functools import reduce

from .DOMINO.src.core import domino
from .DOMINO.src.core import preprocess_slices as sl
from .clustering import clustering

dir_path = os.path.dirname(os.path.realpath(__file__))

def _combine_pvals(pvals):
        chisqprob = lambda chisq, df: chi2.sf(chisq, df)
        s = -2 * np.sum(np.log(pvals))
        
        return chisqprob(s, 2 * len(pvals))

def _check_nodes(target, network_file):
    if type(network_file)==str:
        G = nx.read_edgelist(network_file)
    else:
        G = network_file.networkz
    
    checkedtarget = []
    for i in target:
        if i in set(G.nodes):
            checkedtarget.append(i)

    return checkedtarget

def _affected_domains(targets, ascov):
    # #get genes switching events
    DDInodes=[]

    for g in ascov.iterrows():
        for domain in g[1]['exclusive_domains']:
            tmp = str(g[1]['gene']+"/"+str(domain))
    
            DDInodes.append(tmp)
    
    return DDInodes
  
            
            
def run_domain_domino(target, is_results, name=None, scores = None, network_file = os.path.join(dir_path,"data/network/network_human_PPIDDI.tab"), output_file_path = "slices.txt", 
                run_cluster=None, slice_threshold=0.3, module_threshold=0.05, prize_factor = 0, n_steps=20):
    '''
    Parameters
    -----------
    target :
        clustering object from spycone or gene list in entrez ID
    is_results : DataFrame
        Data Frame of isoform switch detection result
    scores : None
        activity scores of the genes (e.g. p-values from differential expression analysis)

    run_cluster :
        Specify the cluster name if you only want to run a specific cluster

    Network file : str
        default: "data/network/network_human_PPIDDI.tab"

    output file path
        default: output slices file for DOMINO.

    slice_threshold : float
    
    module_threshold : float

    prize_factor : float 

    n_steps : int
    '''
    if name is None:
        name="a"

    if is_results is not None:
        genescores = defaultdict(list)
        for i in range(is_results.shape[0]):
            genescores[is_results.gene.tolist()[i]].append(is_results.adj_pval.tolist()[i])

        for u,v in genescores.items():
            if len(v)>1:
                genescores[u] = _combine_pvals(v)
            else:
                genescores[u] = v[0]
    else:
        genescores = None
    

    sl.create_slices(network_file, output_file_path)
    if isinstance(target, list):
        ##check list    

        DDInodes = _affected_domains(list(map(str,target)), is_results)
        checkedtarget = _check_nodes(DDInodes, network_file)

        ##TODO scores

        ## clusterobj can also be list
        a = defaultdict(list)
        
        tmp, scores = domino.main(DDInodes, network_file, slices_file=output_file_path, slice_threshold=slice_threshold,  module_threshold=module_threshold, prize_factor=prize_factor, n_steps=n_steps)
        a[name].append(tmp)
        a[name].append(scores)

        return a

    elif isinstance(target, clustering) and run_cluster is not None:
        a = defaultdict(list)
        gene_list = []
        for cc in run_cluster:
            gene_list.append(target.genelist_clusters[cc])
        
        gene_list = reduce(lambda x,y: x+y, gene_list)
        checkedtarget = _check_nodes(list(map(str,gene_list)), network_file)

        DDInodes = _affected_domains(checkedtarget, is_results) ##get affected domains
        tmp, scores = domino.main(DDInodes, network_file, slices_file=output_file_path, slice_threshold=slice_threshold,  module_threshold=module_threshold, prize_factor=prize_factor, n_steps=n_steps)

        a[','.join(map(str,run_cluster))].append(tmp)
        a[','.join(map(str,run_cluster))].append(scores)

        return a

    else:

        a=defaultdict(list)
        for u,v in target.genelist_clusters.items():
            #scoresdf.to_csv("/nfs/home/students/chit/lrz_ticone/domino_emp/{}_cluster{}_mod.csv".format(name, u), index=False)
            ### 
            checkedtarget = _check_nodes(list(map(str,v)), network_file)
            DDInodes = _affected_domains(checkedtarget, is_results)

            scores = []
            for gene in DDInodes:
                if genescores is not None:
                    if gene in genescores.keys():
                        scores.append(genescores[gene])
                    else:
                        scores.append(1)

                else:
                    scores.append(1)
            ##fornow emp
            scoresdf = pd.DataFrame(index=list(map(str,v)), dtype=str)
            ## 
            ##TODO scores

            ddimod, scores=domino.main(DDInodes, network_file=network_file, slices_file=output_file_path, slice_threshold=slice_threshold,  module_threshold=module_threshold, prize_factor=prize_factor, n_steps=n_steps)
            
            a[u].append(ddimod)
            a[u].append(scores)

        print("---------Network enrichment Result---------\n")
        for u,v in a.items():
            for e, vv in enumerate(v[0]):
                print(f"Cluster {u} Module {e} has {len(vv)} nodes.")
        print("-----END-----")        
        return a 

def run_domino(target, name=None, is_results=None, scores = None, network_file = os.path.join(dir_path,"data/network/mouse_biogrid_entrez.tab"), output_file_path = "./slices.txt", 
                run_cluster=None, slice_threshold=0.3, module_threshold=0.05, prize_factor = 0, n_steps=20):
    '''
    Parameters
    -----------
    target :
        clustering object from spycone or gene list in entrez ID
        
    is_results (Optional) : DataFrame 
        Data Frame of isoform switch detection result
        
    scores (Optional) : None
        activity scores of the genes (e.g. p-values from differential expression analysis)

    run_cluster (Optional):
        Specify the cluster name if you only want to run a specific cluster

    Network file: Spycone biological network object or str
        default: "data/network/network_human_PPIDDI.tab"

    output_file_path:
        default: output slices file for DOMINO.
        
    slice_threshold : float
    
    module_threshold : float

    prize_factor : float 

    n_steps : int
    '''
    print("start running DOMINO...")
    if name is None:
        name="a"

    if is_results is not None:
        genescores = defaultdict(list)
        for i in range(is_results.shape[0]):
            genescores[is_results['gene_symb'].tolist()[i]].append(is_results.adj_pval.tolist()[i])

        for u,v in genescores.items():
            if len(v)>1:
                genescores[u] = _combine_pvals(v)
            else:
                genescores[u] = v[0]
    else:
        genescores = None
    

    sl.create_slices(network_file, output_file_path)
    # if isinstance(target, list):
    #     ##check list    
    #     checkedtarget = _check_nodes(list(map(str,target)), network_file)

    #     ## clusterobj can also be list
    #     a = defaultdict(list)

    #     tmp, scores = domino.main(list(map(str,checkedtarget)), network_file, slices_file=output_file_path, slice_threshold=slice_threshold,  module_threshold=module_threshold, prize_factor=prize_factor, n_steps=n_steps)
    #     a[name].append(tmp)
    #     a[name].append(scores)

    #     return a

    # elif isinstance(target, clustering) and run_cluster is not None:
    #     a = defaultdict(list)
    #     gene_list = []
    #     for cc in run_cluster:
    #         gene_list.append(target.genelist_clusters[cc])
        
    #     gene_list = reduce(lambda x,y: x+y, gene_list)
    #     checkedtarget = _check_nodes(list(map(str,gene_list)), network_file)

    #     tmp, scores = domino.main(list(map(str, checkedtarget)), network_file, scores=scores, slices_file=output_file_path, slice_threshold=slice_threshold,  module_threshold=module_threshold, prize_factor=prize_factor, n_steps=n_steps)

    #     a[','.join(map(str,run_cluster))].append(tmp)
    #     a[','.join(map(str,run_cluster))].append(scores)

    #     return a

    # else:

    a=defaultdict(list)
    for u,v in target.genelist_clusters.items():
        scores = []
        for gene in target.symbs_clusters[u]:
            if genescores is not None:
                if gene in genescores.keys():
                    scores.append(genescores[gene])
                else:
                    scores.append(1)

            else:
                scores.append(1)
        ##fornow emp
        ##TODO scores
        
        #scoresdf.to_csv("/nfs/home/students/chit/lrz_ticone/domino_emp/{}_cluster{}_mod.csv".format(name, u), index=False)
        ### 
        checkedtarget = _check_nodes(list(map(str,v)), network_file)
        
        tmp, scores = domino.main(list(map(str,checkedtarget)), network_file,  slices_file=output_file_path, slice_threshold=slice_threshold,  module_threshold=module_threshold, prize_factor=prize_factor, n_steps=n_steps)
        a[u].append(tmp)
        a[u].append(scores)

    print("---------Network enrichment Result---------\n")
    for u,v in a.items():
        #for e, vv in enumerate(v[0]):
        print(f"Cluster {u} found {len(v[0])} module(s).")
    print("-----END-----")        
    return a 
