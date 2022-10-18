import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.stats

#from stat.pq_val import _cal_pqvals

def _plot_mean_and_CI(mean, lb, ub, tps, clus, color_mean=None, color_shading=None,  xlab="Time points", ylab="Gene expression"):
    # plot the shaded range of the confidence intervals
    plt.fill_between(range(mean.shape[0]), ub, lb,
                    color=color_shading, alpha=.5)
    # plot the mean on top
    plt.plot(mean, color_mean)
    plt.title('Mean 95% confidence interval of Cluster {}'.format(clus))
    #plt.ylabel("Hyper-similarity score")
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.xticks(np.arange(tps),np.arange(start=1,stop=tps+1,step=1))

def _mean_confidence_interval(data, confidence):
    a = np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

def plotting_clusters(DataSet, confidence=0.95, xlab="Time points", ylab="Gene expression", plot="none"):
    all_g_names = {}
    wilc = []
    t1 = DataSet.ts[1]
    l = DataSet.l[0]
    try:
        t2 = DataSet.group2[1]
    except:
        t2 = None
    object_index = DataSet.object_index[0]
    use_clu = DataSet.use_clu[0]
    symbs1 = DataSet.symbs[0]
    gene_list = DataSet.gene_list[0]
    tps = DataSet.timepts

    if t2 ==None:
        col = ["red", "blue", "green", "yellow", "purple", "white", "brown", "orange", "teal", "silver", "pink"]

        for ii, i in use_clu.items():
            if len(i) != 0:#len(nms):
                if len(i) >= l:
                    tmm = {}
                    for rind in i:
                        tmm[str(gene_list[rind])] = symbs1[rind]
                    all_g_names[str(ii)] = tmm

                    ub1 = []
                    lb1 = []
                    m1 = []
                    for ar in t1[list(i)].T:
                        x,y,z = _mean_confidence_interval(ar, confidence=0.95)
                        m1.append(x)
                        lb1.append(y)
                        ub1.append(z)
                    m1 = np.array(m1)
                    ub1 = np.array(ub1)
                    lb1 = np.array(lb1)

                    wilc.append([t1[list(i)]])
                    # fig = plt.figure(1, figsize=(12, 10))
                    # plot_mean_and_CI_ylim(m1, ub1, lb1, color_mean='red', color_shading='lightcoral')
                    # plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_{}_ylim_conf.pdf".format(j)) 
                    # plt.close()
                    
                    if plot != "all":
                        if type(plot) !="list":
                            plot = [plot]
                        if ii in plot:
                            fig = plt.figure(1, figsize=(12, 10))
                            _plot_mean_and_CI(m1, ub1, lb1, tps, clus=ii, color_mean='red', color_shading='lightcoral', xlab=xlab, ylab=ylab)
                            #plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_{}_conf.png".format(j)) 
                            #plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_{}_conf.pdf".format(j)) 
                            plt.show()
                        else:
                            continue

    
    if t2 != None:
        if len(t2)>0:
            #plt.scatter(d1[:,0], d1[:,1], c="gray")
            col = ["red", "blue", "green", "yellow", "purple", "white", "brown", "orange", "teal", "silver", "pink"]
            j = -1
            for ii, i in use_clu.items():
                if len(i) != 0:#len(nms):
                    #if list(labls).index("3040") in i:
                    #    print("XXXXXXXXXXXXXXXXXXX", i)
                    if len(i) >= l: 
                        j += 1

                        tmm = {}
                        for rind in i:
                            tmm[str(gene_list[rind])] = symbs1[rind]
                        all_g_names[str(ii)] = tmm

                        ub1 = []
                        lb1 = []
                        m1 = []
                        for ar in t1[list(i)].T:
                            x,y,z = _mean_confidence_interval(ar, confidence=0.95)
                            m1.append(x)
                            lb1.append(y)
                            ub1.append(z)
                        m1 = np.array(m1)
                        ub1 = np.array(ub1)
                        lb1 = np.array(lb1)

                        ub2 = []
                        lb2 = []
                        m2 = []

                        for ar in t2[list(i)].T:
                            x,y,z = _mean_confidence_interval(ar, confidence=0.95)
                            m2.append(x)
                            lb2.append(y)
                            ub2.append(z)
                        m2 = np.array(m2)
                        ub2 = np.array(ub2)
                        lb2 = np.array(lb2)

                        wilc.append([t1[list(i)], t2[list(i)]]) #for calculating wilcoxonstats


                        #print([M1[i].flatten(), M2[i].flatten()])
                        
                        if plot != "none":
                            if type(plot)!="list":
                                plot = [plot]
                            if ii in plot:
                                fig = plt.figure(1, figsize=(12, 10))
                                _plot_mean_and_CI(m1, ub1, lb1, tps, clus=ii, color_mean='red', color_shading='lightcoral', xlab=xlab, ylab=ylab)
                                _plot_mean_and_CI(m2, ub2, lb2, tps, clus = ii, color_mean='green', color_shading='lightgreen', xlab=xlab, ylab=ylab)
                                #plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_green_{}_conf.pdf".format(j)) ## png
                                #plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_green_{}_conf.png".format(j)) ## png
                                plt.show()
                            else:
                                continue


    
    DataSet.add_for_graph(all_g_names, wilc)

def _make_simple_graph_plot(gtt1, edg, col, fsize=(40, 40), k=0.2, iterations=90):

        plt.figure(1, figsize=fsize)
        pos = nx.spring_layout(gtt1, seed=42, k=k, iterations=iterations)
        #nx.draw_networkx_nodes(gtt1, pos, nodelist=gtt1.nodes(), node_color="grey",node_size=1000)
        nx.draw_networkx_nodes(gtt1, pos, nodelist=gtt1.nodes(), node_color=col,node_size=1500)


        nx.draw_networkx_edges(gtt1, pos, alpha=0.5, width=3, edge_color="grey")
        nx.draw_networkx_labels(gtt1, pos, font_size=50, font_family='sans-serif', font_color="k")
        plt.axis('off')
        #plt.savefig(name)
        plt.show()

def extract_subnetworks(DataSet, BioNetwork, plot="all"):
    _cal_pqvals(DataSet)

    object_entrez = DataSet.all_g_names[0]
    mtrx = DataSet.mtrx[0]
    object_simval = DataSet.object_simval[0]
    l = DataSet.l[0]
    clu_obj = DataSet.use_clu[0]
    pvals = DataSet.pval[0]

    g = BioNetwork.g()
    mtrx_av = np.average(mtrx)
    print("making graphs")
    #### for bias graph
    
    if plot=="all":
        x=0
        for i,j in object_entrez.items():
            if len(j) >= l:
                print("p-values for cluster {} : {}".format(i, pvals.values()[x]))
                gtmp = g.subgraph(list(j.keys()))
                gtmp = nx.relabel_nodes(gtmp, j)

                _make_simple_graph_plot(gtmp, gtmp.edges(), col="red")
                x+=1

        print("Done making graphs.")
    else:
        if type(plot)!="list":
            plot = [plot]
        for x in plot:
            j = object_entrez[str(x)]
            print("p-values for cluster {} : {}".format(int(x), pvals[x]))
            gtmp = g.subgraph(list(j.keys()))
            gtmp = nx.relabel_nodes(gtmp, j)

            _make_simple_graph_plot(gtmp, gtmp.edges(), col="red")
    