import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import matplotlib.colors as mcolors
import random

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

def plotting_clusters(clusterObj, confidence=0.95, xlab="Time points", ylab="Gene expression", plot="none"):
    all_g_names = {}
    wilc = []
    t1 = clusterObj.DataSet.timeserieslist[0]
    l = clusterObj.DataSet.l
    index_cluster = clusterObj.index_clusters[0]
    symbs1 = clusterObj.DataSet.symbs[0]
    gene_list = clusterObj.DataSet.gene_list[0]
    tps = clusterObj.DataSet.timepts

    
    col = list(mcolors.CSS4_COLORS.values())

    for u,v in index_cluster.items():
        if len(v) != 0:#len(nms):
            tmm = {}
            for rind in v:
                tmm[str(gene_list[rind])] = symbs1[rind]
            all_g_names[str(u)] = tmm

            ub1 = []
            lb1 = []
            m1 = []
            for ar in t1[list(v)].T:
                x,y,z = _mean_confidence_interval(ar, confidence=0.95)
                m1.append(x)
                lb1.append(y)
                ub1.append(z)
            m1 = np.array(m1)
            ub1 = np.array(ub1)
            lb1 = np.array(lb1)

            wilc.append([t1[list(v)]])
            # fig = plt.figure(1, figsize=(12, 10))
            # plot_mean_and_CI_ylim(m1, ub1, lb1, color_mean='red', color_shading='lightcoral')
            # plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_{}_ylim_conf.pdf".format(j)) 
            # plt.close()
            
            if plot != "all":
                if type(plot) !="list":
                    plot = [plot]
                if u in plot:
                    fig = plt.figure(1, figsize=(12, 10))
                    figshow = _plot_mean_and_CI(m1, ub1, lb1, tps, clus=u, color_mean=random.choice(col), color_shading='lightgrey', xlab=xlab, ylab=ylab)
                    #plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_{}_conf.png".format(j)) 
                    #plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_{}_conf.pdf".format(j)) 
                    print(figshow)
                else:
                    continue

    
    
                # #plt.scatter(d1[:,0], d1[:,1], c="gray")
                # col = ["red", "blue", "green", "yellow", "purple", "white", "brown", "orange", "teal", "silver", "pink"]
                # j = -1
                # for ii, i in use_clu.items():
                #     if len(i) != 0:#len(nms):
                #         #if list(labls).index("3040") in i:
                #         #    print("XXXXXXXXXXXXXXXXXXX", i)
                #         if len(i) >= l: 
                #             j += 1

                #             tmm = {}
                #             for rind in i:
                #                 tmm[str(gene_list[rind])] = symbs1[rind]
                #             all_g_names[str(ii)] = tmm

                #             ub1 = []
                #             lb1 = []
                #             m1 = []
                #             for ar in t1[list(i)].T:
                #                 x,y,z = _mean_confidence_interval(ar, confidence=0.95)
                #                 m1.append(x)
                #                 lb1.append(y)
                #                 ub1.append(z)
                #             m1 = np.array(m1)
                #             ub1 = np.array(ub1)
                #             lb1 = np.array(lb1)

                #             ub2 = []
                #             lb2 = []
                #             m2 = []

                #             for ar in t2[list(i)].T:
                #                 x,y,z = _mean_confidence_interval(ar, confidence=0.95)
                #                 m2.append(x)
                #                 lb2.append(y)
                #                 ub2.append(z)
                #             m2 = np.array(m2)
                #             ub2 = np.array(ub2)
                #             lb2 = np.array(lb2)

                #             wilc.append([t1[list(i)], t2[list(i)]]) #for calculating wilcoxonstats


                #             #print([M1[i].flatten(), M2[i].flatten()])
                            
                #             if plot != "none":
                #                 if type(plot)!="list":
                #                     plot = [plot]
                #                 if ii in plot:
                #                     fig = plt.figure(1, figsize=(12, 10))
                #                     _plot_mean_and_CI(m1, ub1, lb1, tps, clus=ii, color_mean='red', color_shading='lightcoral', xlab=xlab, ylab=ylab)
                #                     _plot_mean_and_CI(m2, ub2, lb2, tps, clus = ii, color_mean='green', color_shading='lightgreen', xlab=xlab, ylab=ylab)
                #                     #plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_green_{}_conf.pdf".format(j)) ## png
                #                     #plt.savefig(other_path + "/" + str(j) +"_"+ args.name + "_red_green_{}_conf.png".format(j)) ## png
                #                     plt.show()
                #                 else:
                #                     continue


    
    clusterObj.DataSet.add_for_graph(all_g_names, wilc)