import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def _cal_pqvals(DataSet):
    clu_val = DataSet.use_clu[0]
    graph_nms = DataSet.object_index[0]
    clst_values = DataSet.object_simval[0]
    mtrx = DataSet.mtrx[0]
    l=DataSet.l[0]
    wilc = DataSet.wilc[0]

    pvals = {}
    mtrx_av = np.average(mtrx)
    assert len(graph_nms) == len(clst_values)
    for i in range(len(graph_nms)):
        if len(graph_nms[i]) >= l:
            pvals[clu_val.keys()[i]] = (stats.mannwhitneyu(clst_values[i], mtrx_av.flatten(), alternative="two-sided")[1] )
        


    nr_of_pvs = len(pvals)
#     for i in wilc:
#             pvals.append(stats.wilcoxon(np.average(i[0],axis=0), np.average(i[1],axis=0) )[1])

    qvals_true = multipletests(pvals.values(), method="fdr_bh")[1]

    qvals = qvals_true[:nr_of_pvs]
    p_gr = pvals.items()[:nr_of_pvs] 
    p_gr = dict(p_gr)
    DataSet.add_pval(p_gr)
    DataSet.add_qval(qvals)