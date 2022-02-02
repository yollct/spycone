
import pandas as pd
import time
import numpy as np
from cython.parallel import prange

import itertools as it
cimport cython
cimport numpy as np
from libc.math cimport pow, sqrt
from libc.stdlib cimport malloc

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
cpdef double anti_pearson_corr(double[:] x, double[:] y, int size):
    cdef double r
    r = pearson_corr(x,y,size)

    if r == 0:
        return 0

    return (2-r)/2

@cython.boundscheck(False)
cpdef double[:,:] anti_corr_mat(double[:,:] arr1, double[:,:]arr2):
    cdef double[:,:] result_mat = np.zeros((len(arr1), len(arr2)), dtype="double")
    cdef int size=len(arr1[0])
    assert len(arr1[0]) == len(arr2[0])

    for a in range(len(arr1)):
        for b in range(len(arr2)):
            pr = anti_pearson_corr(arr1[a],arr2[b], size)
            result_mat[a][b] = pr
    
    return result_mat

@cython.boundscheck(False)
cpdef double euclidean_dist(double[:] x, double[:] y, int size):
    cdef double distance = 0
    cdef double diff 
    for i in range(size):
        diff = 0
        diff = x[i] - y[i]
        distance += diff * diff

    distance = sqrt(distance)
    return -distance

@cython.boundscheck(False)
cpdef double[:,:] dist_mat(double[:,:] arr1, double[:,:]arr2):
    cdef double[:,:] result_mat = np.zeros((len(arr1), len(arr2)), dtype="double")
    cdef int size=len(arr1[0])
    assert len(arr1[0]) == len(arr2[0])

    for a in range(len(arr1)):
        for b in range(len(arr2)):
            pr = euclidean_dist(arr1[a],arr2[b], size)
            result_mat[a][b] = pr
    
    return result_mat

@cython.boundscheck(False)
cdef double var_for_pearson(double[:] a, double mean, int size):
    cdef double[:] squares = np.zeros(size)
    cdef double var=0
    cdef double meansq=0
    cdef double sqmean=0

    for i in range(size):
        meansq += pow(a[i], 2)
    meansq /= size
    sqmean = pow(mean,2)
    var = meansq - sqmean

    return var

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
cpdef double pearson_corr(double[:] x, double[:] y, int size):
    cdef double rho
    cdef double xmean=0
    cdef double ymean=0
    cdef double cov=0
    cdef double xsd=0
    cdef double ysd=0
    cdef double[:] xnorm=np.zeros(size)
    cdef double[:] ynorm=np.zeros(size)
    
    cov=0
    #arithmetic mean
    for i in range(size):
        xmean += x[i]
        ymean += y[i] 

    xmean /= size
    ymean /= size

    #covariance
    
    for i in range(size):
        xnorm[i] = x[i] - xmean
        ynorm[i] = y[i] - ymean
        cov += xnorm[i] * ynorm[i]
        
    cov /= size
    #standard deviation

    xsd = var_for_pearson(x, xmean, size)
    ysd = var_for_pearson(y, ymean, size)

   
    rho = cov / sqrt(xsd*ysd)
    if not np.isfinite(rho) or np.isnan(rho):
        return 0
    

    return (1+rho)/2

cpdef double[:,:] corr_mat(np.ndarray arr1, np.ndarray arr2):
    cdef double[:,:] result_mat = np.zeros((len(arr1), len(arr2)), dtype="double")
    cdef int size=len(arr1[0])
    assert len(arr1[0]) == len(arr2[0])

    for a in range(len(arr1)):
        for b in range(len(arr2)):
            pr = pearson_corr(arr1[a],arr2[b], size)
            result_mat[a][b] = pr
    
    return result_mat

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
cpdef tuple fast_calculate_similarities_timeseriesobj_to_medoid(double[:,:,:] alltimeseries, np.ndarray[np.int64_t, ndim=1] timeseriesobj_id, dict medoids, str similarityfunction):
    '''
    Calculate similarities between each timeseries object and medoids. 
    Replicates are taken care of.
        closestmedoid = []
    '''
    cdef int n_allobjects, closestid
    n_allobjects = timeseriesobj_id.shape[0]

    cdef double[:,:] tsdata1 
    cdef double optim_sim 
    cdef dict clusters = {}
    cdef double totalsim =0 

    cdef double[:] second_med_dist = np.zeros(alltimeseries.shape[1])
    cdef double[:] first_med_dist = np.zeros(alltimeseries.shape[1])

    #use parallel
    for idx in range(n_allobjects):
        # print("objid", timeseriesobj_id[objidx])
        objidx = timeseriesobj_id[idx]
        optim_sim = -np.inf
        second_optim_sim = -np.inf
        for medoid_index, medoidarray in medoids.items():
            tsdata1 = alltimeseries[:,objidx,:]
            #closestid, sim = calculate_closest_medoid_for_each_obj(objidx, n_timeserieslist, alltimeseries, obj_to_compare, n_objtocomp, similarityfunction)
            tmp_sim = calculate_timeseriesobj_pairs_similarities(tsdata1, medoidarray, similarityfunction)
            
            if tmp_sim > optim_sim:
                optim_sim = tmp_sim
                optim_medoid_index = medoid_index
            elif tmp_sim > second_optim_sim:
                second_optim_sim = tmp_sim

        if str(optim_medoid_index) in clusters.keys():
            clusters[optim_medoid_index].append(objidx)
        else: 
            clusters.setdefault(str(optim_medoid_index), [objidx])
        
        totalsim += optim_sim
        second_med_dist[idx] = second_optim_sim
        first_med_dist[idx] = optim_sim
        
    return clusters, second_med_dist, first_med_dist, totalsim

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
##for each clusters, update medoids with the best similarity within the cluster
cpdef tuple update_medoids(double[:,:,:] alltimeseriesobj, np.ndarray[np.int64_t, ndim=1]   alltimeseries_id, int[:] cluster_obj, str cluster_medoid, str similarityfunction):
    
    #cdef double[:,:] cluster_obj_timeseriesarray = alltimeseriesobj[:,cluster_obj,:]
    cdef double distance_ts1, best_dist
    cdef int medoid_index
    best_dist = -np.inf

    for ts1 in cluster_obj:
        distance_ts1 = 0
        for ts2 in cluster_obj:
            if ts1 == ts2:
                continue
            distance_ts1 += calculate_timeseriesobj_pairs_similarities(alltimeseriesobj[:,ts1,:], alltimeseriesobj[:,ts2,:], similarityfunction)

        #max dist 
        if distance_ts1 > best_dist:
            best_dist = distance_ts1
            medoid_index = ts1
        
    return medoid_index, best_dist

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef double calculate_medoids_difference(double[:,:,:] current_medoid, double[:,:,:] updated_medoid, str similarityfunction):
    cdef double tmp_dis =0
    cdef double result

    for cm_index in range(current_medoid.shape[0]):
        for um_index in range(updated_medoid.shape[0]):
            tmp_dis += calculate_timeseriesobj_pairs_similarities(current_medoid[cm_index], updated_medoid[um_index], similarityfunction)
    
    return 1-tmp_dis/current_medoid.shape[0]/updated_medoid.shape[0]

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef double[:,:] clusters_pair_similarity(dict prototype, str similarityfunction):
    '''
    Calculate inter cluster similarity using clusters prototype
    here takes the maximum of clusters similarity because the pearson is higher the more similar
    '''

    cdef list clist = list(prototype.keys())
    cdef double cmax = -np.inf
    cdef double tmp_max 
    cdef double[:,:] clusters_pair_sim = np.zeros((len(prototype),len(prototype)-1))

    #loop for each combination of clusters pairs
    n=0
    for o in range(len(clist)):
        k = clist[o]
        m=0
        for l in clist[0:o]+clist[o+1:]:                
            tmp_max = 0
            if similarityfunction=="correlation":
                tmp_max = pearson_corr(np.array(prototype[k]).astype("double"), np.array(prototype[l]).astype("double"), len(prototype[k]))

            if similarityfunction=="euclidean":
                tmp_max = euclidean_dist(np.array(prototype[k]).astype("double"), np.array(prototype[l]).astype("double"), len(prototype[k]))
            clusters_pair_sim[n][m] = tmp_max 
            m+=1
        n+=1
    return clusters_pair_sim


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef double[:] intra_clusters_similarity(double[:,:,:] alltimeseriesobj, dict clusters, dict prototype, str similarityfunction):
    '''
    calculate similarity from all points in a cluster to prototype
    max distance == min correlation
    '''
    cdef double dist_med
    cdef double[:] avg_eachclusters = np.zeros(len(clusters))

    n=0
    for c, o in clusters.items():
        dist_med = 0
        for i in o:
            # if prop is False:
            #     dist_med += calculate_timeseriesobj_pairs_similarities(alltimeseriesobj[:,i,:], prototype[c], similarityfunction)

            dist_med += calculate_timeseriesobj_pairs_similarities(alltimeseriesobj[:,i,:], np.expand_dims(prototype[c], axis=0).astype("double"), similarityfunction)

        avg_eachclusters[n] = dist_med / len(o)
        n+=1
    
    return avg_eachclusters


@cython.boundscheck(False)
cpdef double calculate_timeseriesobj_pairs_similarities(double[:,:] timeseries1, double[:,:] timeseries2, simfunc):
    '''
    actual euclidean distance between two time series objects.
    replicates are tco.
    '''
    cdef double simsum =0
    cdef double count 

    for ts1 in range(timeseries1.shape[0]):
        for ts2 in range(timeseries2.shape[0]):
            if simfunc == "pearson":
                simsum+= pearson_corr(timeseries1[ts1], timeseries2[ts2], timeseries1.shape[1])

            if simfunc == "euclidean":
                simsum+= euclidean_dist(timeseries1[ts1], timeseries2[ts2], timeseries1.shape[1])
         

    return simsum/timeseries1.shape[0]/timeseries2.shape[0]




##hierarchical clustering (scellnetor)
def _coef_of_determination(y, yhat):

        """
        y is the matrix with observed values
        yhat is vector with modelled values (or the average)

        R2 is a statistic that will give some information about the goodness of
        fit of a model. In regression, the R2 coefficient of determination is a
        statistical measure of how well the regression predictions approximate
        the real data points. An R2 of 1 indicates that the regression predictions
        perfectly fit the data.

        In all instances where R2 is used, the predictors are calculated by
        ordinary least-squares regression: that is, by minimizing SSres (residuals).

        Suppose R2 = 0.49. This implies that 49% of the variability of the
        dependent variable has been accounted for, and the remaining 51% of
        the variability is still unaccounted for.
        """

        y = np.array(y)
        yhat = np.array(yhat)
        mean_of_obs = np.mean(y)  # , axis=0)
        tot_sum_of_sqrs = np.sum((y - mean_of_obs) ** 2)
        # reg_sum_of_sqrs = np.sum((yhat-mean_of_obs)**2)
        residuals = 0
        for i in y:
            residuals += np.sum((i - yhat) ** 2)
        R = 1 - (residuals / tot_sum_of_sqrs)
        return R

cdef tuple _sorting_order_of_genes(list list_entrz):
        cdef list true
        cdef list temp =[]
        cdef list rearrange_list = []
        cdef double no, no1
  
        
        for i in list_entrz[0]:
            if list_entrz[0].count(i) > 1:
                print(i, list_entrz[0].count(i))
        true = list(set(list_entrz[0])) ## remove potential replications...
        true.sort()

        ##removing genes from true that are not in all of the other lists
        no = 0
        temp1 = [set(i) for i in list_entrz[1:]]

        for i in true:
            for ii in temp1:
                if i not in ii:
                    no += 1
                    temp.append(i)
                    break

        for i in temp:
            true.remove(i)

        ## Sorting genes and removing genes from the other replicates
        ## that are not in the "true" list
        no1 = -1
        rearrange_list = []
        for j, i in enumerate(list_entrz):
            if true != i:
                no1 += 1
                temp = []
                for ii in true:
                    temp.append(i.index(ii))
                rearrange_list.append(temp)
                print("\nThe order of the genes in", no1, "replicate(s) was changes")
                print("Removed", len(i) - len(true), "genes in total \n")
            else:
                rearrange_list.append([])

        return rearrange_list, true

def _merge_samples(g1, reps1, mid, gene_list1, tps, symbs1=None):
        #G1 = pd.read_table(g1[0], index_col=0)
        tmstps = tps
        G1 = pd.DataFrame(g1.reshape(-1,tmstps))
        gl1 = pd.Series([str(i) for i in list(gene_list1)], name="Entrez")

        if symbs1 is not None:
            sym1 = pd.Series([str(i) for i in list(symbs1)],name="Symbol")
        else: 
            sym1 = gl1

        G1 = pd.concat([gl1, sym1, G1], axis=1)  
        G1 = G1.sort_values(by=['Entrez'])

        #print(G1['Entrez'])
        #G1 = pd.read_table(g1[0]) ## This is for the control data
        #G1 = pd.read_table(g1[0], index_col=0, sep=",")
        lng1 = int(len(G1)/reps1)
        #lng1 = int(len(G1)/reps1)
        lst1 = []
        #entrez1 = G1['EntrezID'][:lng1]  # g1.index
        entrz = []
        symbs = []

        for i in range(reps1):
            temp = G1[i * lng1: (i + 1) * lng1]
            #lst1.append(temp.iloc[:, -tmstps:].values)
            lst1.append(temp.iloc[:, 2:].values)
            #print(G1.columns, lng1, len(temp), len(lst1[0][0]))
            ##entrz.append([str(gn) for gn in G1['EntrezID'][i * lng1: (i + 1) * lng1]])
            #entrz.append([str(gn) for gn in G1.iloc[:,1][i * lng1: (i + 1) * lng1]])
            

        entrez = G1.iloc[0:lng1,0].values
        symbs = G1.iloc[0:lng1, 1].values
        nplst1 = np.array(lst1, dtype='float32').reshape((reps1, lng1, tmstps))
        #lst1 = np.array(lst1, dtype='float32')
        """
        print("sorting the order of genes and removing genes that are not present in all replicates")
        # entrez is a list with lists of all the different entrez ids in the samples
        sorted_lst, entrez1 = sorting_order_of_genes(entrz)

        for srt in range(len(sorted_lst)):
            if len(sorted_lst[srt]) > 0:
                lst1[srt] = lst1[srt][sorted_lst[srt]]

        ###################
        #entrz[0] = entrez1
        for i in range(len(sorted_lst)):
            if len(sorted_lst[i]) > 0:
                entrz[i] = [entrz[i][x] for x in sorted_lst[i]]

            assert entrz[0] == entrz[i]

        ###################
        """
        ### Calculate the average/median for the data matrix
        if mid == 'average':
            GG1 = np.average(nplst1,axis=0)

        if mid == 'median':
            GG1 = np.median(nplst1, axis=0)


        ### Calcualte the coefficient of determinant
        if reps1 == 1:
            r1 = np.ones((lng1,))
        if reps1 > 1:
            r1 = []
            for i in range(lng1):
                r1.append( _coef_of_determination(nplst1[:,i], GG1[i]) )

        return GG1, entrez, r1, symbs  #entrez1
        
def _all_set_same_length(t1, r1, gene_list, lst_g, adj, n_group = int, t2=None,r2=None, gene_list2=None):
    if t2==None:        
        
        ## Removing genes from gene_list if it contains genes that are not in the PPI network\
        finalgenelist_index = np.where(pd.Series(gene_list).isin(lst_g))[0]
        print(finalgenelist_index.shape)
        finalgenelist = np.array(gene_list)[finalgenelist_index]

        ## Removing genes from the network if it contains genes that are not in either gene lists
        #dlg = np.where(pd.Series(lst_g).isin(temp1)==True)[0]
        
        finalnetwork_index = np.where(pd.Series(lst_g).isin(gene_list[finalgenelist]))[0]
        finalinnetwork = np.array(lst_g)[finalnetwork_index]
        adj = adj[finalnetwork_index, finalnetwork_index]


        t1 = np.delete(t1, (finalgenelist_index), axis=0)
        r1 = np.delete(np.array(r1), (finalgenelist_index), axis=0)


        print("everything has the same lenghts:")
        print(len(r1), len(finalgenelist), len(adj), len(t1), len(finalinnetwork))
        

        return t1, finalgenelist, r1, adj, finalinnetwork

    

def _make_combined_simmatrix_euclidean(d1,lst_g, adj, r1, reps1, corr_type=True, dist="euclidean"):
    
    import sklearn.metrics
    from scipy.spatial import distance
    #from scipy.spatial.distance import pdist
    from scipy.spatial.distance import cdist
    
    adj_mtrx = adj

    ### Get the coefficients for calculating the similarity function (euclidean functions)
    T = d1.shape[1]

    ### Normalize the matrix
    cdef np.ndarray M = np.array(d1).copy()
    M = M - np.min(M)
    cdef np.ndarray M1 = M/np.max(M)
    ### Calculate the similarity matrix
    slk = sklearn.metrics.pairwise.euclidean_distances(np.array(M1),np.array(M1))
    slk = np.array(slk, dtype='double')

    xmax = np.max(slk)
    xmin = np.min(slk)
    ME = T*(xmax-xmin)
    slk = (ME - slk)
    slk = slk/ME
    
    
    ### Normalize the matrix again
#         slk = slk - np.min(slk)
#         slk1 = slk / np.max(slk)
    slk = abs(slk - 1)
    
    #print('size of simmatrix: ', slk.strides)
    return slk, adj_mtrx, lst_g

def _make_combined_simmatrix_pearson(d1, lst_g, adj, r1, reps1, gene_list, corr_type=True, dist="pearson"):
    import sklearn.metrics
    from scipy.spatial import distance
    #from scipy.spatial.distance import pdist
    from scipy.spatial.distance import cdist
    from scipy.stats import pearsonr
    adj_mtrx = adj
    
    ### Normalize the matrix
    M = np.array(d1).copy()
    M = M - np.min(M)
    M1 = M/np.max(M)
    
    ncovs = []

    for i,j  in enumerate(M1):
        c = np.cov(j)
        ncovs.append(c)
    
    #stores indices with 0 covariance for delete
    noncov = np.where(np.array(ncovs)==0)[0]
    M1_cov = np.delete(M1, (noncov), axis=0)
    
    #delete also in gene list and other things 
    lst_g = np.delete(lst_g, (noncov))
    adj_mtrx = np.delete(adj_mtrx, (noncov), axis=0)
    adj_mtrx = np.delete(adj_mtrx, (noncov), axis=1)
    gene_list = np.delete(gene_list, (noncov))
    r1 = np.delete(r1, (noncov), axis=0)
    ### Calculate the similarity matrix
    #generate pairwise pearson correlation matrix
    corre = np.corrcoef(M1_cov)
    corre = np.nan_to_num(corre)
    slk = corre
    slk = np.array(slk, dtype=float)
    adj_mtrx = np.array(adj_mtrx, dtype='uint8')
    slk = slk + 1
    slk = slk/2
        
        
        ### Normalize the matrix again
        # slk = slk - np.min(slk)
        # slk1 = slk / np.max(slk)
        # slk1 = abs(slk1 - 1)
    
    print(adj_mtrx.shape, len(lst_g), len(gene_list), len(r1))
    #print('size of simmatrix: ', slk.shape)
    return slk, adj_mtrx, lst_g, gene_list, r1


#kemdoids
cimport cython

import numpy as np
cimport numpy as np
from cython cimport floating, integral

@cython.boundscheck(False)  # Deactivate bounds checking
cpdef _compute_optimal_swap( floating[:,:] D,
                           int[:] medoid_idxs,
                           int[:] not_medoid_idxs,
                           floating[:] Djs,
                           floating[:] Ejs,
                           int n_clusters):
    """Compute best cost change for all the possible swaps."""

    # Initialize best cost change and the associated swap couple.
    cdef (int, int, floating) best_cost_change = (1, 1, 0.0)
    cdef int sample_size = len(D)
    cdef int i, j, h, id_i, id_h, id_j
    cdef floating cost_change
    cdef int not_medoid_shape = sample_size - n_clusters
    cdef bint cluster_i_bool, not_cluster_i_bool, second_best_medoid
    cdef bint not_second_best_medoid

    # Compute the change in cost for each swap.
    for h in range(not_medoid_shape):
        # id of the potential new medoid.
        id_h = not_medoid_idxs[h]
        for i in range(n_clusters):
            # id of the medoid we want to replace.
            id_i = medoid_idxs[i]
            cost_change = 0.0
            # compute for all not-selected points the change in cost
            for j in range(not_medoid_shape):
                id_j = not_medoid_idxs[j]
                cluster_i_bool = D[id_i, id_j] == Djs[id_j]
                not_cluster_i_bool = D[id_i, id_j] != Djs[id_j]
                second_best_medoid = D[id_h, id_j] < Ejs[id_j]
                not_second_best_medoid = D[id_h, id_j] >= Ejs[id_j]

                if cluster_i_bool & second_best_medoid:
                    cost_change +=  D[id_j, id_h] - Djs[id_j]
                elif cluster_i_bool & not_second_best_medoid:
                    cost_change +=  Ejs[id_j] - Djs[id_j]
                elif not_cluster_i_bool & (D[id_j, id_h] < Djs[id_j]):
                    cost_change +=  D[id_j, id_h] - Djs[id_j]

            # same for i
            second_best_medoid = D[id_h, id_i] < Ejs[id_i]
            if  second_best_medoid:
                cost_change +=  D[id_i, id_h]
            else:
                cost_change +=  Ejs[id_i]

            if cost_change < best_cost_change[2]:
                best_cost_change = (id_i, id_h, cost_change)

    # If one of the swap decrease the objective, return that swap.
    if best_cost_change[2] < 0:
        return best_cost_change
    else:
        return None




cpdef _build( floating[:, :] D, int n_clusters):
    """Compute BUILD initialization, a greedy medoid initialization."""

    cdef int[:] medoid_idxs = np.zeros(n_clusters, dtype = np.intc)
    cdef int sample_size = len(D)
    cdef int[:] not_medoid_idxs = np.zeros(sample_size, dtype = np.intc)
    cdef int i, j,  id_i, id_j

    medoid_idxs[0] = np.argmin(np.sum(D,axis=0))
    not_medoid_idxs = np.delete(not_medoid_idxs, medoid_idxs[0])

    cdef int n_medoids_current = 1

    cdef floating[:] Dj = D[medoid_idxs[0]].copy()
    cdef floating cost_change
    cdef (int, int) new_medoid = (medoid_idxs[0], 0)
    cdef floating cost_change_max

    for _ in range(n_clusters -1):
        cost_change_max = 0
        for i in range(sample_size - n_medoids_current):
            id_i = not_medoid_idxs[i]
            cost_change = 0
            for j in range(sample_size - n_medoids_current):
                id_j = not_medoid_idxs[j]
                cost_change +=   max(0, Dj[id_j] - D[id_i, id_j])
            if cost_change >= cost_change_max:
                cost_change_max = cost_change
                new_medoid = (id_i, i)


        medoid_idxs[n_medoids_current] = new_medoid[0]
        n_medoids_current +=  1
        not_medoid_idxs = np.delete(not_medoid_idxs, new_medoid[1])


        for id_j in range(sample_size):
            Dj[id_j] = min(Dj[id_j], D[id_j, new_medoid[0]])
    return np.array(medoid_idxs)