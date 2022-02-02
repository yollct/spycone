# cython: infer_types=True

cimport cython
import numpy as np
cimport numpy as np
from cython cimport floating
from libcpp cimport bool


DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


@cython.boundscheck(False)  
cdef tuple _tsis_scores(floating idx, np.ndarray[DTYPE_t, ndim=3] array1, np.ndarray[DTYPE_t, ndim=3] array2):
    """implement tsis 
    return for each pair of isoforms that switch d
    """
    #look for switch points 
    #take before and after switch points

    
    pass

@cython.boundscheck(False)  
cpdef tuple _isoform_differential_usage( np.ndarray[DTYPE_t, ndim=3] normdf, int maj):
    
    cdef Py_ssize_t n_iso = normdf.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=1] iso_ratio = np.zeros(n_iso, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] iso_diff_value = np.zeros(n_iso, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] p_val = np.zeros(n_iso, dtype=DTYPE)

    cdef int maj_higher, min_higher ,x
    cdef DTYPE_t diff_value
    cdef DTYPE_t pp

    for x in range(n_iso):
        if maj == x:
            iso_ratio[x] = 0
            iso_diff_value[x] = 0
            continue
        
        else:
            maj_higher = np.sum(normdf[:,maj,:] > normdf[:,x,:]) 
            min_higher = np.sum(normdf[:,x,:] > normdf[:,maj,:]) 
            diff_value = np.amax(normdf[:,x,:] - normdf[:,maj,:])

            pp = 0
            #calculate ratio
            if min_higher != 0:
                if maj_higher != 0:
                    
                    iso_ratio[x] = min_higher/maj_higher

                else:
                    iso_ratio[x]=0
            else:
                iso_ratio[x] = 0


            iso_diff_value[x] = diff_value
            p_val[x] = pp
    
    return iso_ratio, p_val, maj, iso_diff_value

@cython.boundscheck(False)  
cpdef np.ndarray _isoform_usage(np.ndarray[DTYPE_t, ndim=3] eachgene, bint norm):
    #eachgene 3D array [replicates, isoforms, timepoints]
    cdef:
        np.ndarray[DTYPE_t, ndim=3] result = np.zeros((eachgene.shape[0], eachgene.shape[1], eachgene.shape[2]-1), dtype=DTYPE)
        np.ndarray[DTYPE_t, ndim=2] sizefactor 
        np.ndarray[DTYPE_t, ndim=3] normdf
        Py_ssize_t ngene = eachgene.shape[2]
        np.ndarray[DTYPE_t, ndim=1] finalresult

    if norm:
        #calculate the size factors 
        sizefactor = np.sum(eachgene, axis=1)
        sizefactor[sizefactor==0] = 1

        #normalize to transcript abundance for each replicates
        normdf = np.array([eachgene[x,:,:] / sizefactor[x] for x in range(eachgene.shape[0])])
    
    else:
        normdf = eachgene
    
    #subtract usage from next time point
    for i in range(ngene-1):
        result[:,:,i
        ] = abs(normdf[:,:,i+1]-normdf[:,:,i])
 
    #sum up and take median for one gene across all replicates
    finalresult = np.median(np.sum(result, axis=1),axis=0)
    
    
    return finalresult