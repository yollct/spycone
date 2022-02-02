import pandas as pd
import numpy as np
import time
import random
import warnings
import itertools as it
cimport cython
cimport numpy as np

from ..DataSet import DataSet
from ..BioNetwork import BioNetwork




cdef class shuffling():
    cdef double[:,:,:] timeserieslist
    # cdef str shuffetype 
    cdef double[:,:,:] shuffle_timeseries
    def __cinit__(self, timeserieslist=None, shuffle_timeseries=None):
        self.timeserieslist = timeserieslist
        # self.shuffletype = shuffletype
        self.shuffle_timeseries = shuffle_timeseries
  


    # cdef double[:,:,:] shufflearray(self, double[:,:,:] shuffle):
    #     cdef list choices = np.arange(shuffle.shape[1])
    #     cdef int index
    #     cdef double[:,:] tmp
        
    #     for n in range(shuffle.shape[1]):
    #         index = random.choice(choices)
    #         choices = np.delete(choices, np.argwhere(choices==index))
    #         tmp = shuffle[:,index,:]
    #         shuffle[:,index,:] = shuffle[:,n,:]
    #         shuffle[:,n,:] = tmp
        
    #     return shuffle      

    cdef double[:,::1] _permute_row(self, double[:,:] row_to_shuffle):
        cdef int timepoints = row_to_shuffle.shape[0]
        cdef int replicates= row_to_shuffle.shape[1]
        cdef double[:,::1] newpattern = np.empty((replicates, 1, timepoints))

        for r in row_to_shuffle:
            np.random.shuffle(r)

        newpattern = row_to_shuffle.copy()

        return newpattern
                      
    # cpdef double[:,:,:] shuffle_timeseries(self):
    #     cdef double[:,:,:] alltimeseriesobj = self.clusterobj.DataSet.timeserieslist
    #     cdef double[:,:,:] shuffle, shuffled

    #     shuffle = alltimeseriesobj.copy()

    #     shuffled = self._shufflearray(shuffle)

    #     self.shuffle_timeseries = shuffled

    #     return shuffled 

    cpdef double[:,:,:] shuffle_dataset_rowwise(self):
        '''
        shuffle time series within one time series 
        '''
        
        cdef double[:,:,:] alltimeseriesobj = self.timeserieslist
        cdef double[:,:] eachrow
        cdef double[:] shuffled_rows 

        for X in alltimeseriesobj:
            for Y in X:
                np.random.shuffle(Y)

        return alltimeseriesobj


    cpdef double[:,:,:] shuffle_globally(self):
        '''
        shuffle time series within the dataset
        '''
  
        cdef double[:,:,:] alltimeseriesobj = self.timeserieslist
        cdef Py_ssize_t objectsize = alltimeseriesobj.shape[1]
        cdef Py_ssize_t timepoints = alltimeseriesobj.shape[2]
        cdef Py_ssize_t replicates = alltimeseriesobj.shape[0]
        cdef int randomid, randomobj, randomtp

        random.seed()
        
        cdef double[:,:] eachrow = np.zeros((objectsize, timepoints), dtype="double")
        cdef double[:,:] shuffled_rows = np.zeros((objectsize, timepoints), dtype="double")
        cdef double[:,:,:] shuffled = np.zeros((replicates, objectsize, timepoints), dtype="double")
        cdef long[:] randomids = np.random.randint(objectsize*timepoints, size=objectsize*timepoints)

        #cython

        
        for r in range(replicates):
            eachrow = np.zeros((objectsize, timepoints), dtype="double")
            for o in range(objectsize):
                for tp in range(timepoints):
                    randomid = random.choice(randomids)
                    randomobj = randomid // timepoints
                    randomtp = randomid % timepoints
                    eachrow[o][tp] = alltimeseriesobj[r,randomobj,randomtp]

                shuffled_rows[o] = eachrow[o]


            shuffled[r] = shuffled_rows
     
        #multiprocessing here
        # eachrow = [alltimeseriesobj[:,i,:] for i in range(alltimeseriesobj.shape[1])]
       
        # shuffled_rows = []
        # with ProcessPoolExecutor() as exe:
        #     for row in exe.map(self._permute_row, eachrow):
        #         shuffled_rows.append(row)
        
        # #swap axes of the np.array for 0 to 1
        # shuffled_rows = np.swapaxes(shuffled_rows, 0, 1)

        # shuffled = self._shufflearray(shuffled_rows)

        self.shuffle_timeseries = shuffled 

        return shuffled

    # cpdef create_random_dataset(self):
    #     alltimeseriesobj = self.clusterobj.DataSet.timeserieslist

    #     minval = np.min(alltimeseriesobj)
    #     maxval = np.max(alltimeseriesobj)

    #     randomdataset = np.random.uniform(minval, maxval, size=(alltimeseriesobj.shape))

    #     self.shuffle_timeseries = randomdataset

    #     return randomdataset

  
        

        



        

