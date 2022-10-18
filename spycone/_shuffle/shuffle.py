import numpy as np
import random

class shuffling():
    def __init__(self, timeserieslist=None, shuffle_timeseries=None):
        self.timeserieslist = timeserieslist
        # self.shuffletype = shuffletype
        self.shuffle_timeseries = shuffle_timeseries
  
    def shuffle_dataset_rowwise(self):
        '''
        shuffle time series within one time series 
        '''
        alltimeseriesobj = self.timeserieslist

        for X in alltimeseriesobj:
            for Y in X:
                np.random.shuffle(Y)

        return alltimeseriesobj

    def shuffle_globally(self):
        '''
        shuffle time series within the dataset
        '''
  
        alltimeseriesobj = self.timeserieslist
        objectsize = alltimeseriesobj.shape[1]
        timepoints = alltimeseriesobj.shape[2]
        replicates = alltimeseriesobj.shape[0]

        eachrow = np.zeros((objectsize, timepoints), dtype="double")
        shuffled_rows = np.zeros((objectsize, timepoints), dtype="double")
        shuffled = np.zeros((replicates, objectsize, timepoints), dtype="double")
        randomids = np.random.randint(objectsize*timepoints, size=objectsize*timepoints)

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
