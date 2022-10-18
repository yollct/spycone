import pandas as pd
import numpy as np
import sys


def discretization_with_steps(min, max, neg, pos) -> list:
    '''
    min max of the time series object
    neg = negative steps 
    pos = positive steps
    
    '''
    totalsteps = 1
    if min <0:    
        totalsteps += neg
    if max >0:
        totalsteps += pos
    
    discretizevalues = np.zeros(totalsteps)
    index = 0
    if min <0:
        for i in range(neg):
            discretizevalues[index] = (min + (abs(min)/neg)*i)
            index += 1
    discretizevalues[index] = 0
    index += 1
    if max >0:
        for i in range(pos):
            discretizevalues[index] = (max / pos) * i+1
            index+=1
    
    return np.array(discretizevalues)

def discretize_replicates(tsrep, discretizevalues):
    '''
    Handle the replicates
    '''
    newpattern = np.zeros((tsrep.shape))
    newpattern = [discretize_timeseries(i, discretizevalues) for i in tsrep]

    return newpattern


def discretize_timeseries(tsdata, discretizevalues):
    '''
    Discretize the timeseries data, 
    '''
    pattern = np.array(tsdata, dtype="float")
    newvalue = 0
    for i in range(tsdata.shape[0]):
        mindistance = sys.float_info.max
        for j in range(discretizevalues.shape[0]):
            distance = abs(pattern[i]-discretizevalues[j])
            if mindistance > distance:
                mindistance = distance
                newvalue = discretizevalues[j]
        
        pattern[i] = newvalue
    
    return pattern

