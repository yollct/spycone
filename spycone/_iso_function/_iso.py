import numpy as np
import pandas as pd
from functools import reduce
from itertools import combinations

def _iso_switch_between_arrays(arr1, arr2, orgarr1, orgarr2):
        ab = np.sign(arr1-arr2)

        ##look for switch points
        points=[[x for x in range(ab.shape[1]) if ab[r,x] != np.roll(ab[r], 1)[x] and x!=0] for r in range(ab.shape[0])]

        #points: switch points of all replicates
        ##check if all r have switch points
        #then take the points where all appeared
        def multiple_intersect(points):
            return np.unique(list(reduce(lambda x, y: set(x).intersection(set(y)), points)))

        #calculate differences
        #[(abs(a[x-1] - a[x]) + abs(b[x-1] -b[x]))/2 for x in points], a,b
        if any(points):
            final_sp = multiple_intersect(points) 
            if not any(final_sp):
                final_sp = [list(set(x[0]).intersection(set(x[1]))) for x in combinations(points,2)]
                final_sp = np.unique(reduce(lambda x,y: x+y, final_sp))
                if not any(final_sp):
                    final_sp = reduce(lambda x,y : x+y, points)
        

            final_sp = np.sort(final_sp)
            # mean of max difference at the time points between 2 arrays for all replicates
            allsp_diff = [np.mean([np.max([abs(arr1[r,x-1] - arr2[r,x-1]), abs(arr2[r,x] -arr1[r,x])]) for r in range(arr1.shape[0])]) for x in final_sp]
            
        else:
            final_sp = []
            allsp_diff=[]

        #calculate corr
        cors = [(1-np.corrcoef(orgarr1[r], orgarr2[r])[0,1])/2 for r in range(arr1.shape[0])]
        
        #for each switch point, the max switch diff before and after 
        #allsp : all the before and after switch differences for all switch points
        return len(final_sp)/arr1.shape[1], allsp_diff, final_sp, np.mean(cors)