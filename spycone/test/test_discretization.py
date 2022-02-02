import pandas as pd
import numpy as np
import time
from scipy.spatial.distance import cdist
from numba import njit, prange
import multiprocessing as mp
from sklearn.metrics.pairwise import euclidean_distances
import sys

from ticone.preprocessing._discretization import discretization_with_steps, discretize_timeseries, discretize_replicates

inputts = np.array([[-3, -2.5, -0.5, 0, 0.49, 2.51]])

discretizeprototype = discretization_with_steps(-3,3,3,3)
discretizepattern = discretize_replicates(inputts, discretizeprototype)

print(discretizepattern)