import numpy as np
from sklearn.cluster import *
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import pairwise_distances, silhouette_score
from tslearn.metrics import cdist_soft_dtw
from tslearn.clustering import TimeSeriesKMeans
from tslearn.clustering import silhouette_score as tssilhouette_score
from sklearn.metrics import davies_bouldin_score

