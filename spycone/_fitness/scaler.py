import pandas as pd
import numpy as np
import time
import random
import warnings

from ..DataSet import DataSet
from ..BioNetwork import BioNetwork

def tanhnormalizer_scaler(feature, seed=1222):
    scaled_feature = []
    smean = np.mean(feature)
    svar = np.var(feature)
    
    for i in feature:
        scaled_feature.append(0.5 * (np.tanh(0.01 * i-smean)/svar)+1)

    return scaled_feature

