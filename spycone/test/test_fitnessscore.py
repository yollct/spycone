import pandas as pd
import numpy as np
import time
import random
import warnings

from ..inputdata.DataSet import DataSet
from ..inputdata.BioNetwork import BioNetwork
from ..scaler import tanhnormalizer_scaler
from ..feature.get_features import featureObj
from ..fitness.fitness_score import fitness_score

feature_test = featureObj(feature_store=[[]])

fitness_test = fitness_score(feature_test)
value = fitness_test.calculate_fitness()
print(value)