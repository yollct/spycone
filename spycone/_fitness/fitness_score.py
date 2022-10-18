import numpy as np

class fitness_score():
    def __init__(self, featureObj):
        self.featureObj = featureObj
        

    def calculate_fitness(self):
        object_type = self.featureObj.feature_type
        result = 0
        feature_store = self.featureObj.feature_store
        final_fvalue=[]

        for f in range(len(feature_store)):
            feature = feature_store[f]
            final_fvalue.append(feature)

        final_fvalue = np.array(final_fvalue)

        return final_fvalue

