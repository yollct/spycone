from ticone.clustering import similarities
from similarities import pearson_corr
from array import array

def main():
    data1 = array('d', [0,1,2,3,4])
    data2 = array('d', [0,-1,-2,-3,-4])
    
    corr = pearson_corr(data1, data2, 5)  
    print(corr)  
    assert corr == 0

    data2 = array('d', [0, 2, 4, 6, 8])
    corr = pearson_corr(data1, data2, 5)
    
    assert corr == 1

    data2 = array('d', [0, 2, 0, 2, 0])
    corr = pearson_corr(data1, data2, 5)

    assert corr == 0.5