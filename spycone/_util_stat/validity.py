from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score


def jaccard(c,cnorm):
    intersect = len(set(c).intersection(set(cnorm)))
    if intersect == 0:
        return 0
    total = float(intersect) / len(set(c).union(set(cnorm)))
    return total

def AMIscore(clusterObj, groundtruth):
    '''
    groundtruth: dataframe with two columns names: cluster and object, which is the correct classification of the objects.
    '''
    #convert groundtruth to dict
    '''
    Compute the similarity measure between two clusterings
    '''
    clusters = clusterObj.symbs_clusters
    testdict = {}
    keys = groundtruth['cluster'].unique()
    for key in keys:
        testdict[key] = list(groundtruth[groundtruth['cluster']==key].object)

    true = []
    for c,o in testdict.items():
        for i in o:
            true.append(c)

    pred = []
    for obj in list(groundtruth.object):
        tmp = [pred.append(k) for k, v in clusters.items() for i in v if i == obj]

    return adjusted_mutual_info_score(true, pred)


def adjusted_rand_index(clusterObj, groundtruth):
    '''
    Compute the similarity measure between two clusterings
    '''
    clusters = clusterObj.symbs_clusters
    testdict = {}
    keys = groundtruth['cluster'].unique()
    for key in keys:
        testdict[key] = list(groundtruth[groundtruth['cluster']==key].object)

    true = []
    for c,o in testdict.items():
        for i in o:
            true.append(c)

    pred = []
    for obj in list(groundtruth.object):
        tmp = [pred.append(k) for k, v in clusters.items() for i in v if i == obj]

    return adjusted_rand_score(true, pred)