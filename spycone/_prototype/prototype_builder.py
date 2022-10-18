import numpy as np

def _prototype_builder(timeseriesobj, clusters, prototype_func):
    '''
    ##calculate prototype of clusters
    #medoid
    ##median or average

    return dict:
    {'34': array([-0.14823925, -0.74057195,  0.41986664, -0.24045551]),
    '146': array([-0.20850356,  0.29922031,  0.74922338,  0.31151074]),
    '195': array([ 0.26338935,  0.19041125,  0.22270178, -0.19392188]),
    '27': array([-0.27389526, -0.6797392 ,  0.98956317,  0.29812246]),
    '78': array([ 0.19380661, -0.91636647,  0.44639439, -0.18150735]),
    '168': array([-1.6519947 , -0.25548007, -0.3339432 , -0.28108202])}
    '''
    cluster_pro = {}
    len_ts = timeseriesobj.shape[0]
    mean_items = 0

    c=0
    for m, n in clusters.items():
        timeseriesobjdf_list = []
        if len(n)>=1:
            timeseriesobj_list = np.asarray(timeseriesobj)[:,list(n),:]
            
            if prototype_func == "mean":
                try:
                    mean_items = np.mean(np.mean(timeseriesobj_list, axis=0), axis=0)
                except:
                    mean_items = np.mean(timeseriesobj_list, axis=0)
            if prototype_func == "median":
                try:
                    mean_items = np.median(np.median(timeseriesobj_list, axis=0), axis=0)
                except:
                    mean_items = np.median(timeseriesobj_list, axis=0)
        
            cluster_pro.setdefault(m, np.array(mean_items, dtype="double"))
        
        
        else:
            continue
        c+=1
    
    for i in clusters.keys():
        if i not in cluster_pro.keys():
            print(i)

    return cluster_pro