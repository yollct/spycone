import numpy as np

def remove_object_notin_network(DataSet, BioNetwork):
    index_to_remove = []
    genelist = DataSet.gene_id
    tmpb = BioNetwork.lst_g()
    # for i in genelist:
    #     if i not in tmpb:
    #         index_to_remove.append([idd for idd, ob in enumerate(genelist) if ob == i][0])

    index_to_remove = np.where(np.isin(genelist, tmpb, invert=True))[0]

    
    DataSet._remove_objects(index_to_remove)
    
    return len(index_to_remove)


def remove_nodes_notin_dataset(DataSet, BioNetwork):
    nodes_to_remove = []
    tmpg = DataSet.gene_id
    nodelist= BioNetwork.lst_g()
    # for i in nodelist:
    #     if i not in tmpg:
    #         nodes_to_remove.append(i)

    nodes_to_remove = nodelist[np.array(np.isin(BioNetwork.lst_g(), DataSet.gene_id, invert=True))]

    BioNetwork._removing_nodes(list(nodes_to_remove))

    return len(nodes_to_remove)

def remove_low_variance(DataSet):
    X = DataSet.ts[0]
    

    rowsums = list(map(int, np.var(X, axis=1)))
    filtered_index = []
    for x,y in enumerate(rowsums):
        if y == 0:
            filtered_index.append(int(x))
    
    DataSet._remove_objects(filtered_index)

    return(len(filtered_index))

def filter_with_cutoff(DataSet, cutoff):
    X = DataSet.ts[0]

    rowsums = np.mean(X, axis=1)
    if cutoff==0:
        filtered_index = np.where(rowsums==cutoff)[0]
    else:
        filtered_index = np.where(rowsums<cutoff)[0]
    
    filtered_genes = [DataSet.gene_id[x] for x in filtered_index]
    DataSet._remove_objects(filtered_index)
    return(len(filtered_index), filtered_genes)


def preprocess(DataSet, BioNetwork=None, remove_low_var=False, cutoff=0):

    """Preprocess data, remove objects without expression along all timepoints.

    Parameters
    ----------

    DataSet : Dataset object.

    BioNetwork : BioNetwork object. If provided, objects that is not in the network will be removed.

    remove_low_var : (boolean) default=False. If true, objects with variance 0 will be removed.

    cutoff : default=0. If given, objects will mean expression across all timepoints lower than the cutoff will be removed.

    Returns 
    ------
    None, changes made directly in the DataSet object


    """
    print(f"Input data dimension: {DataSet.timeserieslist.shape}")
    
    x, filtered_genes = filter_with_cutoff(DataSet, cutoff)
    if cutoff >0:
        print(f"Removed {x} objects lower than {cutoff}")
    else:
        print(f"Removed {x} with 0 values.")
    
    if remove_low_var:
        y= remove_low_variance(DataSet)
        print("Removed {} objects with 0 variance.".format(y))
    
    if BioNetwork is not None:
        i = remove_nodes_notin_dataset(DataSet, BioNetwork)
        print("Removed {} objects from dataset that are not in the network".format(i))
        j = remove_object_notin_network(DataSet, BioNetwork)
        print("Removed {} nodes from network that are not in the dataset (included the genes with lower than cutoff expression).".format(j))

    print("Filtered data: {}".format(DataSet.timeserieslist.shape))


