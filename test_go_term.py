import sys
import numpy as np
import os
import pandas as pd

sys.path.insert(0, os.path.abspath('.'))
from spycone.go_terms import list_genesets, clusters_gsea, modules_gsea
from spycone.DataSet import DataSet as dataset
from spycone.clustering import clustering


genelist = ["ENSG00000126353","ENSG00000108839","ENSG00000132965","ENSG00000179148","ENSG00000171824","ENSG00000105329","ENSG00000021762","ENSG00000188676","ENSG00000118322","ENSG00000180113","ENSG00000110074","ENSG00000085563","ENSG00000132326","ENSG00000168389","ENSG00000197580","ENSG00000129484","ENSG00000206503","ENSG00000145888","ENSG00000198691","ENSG00000167972"]
tsmat = pd.DataFrame(np.random.rand(20,21))
species = "hsapiens"

dset = dataset(tsmat,
            species=9606,
            keytype="ensemblgeneid",
            reps1=7,
            timepts=3,
            gene_id=genelist,
            symbs=genelist)

asclu=clustering(dset, input_type="expression", algorithm="hierarchical", linkage="ward", metric="euclidean", n_clusters=3)
clu = asclu.find_clusters()

go_res=clusters_gsea(dset, species="hsapiens", method="gsea", cutoff=1)
print(go_res)