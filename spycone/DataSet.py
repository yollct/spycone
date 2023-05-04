import numpy as np
import pandas as pd
from collections import defaultdict
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
from ._preprocessing._discretization import discretization_with_steps, discretize_timeseries, discretize_replicates

def _map_to_gene_name(keys, keytype, species, transid_type = None):
    keymap = {'entrezgeneid': 'NCBI gene (formerly Entrezgene) ID', 'ensemblgeneid': 'Gene stable ID', 'ensembltransid' : 'Transcript stable ID', 'ensembltransidver': "Transcript stable ID version"}
    ann = pd.read_csv(os.path.join(dir_path, "data/annotation/{}_ann.csv".format(species)), index_col=keymap[keytype], dtype=str)
    ann.index = list(map(str, list(map(lambda x: int(x) if not np.isnan(x) else x, ann.index))))
    ann_dict = ann.to_dict()
    
    symbs = []
    for k in keys:
        try:
            symbs.append(ann_dict['Gene name'][str(k)])
        except:
            symbs.append('n')

    return symbs

def _map_trans_to_gene(keys, keytype, species, transid_type = None):
    keymap = {'entrezgeneid': 'NCBI gene (formerly Entrezgene) ID', 'ensemblgeneid': 'Gene stable ID', 'ensembltransid' : 'Transcript stable ID', 'ensembltransidver': "Transcript stable ID version"}
    ann = pd.read_csv(os.path.join(dir_path, "data/annotation/{}_ann.csv".format(species)), index_col=keymap[transid_type], dtype=str)
    ann_dict = ann.to_dict()

    a = []
    for k in keys:
        try:
            a.append(str(ann_dict[keymap[keytype]][k]))
        except:
            a.append('n')

    return a

##dataset and settings
class dataset():
    """ 
    Input dataset.

    Parameters
    ----------

    ts 
        matrix or dataframe of the time series dataset.

    gene_id 
        list of gene id. (id that matched the biological network, default network : entrez ID). If gene_id is not given, entrez gene id will be mapped. 

    transcript_id 
        list of transcript id if the expression matrix is in transcript-level.

    species
        Specifying species ID for annotation (human : 9606, mouse: 10090)

    reps1 
        number of replicates 

    timepts 
        number of time points

    gtf : (optional)
        provide corresponding gtf file for mapping gene names 

    symbs : (optional)
        list of gene symbols or gene names (can be automatically mapped for human and mouse)

    discretization_steps : (optional)
        discretize the expression matrix according to the number of steps given


    Attributes
    -----------
    timeserieslist 
        automatically generated 3D-array for analysis

    Methods
    ----------
    __copy__ 
        copy function for the dataset object

    remove_objects 
        remove the given index of the object

    """
    SPECIES = [9606, 10090]
    KEYTYPE = ['entrezgeneid', 'ensemblgeneid', 'ensembltransid', 'ensembltransidver']
    def __init__(
        self, 
        ts, 
        species,
        #keytype,
        reps1,
        timepts,
        gtf=None,
        gene_id=None,
        transcript_id = None, 
        timeserieslist = None,
        symbs = None,
        discretization_steps = None):
        """
        ts : matrix of the time series data.
        gene_list : list of the gene id
        reps1 : Number of replicates for ts.
        timepts : Number of time points.
        """
        self.ts = [ts]
        self.reps1 = reps1
        self.gene_id = gene_id
        self.species=species
        #self.keytype=keytype
        self.timepts = timepts
        self.timeserieslist = []
        self.symbs = symbs
        self.transcript_id = transcript_id
        self.discretization_steps = discretization_steps

        # if not isinstance(self.reps1, list):
        #     self.reps1 = [self.reps1]
        #check attributes
        if not hasattr(self.ts, "shape"):
            self.ts[0] = np.array(self.ts[0], dtype="double")

        if self.timepts*self.reps1 != self.ts[0].shape[1]:
            raise ValueError("Number of columns is not the same as number of time points.")

        if self.species not in self.SPECIES:
            raise ValueError("Please provide a supported species ID.")

        # if self.keytype not in self.KEYTYPE:
        #     raise ValueError("Please provide a supported key type.")

        ######gene list check
        if self.gene_id is None and self.transcript_id is None:
            raise ValueError("Please provide at least gene ID or transcript ID.")

        TRANS_ID_TYPE = None
        if self.transcript_id is not None:
            if "." in self.transcript_id[0]:
                TRANS_ID_TYPE = "ensembltransidver"
            else:
                TRANS_ID_TYPE = "ensembltransid"

        # if self.transcript_id is not None and self.gene_id is None:
        #     self.gene_id = _map_trans_to_gene(self.transcript_id, self.keytype, self.species, transid_type=TRANS_ID_TYPE)

        # if self.symbs is None:
        #     self.symbs = _map_to_gene_name(self.gene_id, self.keytype, self.species, transid_type=TRANS_ID_TYPE)

        if not isinstance(self.gene_id, list):
            self.gene_id = [str(x) for x in self.gene_id]
            
        if not isinstance(self.symbs, (np.ndarray)):
            self.symbs = np.array(self.symbs)
        # if (self.ts[0].shape[0] != self.gene_list[0] | self.gene_list[0] != self.symbs1[0]).any():
        #     raise ValueError("Input length do not match.") 

        if pd.isna((np.array(self.gene_id))).all():
            raise ValueError("The gene ID list contains NaN values.")

        if pd.isna(self.ts).all():
            raise ValueError("The count matrix contains NaN values.")
        

        tmp = []
        
        if self.reps1 > 1:
            tp = self.timepts
            ts = self.ts[0]
            rep = int(self.reps1)

            for r in range(rep):
                start=tp*r
                end = start+tp
                tmp.append(ts[:,start:end])
            
            self.timeserieslist = np.array(tmp, dtype="double")
        else:
            self.timeserieslist = np.expand_dims(ts, axis=0)

        ##discretization
        if self.discretization_steps is not None:
            data_min = np.min(self.timeserieslist)
            data_max = np.max(self.timeserieslist)
            neg = self.discretization_steps[0]
            pos = self.discretization_steps[1]
            discretizeprototype = discretization_with_steps(data_min, data_max, neg, pos)
            

            for i in range(self.timeserieslist.shape[1]):
                tmp_ts = self.timeserieslist[:,i,:]
                discretizepattern = discretize_replicates(tmp_ts, discretizeprototype)
                self.timeserieslist[:,i,:] = tmp_ts 
        
        if self.transcript_id is not None:
            self._get_gene_level()

    @property
    def shape(self):
        return self.timeserieslist.shape

    def __copy__(self):
        return DataSet(ts=self.ts[0],
        gene_id=self.gene_id,
        reps1=self.reps1,
        timepts=self.timepts,
        species=self.species,
        keytype=self.keytype,
        transcript_id=self.transcript_id,
        timeserieslist=self.timeserieslist,
        symbs=self.symbs)

    def _get_gene_level(self):
        genes = defaultdict(list)
        for a in range(len(self.gene_id)):
            genes[self.gene_id[a]].append(a)

        self.genelevel_expression = np.empty(shape=(self.reps1,0,self.timepts))
        self.genelevel_id = []
        self.genelevel_symb = []
        for u,v in genes.items():
            tmp = np.expand_dims(np.sum(self.timeserieslist[:,v,:], axis=1),axis=1)
            self.genelevel_expression = np.concatenate([self.genelevel_expression, tmp], axis=1)
            self.genelevel_id.append(u)
            self.genelevel_symb.append(self.symbs[v[0]])
        
    def _remove_objects(self, objects_to_remove):
        if not isinstance(objects_to_remove, list):
            objects_to_remove = list(objects_to_remove)

        #remove objects from entries
        if self.transcript_id is not None:
            self.transcript_id = np.delete(np.array(self.transcript_id), objects_to_remove)

        self.ts[0] = np.delete(np.array(self.ts[0]), objects_to_remove, axis=0)
        tmp_genelist = np.delete(np.array(self.gene_id), objects_to_remove)
        self.gene_id = tmp_genelist
        tmp_symbs = np.delete(np.array(self.symbs), objects_to_remove)
        self.symbs = tmp_symbs
        tmp_timeseries_id = np.delete(np.arange(np.array(self.timeserieslist).shape[1]), objects_to_remove, axis=0)
        tmp_timeseries = [[self.timeserieslist[x][y] for y in tmp_timeseries_id] for x in range(np.array(self.timeserieslist).shape[0])]
        self.timeserieslist = np.array(tmp_timeseries, dtype="double")

        #recalculate gene level expression
        self._get_gene_level()

    def _add_clustering_object(self, clusterobj):
        self.clusterobj = clusterobj

    def _add_iso_switch_obj(self, isoobj):
        self.isoobj = isoobj

    
