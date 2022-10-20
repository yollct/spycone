from .BioNetwork import BioNetwork
from .DataSet import dataset
from .clustering import clustering
from .iso import iso_function
from .preprocess import preprocess
from .visualize import vis_all_clusters,  switch_plot, gsea_plot, vis_modules, vis_better_modules, compare_mod
from .go_terms import clusters_gsea, modules_gsea, list_gsea
#from .connectivity import connectivity
from .run_domino import run_domino, run_domain_domino
from .DOMINO.src.core import domino
from .splicingfactor import SF_coexpression, SF_motifsearch
#from ._NEASE import nease