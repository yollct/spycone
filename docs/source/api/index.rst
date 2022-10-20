API
====

Import:

.. code-block:: python

    import spycone as spy


Import dataset
--------
.. currentmodule:: spycone

.. autoclass:: dataset

Import and store biological network
------------
.. autoclass:: BioNetwork

Preprocessing
-----------

.. autoclass::
    preprocess
        

Isoform-level function
-------------

.. autoclass:: iso_function


Perfom clustering
------------

.. autoclass:: clustering

GO terms enrichment analysis
----------

.. autofunction:: list_gsea
.. autofunction:: clusters_gsea
.. autofunction:: modules_gsea
.. autofunction:: list_genesets
    
 
Run DOMINO
-----------

.. autofunction:: run_domino
.. autofunction:: run_domain_domino


Visualization
-----------

.. autofunction:: vis_all_clusters
.. autofunction:: switch_plot
.. autofunction:: gsea_plot
.. autofunction:: vis_modules
.. autofunction:: vis_better_modules


Splicing factor analysis
-------------
.. autofunction:: SF_coexpression
.. autofunction:: SF_motifsearch

    