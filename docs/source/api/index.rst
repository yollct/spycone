API
====

Import:

.. code-block:: python

    import spycone as spy

.. toctree::
   :caption: Contents:

Import dataset
--------------
.. currentmodule:: spycone

.. autoclass::
    dataset


Import and store biological network
------------

.. autoclass::
    BioNetwork

Preprocessing 
-----------

.. autoclass::
    preprocess


Isoform-level function
-------------

.. autoclass:: iso_function


Clustering alogrithms
------------

.. autoclass::
    clustering


Perform GO term enrichment
----------

.. autofunction:: list_gsea
.. autofunction:: clusters_gsea
.. autofunction:: modules_gsea
    
 
Run DOMINO 
-----------

.. autofunction:: run_domino
.. autofunction:: run_domain_domino


Visualization functions
-----------

.. autofunction:: vis_all_clusters
.. autofunction:: switch_plot
.. autofunction:: gsea_plot


    