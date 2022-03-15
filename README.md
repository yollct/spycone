# Spycone - A python package for time series biological data

Spycone is a python package that provides systematic analysis of time course transcriptomics data. Spycone uses gene or isoform expression as an input. Spycone features a novel method for IS detection and employs the sum of changes of all isoforms relative abundances (total isoform usage) across time points. Spycone provides downstream analysis such as clustering by total isoform usage, i.e. grouping genes that are most likely to be coregulated, and network enrichment, i.e. extracting subnetworks or pathways that are over-represented by a list of genes. These analyses can be coupled with gene set enrichment analysis and visualization.

# prerequisite

Spycone is dependent on the pcst_fast library, which is not available through pip install. Please go to [the github page](https://github.com/fraenkel-lab/pcst_fast) or follow this instruction.
```
git clone https://github.com/fraenkel-lab/pcst_fast
cd pcst_fast
pip install .
```

# Installation
```
git clone https://github.com/yollct/spycone
cd spycone
pip install .
```


For more information, please check our documentation [https://spycone.readthedocs.io/en/latest/index.html](https://spycone.readthedocs.io/en/latest/index.html).
