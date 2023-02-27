[![PyPI version](https://badge.fury.io/py/spycone.svg)](https://badge.fury.io/py/spycone) 
[![Documentation Status](https://readthedocs.org/projects/spycone/badge/?version=latest)](https://spycone.readthedocs.io/en/latest/?badge=latest)

# Spycone - SPlicing-aware time-COurse Network Enricher

Spycone is a python package that provides systematic analysis of time course transcriptomics data. Spycone uses gene or isoform expression as an input. Spycone features a novel method for IS detection and employs the sum of changes of all isoforms relative abundances (total isoform usage) across time points. Spycone provides downstream analysis such as clustering by total isoform usage, i.e. grouping genes that are most likely to be coregulated, and network enrichment, i.e. extracting subnetworks or pathways that are over-represented by a list of genes. These analyses can be coupled with gene set enrichment analysis and visualization.

The paper is now accepted in [Bioinformatics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac846/6965022). 

# prerequisite

Spycone is dependent on the pcst_fast library, which is not available through pip install. Please go to [the github page](https://github.com/fraenkel-lab/pcst_fast) or run this command.
```
pip install https://github.com/fraenkel-lab/pcst_fast/archive/refs/tags/1.0.7.tar.gz
```

# Installation
It is available as a pypi package:
```
pip install spycone
```

Or alternatively install the latest development:
```
git clone https://github.com/yollct/spycone
cd spycone
pip install .
```

For more information, please check our documentation [https://spycone.readthedocs.io/en/latest/index.html](https://spycone.readthedocs.io/en/latest/index.html).
