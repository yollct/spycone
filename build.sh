#!/bin/bash

cd spycone
cd _connectivity
python setup.py build_ext --inplace
cd ..
# cd similarity
# python setup.py build_ext --inplace
# cd ..
cd _clustering
python setup.py build_ext --inplace
cython similarities.pyx 
cd ..
# cd _util_stat
# python setup.py build_ext --inplace
# cython compute_pvalues.pyx -a
# cd ..
cd _shuffle
python setup.py build_ext --inplace
cython shuffle.pyx 
cd ..

cd _iso_function
python setup.py build_ext --inplace
cython _iso_fast.pyx 
cd ..
cd ..
