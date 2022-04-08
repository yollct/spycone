from setuptools import setup, find_packages, dist
import os
from distutils.core import Extension
from distutils.command.sdist import sdist as _sdist
from setuptools import Extension
import sys


def setup_package():
    metadata = dict(name='spycone',
    version='0.0.2',
    description='A splicing-aware time course network enricher',
    url='https://github.com/yollct/spycone.git',
    author='Chit Tong Lio',
    author_email='ct7298@gmail.com',
    license='GPLv3',
    classifiers=[
      "Programming Language :: Python :: 3",
      "Programming Language :: Python :: 3.7",
      "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
      "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires=[
      'pandas>=1.0.1',
      'numpy>=1.18.1',
      'seaborn',
      'scikit-learn>=0.23.2',
      'scikit-learn-extra>=0.1.0',
      'networkx>=2.4',
      'matplotlib>=3.1.3',
      'gseapy>=0.10.4',
      'scipy>=1.4.1',
      'statsmodels>=0.11.0',
      'tslearn>=0.5.1.0',
      'python-louvain',
      'plotly>=4.14.3',
      'biopython',
      'gtfparse', 
      'pyvis>=0.1.9', 
      'joblib',
      'nease'])
  
    import numpy
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
    import pybind11
    from setuptools import setup

    class get_pybind_include(object):
      """Helper class to determine the pybind11 include path

      The purpose of this class is to postpone importing pybind11
      until it is actually installed, so that the ``get_include()``
      method can be invoked. """

      def __init__(self, user=False):
          self.user = user

      def __str__(self):
          import pybind11
          return pybind11.get_include(self.user)

    cmdclass = {}
    ext_modules = []
    ext_modules += [
    #     #Extension("spycone._clustering.similarities", sources=["spycone/_clustering/similarities.c"], include_dirs=[numpy.get_include()]),
    #     #Extension("spycone._iso_function._iso_fast", sources=["spycone/_clustering/similarities.c"], include_dirs=[numpy.get_include()]),
          # Extension("spycone._connectivity.connectivity_task", ["spycone/_connectivity/connectivity_task.c"],include_dirs = [numpy.get_include()]),
          # Extension("spycone._shuffle.shuffle", ["spycone/_shuffle/shuffle.c"],include_dirs = [numpy.get_include()]),
        #   Extension(
        # 'pcst_fast',
        # sources=['spycone/DOMINO/src/pcst_fast_py/src/pcst_fast_pybind.cc', 'spycone/DOMINO/src/pcst_fast_py/src/pcst_fast.cc'],
        # include_dirs=[
        #     # Path to pybind11 headers
        #     get_pybind_include(),
        #     get_pybind_include(user=True)
        # ],
        # language='c++'
    
    ]
    
    cmdclass.update({'build_ext': build_ext})

    metadata['setup_requires']=['numpy']
    metadata['ext_modules'] = ext_modules

    setup(**metadata)


if __name__ == "__main__":
  setup_package()
