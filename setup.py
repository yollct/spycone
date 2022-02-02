from setuptools import setup, find_packages, dist
import os
from distutils.core import Extension
from distutils.command.sdist import sdist as _sdist
from setuptools import Extension
import sys




def configuration(parent_package="", top_path=None):
    if os.path.exists("MANIFEST"):
        os.remove("MANIFEST")

    from numpy.distutils.misc_util import Configuration
    import numpy
    config = Configuration("spycone", parent_package, top_path)

    # Avoid useless msg:
    # "Ignoring attempt to set 'name' (from ... "
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True,
    )

    # Cython is required by config.add_subpackage for templated extensions
    # that need the tempita sub-submodule. So check that we have the correct
    # version of Cython so as to be able to raise a more informative error
    # message from the start if it's not the case.

    config.add_subpackage("spycone/DOMINO")
    config.add_subpackage("spycone/DOMINO/src/pcst_fast_py")
    config.add_subpackage("spycone/_feature")
    config.add_subpackage("spycone/_fitness")
    config.add_subpackage("spycone/_NEASE")
    config.add_subpackage("spycone/_preprocessing")
    config.add_subpackage("spycone/_prototype")
    config.add_subpackage("spycone/_util_stat")
    config.add_subpackage("spycone/_visualization")

    config.add_extension("spycone._clustering.similarities", sources=["spycone/_clustering/similarities.c"], include_dirs=[numpy.get_include()])
    config.add_extension("spycone._iso_function._iso_fast", sources=["spycone/_clustering/similarities.c"], include_dirs=[numpy.get_include()])
    config.add_extension("spycone._connectivity.connectivity_task", ["spycone/_connectivity/connectivity_task.c"],include_dirs = [numpy.get_include()])
    config.add_extension("spycone._shuffle.shuffle", ["spycone/_shuffle/shuffle.c"],include_dirs = [numpy.get_include()])

    return config

def setup_package():
    metadata = dict(name='spycone',
    version='0.0.1',
    description='A splicing-aware time course network enricher',
    url='https://gitlab.lrz.de/yollct/ticone_pkg',
    author='Chit Tong Lio',
    author_email='chit-tong.lio@wzw.tum.de',
    license='GPLv3',
    classifiers=[
      "Programming Language :: Python :: 3",
      "Programming Language :: Python :: 3.8",
      "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
      "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.8',
    setup_requires=['pybind11>=2.1.0'],
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
      'pybind11==2.1.0'])
  
    import numpy
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
    from setuptools import setup

    cmdclass = {}
    ext_modules = []
    ext_modules += [
        Extension("spycone._clustering.similarities", sources=["spycone/_clustering/similarities.c"], include_dirs=[numpy.get_include()]),
        Extension("spycone._iso_function._iso_fast", sources=["spycone/_clustering/similarities.c"], include_dirs=[numpy.get_include()]),
        Extension("spycone._connectivity.connectivity_task", ["spycone/_connectivity/connectivity_task.c"],include_dirs = [numpy.get_include()]),
        Extension("spycone._shuffle.shuffle", ["spycone/_shuffle/shuffle.c"],include_dirs = [numpy.get_include()]),
    ]
    
    cmdclass.update({'build_ext': build_ext})

    metadata['setup_requires']=['numpy']
    metadata['ext_modules'] = ext_modules
    metadata["configuration"] = configuration

    setup(**metadata)


if __name__ == "__main__":
  

  setup_package()
