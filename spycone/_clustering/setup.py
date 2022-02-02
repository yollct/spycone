from setuptools import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import os


###sklearn cythonize
def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    import numpy
    config = Configuration("_clustering", parent_package, top_path)

    config.add_extension(
        "similarities",
        sources = ['similarities.pyx'],
        include_dirs = [numpy.get_include()]
    )

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration(top_path="").todict())


