from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
import numpy

extensions = [
    Extension("shuffle", 
    ["shuffle.pyx"],
    include_dirs = [numpy.get_include()]
    )
]

setup(
    ext_modules = cythonize(extensions),
    cmdclass = {'build_ext':build_ext}
    
)