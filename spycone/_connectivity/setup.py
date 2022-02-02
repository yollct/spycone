from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
import numpy

extensions = [
    Extension("connectivity_task", 
    ["connectivity_task.pyx"],
    include_dirs = [numpy.get_include()]
    )
]

setup(
    ext_modules = cythonize(extensions, gdb_debug=True),
    cmdclass = {'build_ext':build_ext}
    
)
