from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("compute_pvalues", 
    ["compute_pvalues.pyx"]#include_dirs = [numpy.get_include()]
    )
]

setup(
    ext_modules = cythonize(extensions, gdb_debug=True)
    #cmdclass = {'build_ext':build_ext}
)
