from numpy.distutils.core import setup
import os
import sys

# def configuration(parent_package="", top_path=None):
#     from numpy.distutils.misc_util import Configuration
#     import numpy

#     libraries = []
#     if os.name == "posix":
#         libraries.append("m")

#     config = Configuration("spycone", parent_package, top_path)

#     config.add_subpackage("_clustering")
#     config.add_subpackage("_connectivity")
#     config.add_subpackage("_clustering")
#     config.add_subpackage("_clustering")

#     return config

# if __name__ == "__main__":
#     setup(**configuration(top_path="").todict())