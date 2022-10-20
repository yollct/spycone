# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / 'extensions')]
sys.path.insert(0, os.path.abspath('../../'))
for x in os.walk('../../spycone'):
  sys.path.insert(0, x[0])
on_rtd = os.environ.get('READTHEDOCS') == 'True'

import spycone 

# -- Project information -----------------------------------------------------

project = 'spycone'
copyright = '2021, Chit Tong Lio'
author = 'Chit Tong Lio'

# The full version, including alpha/beta/rc tags
release = '0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    "sphinx_rtd_theme",
    "nbsphinx",
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.ifconfig', 'sphinx.ext.viewcode',
    'sphinx.ext.inheritance_diagram',
    ]

autosummary_generate = True
autodoc_member_order = 'bysource'
autodoc_default_flags = ['members']
# napoleon_google_docstring = False
# napoleon_numpy_docstring = True
# napoleon_include_init_with_doc = False
# napoleon_use_rtype = True  # having a separate entry generally helps readability
# napoleon_use_param = True
# napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False
api_dir = 'api'  # function_images

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '.idea', '.DS_Store', '.ipynb_checkpoints', '__pycache__', 'data']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


suppress_warnings = [
    'nbsphinx.localfile',
    'nbsphinx.gallery',
    'nbsphinx.thumbnail',
    'nbsphinx.notebooktitle',
    'nbsphinx.ipywidgets',
]