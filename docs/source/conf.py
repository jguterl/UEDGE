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
sys.path.insert(0,'/Users/guterl/opt/anaconda3/lib/python3.7/site-packages/')
sys.path.insert(0,'/Users/guterl/opt/anaconda3/lib/python3.7/site-packages/')
sys.path.insert(0,'/home/jguterl/UEDGE/docs/source/')
#sys.path.insert(0,'/home/jguterl/UEDGE/docs/source/DocPy/')
#sys.path.insert(0, os.path.abspath('../../pyscripts/'))
print(sys.path)

# -- Project information -----------------------------------------------------

project = 'UEDGE'
copyright = '2020, J.Guterl'
author = 'J.Guterl'

# The full version, including alpha/beta/rc tags
release = '0.0.'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon','sphinx.ext.viewcode','sphinxfortran.fortran_domain','sphinxfortran.fortran_autodoc']

fortran_src=['/Users/guterl/git/UEDGE/bbb/boundary.F','/Users/guterl/git/UEDGE/bbb/boundary.F90']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

add_function_parentheses = True
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_theme_path = ["_themes", ]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']