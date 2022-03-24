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
sys.path.insert(0, os.path.abspath('..'))
#sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'dNami'
copyright = '2019, Conservatoire National des Arts et Métiers and Okinawa Institute of Science and Technology School Corporation'
author = 'Nicolas Alferez, S2T Unit'

latex_engine = 'xelatex'
latex_use_xindy = False
# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
  "sphinx_rtd_theme",
  "sphinx.ext.autodoc",
  "sphinx.ext.coverage", 
  "sphinx.ext.napoleon",
  "sphinxcontrib.bibtex",
]

# Citations file
bibtex_bibfiles = ['refs.bib']

#"pydata_sphinx_theme"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','./usage/examples.rst']

#exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','./usage/clusters.rst',
#        './usage/deigo.rst','./usage/fugaku.rst','./usage/fugaku_shared_storage.rst',
#       './usage/ko_fugaku.rst']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'pydata_sphinx_theme'
html_theme = 'sphinx_rtd_theme'
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# -- Make sure figures are numbered in references
numfig = True

autodoc_mock_imports = ["numpy","dnami","dnamiF","dnami_io","dnami_mpi"]
