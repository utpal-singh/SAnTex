# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import sphinx_book_theme
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))


project = 'santex'
copyright = '2024, Utpal Singh, Sinan Ozaydin, Vasileios Chatzaras, Patrice Rey'
author = 'Utpal Singh, Sinan Ozaydin, Vasileios Chatzaras, Patrice Rey'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx_design',
    'myst_nb',
    'jupyterlite_sphinx',
    'sphinx.ext.viewcode',
]




templates_path = ['_templates']
exclude_patterns = []

autosummary_generate = True
autosummary_imported_members = True
numpydoc_show_class_members = True
class_members_toctree = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
