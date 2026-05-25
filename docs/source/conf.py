# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'santex'
copyright = '2024, Utpal Singh'
author = 'Utpal Singh, Sinan Özaydin, Vasileios Chatzaras, Patrice Rey'
release = '1.2.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['_build', '**.ipynb_checkpoints']

extensions = [
   'myst_parser',
   'nbsphinx',
   'sphinx_copybutton',
   'sphinx.ext.duration',
   'sphinx.ext.doctest',
   'sphinx.ext.autodoc',
   'sphinx.ext.autosummary',
   'sphinx.ext.mathjax',
]

# MyST-Parser settings — allow Markdown math and anchored headings
myst_enable_extensions = [
    "dollarmath",
    "colon_fence",
]
myst_heading_anchors = 3


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

def setup(app):
    app.add_css_file('custom.css') 
