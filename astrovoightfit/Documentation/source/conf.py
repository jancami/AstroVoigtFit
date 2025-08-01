# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AstroVoigtFit'
copyright = '2025, Edibles'
author = 'Edibles'
release = '2025'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',        # Enables the autofunction directive
    'sphinx.ext.napoleon',       # Supports Google/NumPy-style docstrings
    'sphinx.ext.viewcode',       # Adds source code links
    'sphinx.ext.mathjax',        # Keep this if you're using math
    'sphinxemoji.sphinxemoji',   # Optional emoji support
    'sphinx_rtd_theme'           # Your selected HTML theme
]


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'English'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


import os
import sys
sys.path.insert(0, os.path.abspath('../../utils'))