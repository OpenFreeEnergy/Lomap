# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

project = 'lomap'
copyright = '2023, Lomap developers'
author = 'Lomap developers'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx_toolbox.collapse',
    'sphinx.ext.autosectionlabel',
    'sphinx_design',
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3.9", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    "rdkit": ("https://www.rdkit.org/docs", None),
    "openff.units": ("https://docs.openforcefield.org/units/en/stable", None),
    "gufe": ("https://gufe.readthedocs.io/en/latest/", None),
}


autoclass_content = 'both'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['ofe_sphinx_theme']
html_theme_options = {
    "accent_color": "FeelingBlue",
}
html_css_files = []
