# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import pathlib
import sys

os.environ.setdefault('NUMBA_DISABLE_JIT', '1')

sys.path.insert(0, str(pathlib.Path("..").resolve()))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'atmPy'
copyright = '2023, Hagen Telg'
author = 'Hagen Telg'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'nbsphinx']

autodoc_mock_imports = [
    'atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cal_facility',
    'atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cal_facility.lab',
    'netCDF4',
]

nbsphinx_execute = 'never'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'#'alabaster'
html_static_path = ['_static']
