# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath('../HMcode2020Emu/HMcode2020Emu'))
autodoc_mock_imports = ['HMcode2020Emu']


# -- Project information -----------------------------------------------------

project = 'HMcode2020Emu'
copyright = '2023, M. Tsedrik'
author = 'M. Tsedrik'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx_rtd_theme'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'groundwork'
html_theme = 'sphinx_rtd_theme'
# html_theme = 'p-main_theme'
# html_theme = 'bootstrap'

# only for p-main_theme, from here
# import os
# from PSphinxTheme import utils
# p, html_theme, needs_sphinx = utils.set_psphinxtheme(html_theme)
# html_theme_path = p
# util here

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']