# docs/conf.py

#import os
#import sys

#sys.path.insert(
#    0,
#    os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
#)

project = "indelinside"
copyright = "2026, Kohei Hagiwara"
author = "Kohei Hagiwara"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
]

autosummary_generate = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"

html_static_path = ["_static"]

autodoc_member_order = "bysource"
napoleon_google_docstring = True
napoleon_numpy_docstring = True
