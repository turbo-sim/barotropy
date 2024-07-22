import os
import sys

# Add package root dir to path
sys.path.insert(0, os.path.abspath(".."))

# Define project metadata
project = "barotropy"
copyright = "2023, Sustainable Thermal Power DTU"
author = "Roberto Agromayor"
release = "v0.0.0"

# Define extensions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "numpydoc",
    "sphinx.ext.todo",
    # 'sphinx_tabs.tabs',
    "sphinx_togglebutton",
    "sphinx_design",
]

todo_include_todos = True

# Add bibliography file
bibtex_bibfiles = ["source/bibliography.bib"]
bibtex_default_style = "alpha"
bibtex_reference_style = "author_year"

# Exclude unnecessary files
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
exclude_patterns.extend(["source/PySolverView.rst"])


# Define theme
# html_theme = "sphinx_book_theme"
# html_theme = 'pydata_sphinx_theme'
html_theme = 'sphinx_rtd_theme'


