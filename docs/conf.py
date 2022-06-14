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


# -- Project information -----------------------------------------------------

project = 'rotsim2d'
copyright = '2020-2022, Grzegorz Kowzan'
author = 'Grzegorz Kowzan'

with open('citations.txt', 'rt') as f:
    rst_prolog = f.read()

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
    'sphinx.ext.mathjax',
    'sphinx.ext.autosummary',
    'IPython.sphinxext.ipython_directive',
    'IPython.sphinxext.ipython_console_highlighting',
    'sphinx_copybutton'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
exclude_patterns = []

add_module_names = False
python_use_unqualified_type_names = True

# sphinx-copybutton configurations
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
copybutton_copy_empty_lines = False

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
# html_sidebars = {
#     '**': [
#         'about.html',
#         'navigation.html',
#         'searchbox.html',
#     ]
# }
# html_show_sourcelink = True
html_theme_options = dict(
    repository_url="https://github.com/gkowzan/rotsim2d",
    repository_branch="master",
    path_to_docs="docs",
    use_edit_page_button=True,
    use_repository_button=True,
    use_issues_button=True,
    home_page_in_toc=False
)

html_logo = "_static/rotsim2d_logo.png"
html_favicon = "_static/rotsim2d_logo.ico"

# -- Intersphinx
tls_cacerts = '/etc/ssl/certs/'
tls_verify = False
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'matplotlib': ('https://matplotlib.org/', None),
    'xarray': ('https://xarray.pydata.org/en/stable/', None),
    'molspecutils': ('https://molspecutils.readthedocs.io/en/latest/', None),
    'anytree': ('https://anytree.readthedocs.io/en/latest/', None)
}

# -- pygments
pygments_style = 'sphinx'

# -- napoleon
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_special_with_doc = False
napoleon_preprocess_types = True
napoleon_attr_annotations = True

# -- autodoc and autosummary
autosummary_generate = True
autodoc_typehints = 'description'
autodoc_typehints_description_target = "documented"
autoclass_content = 'class'
autodoc_preserve_defaults = True
autodoc_default_options = {
    "show-inheritance": True
}
