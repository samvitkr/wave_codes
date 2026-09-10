# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "alps"
copyright = ""
author = ""

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "breathe",
    "sphinx.ext.mathjax",
    "sphinx.ext.graphviz",
    'myst_parser',
    "exhale",
]

# Setup the breathe extension
breathe_projects = {"alps": "./_doxygen/xml"}
breathe_default_project = "alps"

# Setup the MyST parser
myst_heading_anchors = 3 # generate anchors for h1, h2, and h3 level headings
myst_ref_domains = ["cpp"]
myst_enable_extensions = [ 'dollarmath' ]

graphviz_output_format = "svg"

import textwrap

doxy_configs = textwrap.dedent(
    """
    TAB_SIZE           = 2
    EXTRACT_ALL            = YES
    EXTRACT_PRIVATE        = YES
    EXTRACT_PRIV_VIRTUAL   = NO
    EXTRACT_PACKAGE        = YES
    EXTRACT_STATIC         = YES
    EXTRACT_LOCAL_CLASSES  = YES
    EXTRACT_ANON_NSPACES   = NO
    RESOLVE_UNNAMED_PARAMS = YES
    HIDE_FRIEND_COMPOUNDS  = NO
    GENERATE_TODOLIST      = NO
    INPUT = ../src
    EXCLUDE_SYMBOLS        = Kokkos, HighFive
    SEARCH_INCLUDES        = YES
    MACRO_EXPANSION        = YES
    INCLUDE_PATH = ../thirdparty/kokkos/core/src
    INCLUDE_PATH += ../thirdparty/kokkos/containers/src
    INCLUDE_PATH += ../thirdparty/enum/headers
    PREDEFINED += KOKKOS_ENABLE_CUDA
    PREDEFINED += KOKKOS_ENABLE_OPENMP
    PREDEFINED += KOKKOS_INLINE_FUNCTION=inline
    PREDEFINED += KOKKOS_FUNCTION=
    PREDEFINED += KOKKOS_LAMBDA=[=]
    CLASS_GRAPH            = YES
    COLLABORATION_GRAPH    = YES
    GROUP_GRAPHS           = YES
    CALL_GRAPH             = NO
    CALLER_GRAPH           = NO
    INCLUDE_GRAPH          = NO
    INCLUDED_BY_GRAPH      = NO
    DIRECTORY_GRAPH        = NO
    TEMPLATE_RELATIONS     = YES
    HIDE_UNDOC_RELATIONS   = NO
    HAVE_DOT               = YES
    DOT_IMAGE_FORMAT       = svg
"""
)


# custom specifications for Breathe directives
def specificationsForKind(kind):
    """
    For a given input ``kind``, return the list of reStructuredText specifications
    for the associated Breathe directive.
    """
    # Change the defaults for .. doxygenclass:: and .. doxygenstruct::
    if kind == "class" or kind == "struct":
        return [
            ":members:",
            ":protected-members:",
            ":private-members:",
            ":undoc-members:",
            ":allow-dot-graphs:",
        ]
    # An empty list signals to Exhale to use the defaults
    else:
        return []


from exhale import utils

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder": "./api",
    "rootFileName": "EXCLUDE",
    "classHierarchyFilename": "class_view_hierarchy.rst",
    "fileHierarchyFilename": "file_view_hierarchy.rst",
    "unabridgedApiFilename": "unabridged_api.rst",
    "doxygenStripFromPath": "../src",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle": "API reference",
    "fullApiSubSectionTitle": "Full reference",
    # optional arguments
    "createTreeView": True,
    "minifyTreeView": False,
    "fullToctreeMaxDepth": 1,
    "contentsDirectives": False, # furo theme does not need contents list at the top of pages
    "customSpecificationsMapping": utils.makeCustomSpecificationsMapping(
        specificationsForKind
    ),
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": doxy_configs,
}

# Tell sphinx what the primary language being documented is.
primary_domain = "cpp"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "cpp"


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", ".venv"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

pygments_style = "default"
pygments_dark_style = "monokai"

html_theme = "furo"
html_static_path = ["_static"]
# html_theme_options = {
#     "navigation_depth": 3,
# }