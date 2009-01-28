# -*- coding: utf-8 -*-
#
# Python documentation build configuration file
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).

import sys, os, time
sys.path.extend([os.path.abspath('..'), os.path.abspath('tools/jinja2')])

# General configuration
# ---------------------

extensions = ['sphinx.ext.refcounting', 'sphinx.ext.coverage',
              'sphinx.ext.doctest', 'sphinx.ext.pngmath']
templates_path = ['tools']

# General substitutions.
project = 'simuPOP'
copyright = '2004-%s, Bo Peng' % time.strftime('%Y')

# The default replacements for |version| and |release|.
#
import simuPOP_version
version, release = simuPOP_version.SIMUPOP_VER, simuPOP_version.SIMUPOP_VER

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# List of files that shouldn't be included in the build.
unused_docs = [
]

# Relative filename of the reference count data file.
#refcount_file = 'data/refcounts.dat'

# If true, '()' will be appended to :func: etc. cross-reference text.
# This is set to false because the document already had (), and sometimes
# with options
add_function_parentheses = False

master_doc = 'index'

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

pngmath_use_preview = True

# Options for HTML output
# -----------------------

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Custom sidebar templates, filenames relative to this file.
html_sidebars = {
    'index': 'indexsidebar.html',
}

# Additional templates that should be rendered to pages.
html_additional_pages = {
    #'download': 'download.html',
    'index': 'indexcontent.html'
}

# Output an OpenSearch description file.
html_use_opensearch = '' #'http://docs.python.org/dev'

# Additional static files.
html_static_path = ['tools']

html_style = 'simuPOP.css'

# Output file base name for HTML help builder.
htmlhelp_basename = 'python' + release.replace('.', '')

# Split the index
html_split_index = True


# Options for LaTeX output
# ------------------------

# The paper size ('letter' or 'a4').
latex_paper_size = 'a4'

# The font size ('10pt', '11pt' or '12pt').
latex_font_size = '10pt'

# Documents to append as an appendix to all manuals.
latex_appendices = [] #'glossary', 'about', 'license', 'copyright']

# Options for the coverage checker
# --------------------------------

# The coverage checker will ignore all modules/functions/classes whose names
# match any of the following regexes (using re.match).
coverage_ignore_modules = [
]

coverage_ignore_functions = [
]

coverage_ignore_classes = [
]

# Glob patterns for C source files for C API coverage, relative to this directory.
coverage_c_path = [
    '../src/*.h',
]

# Regexes to find C items in the source files.
coverage_c_regexes = {
    'cfunction': (r'^PyAPI_FUNC\(.*\)\s+([^_][\w_]+)'),
    'data': (r'^PyAPI_DATA\(.*\)\s+([^_][\w_]+)'),
    'macro': (r'^#define ([^_][\w_]+)\(.*\)[\s|\\]'),
}

# The coverage checker will ignore all C items whose names match these regexes
# (using re.match) -- the keys must be the same as in coverage_c_regexes.
coverage_ignore_c_items = {
#    'cfunction': [...]
}
