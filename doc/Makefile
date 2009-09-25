#
# Makefile for Python documentation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# You can set these variables from the command line.

ifeq ($(HOSTNAME), bp8)
PYTHON       = $(HOME)/Python26/bin/python
else
PYTHON       = python
endif
SVNROOT      = http://svn.python.org/projects
SPHINXOPTS   =
PAPER        =
SOURCES      =
DISTVERSION  =

ALLSPHINXOPTS = -b $(BUILDER) -d build/doctrees -D latex_paper_size=$(PAPER) \
                $(SPHINXOPTS) . build/$(BUILDER) $(SOURCES)

.PHONY: help checkout update build html htmlhelp clean coverage dist

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  text      to make plain text files"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  coverage  to check documentation coverage for library and C API"
	@echo "  dist      to create a \"dist\" directory with archived docs for download"

checkout:
	@if [ ! -d tools/sphinx ]; then \
	  echo "Checking out Sphinx..."; \
	  svn checkout $(SVNROOT)/doctools/trunk/sphinx tools/sphinx; \
	fi
	@if [ ! -d tools/docutils ]; then \
	  echo "Checking out Docutils..."; \
	  svn checkout $(SVNROOT)/external/docutils-0.5/docutils tools/docutils; \
	fi
	@if [ ! -d tools/jinja2 ]; then \
	  echo "Checking out Jinja2..."; \
	  hg  clone http://dev.pocoo.org/hg/jinja2-main tools/jinja2; \
	fi
	@if [ ! -d tools/pygments ]; then \
	  echo "Checking out Pygments..."; \
	  svn checkout $(SVNROOT)/external/Pygments-0.11.1/pygments tools/pygments; \
	fi
	@if [ ! -d tools/converter ]; then \
	  echo "Checking out converter..."; \
	  svn checkout $(SVNROOT)/doctools/converter/converter tools/converter; \
	fi

update: checkout
	svn update tools/sphinx
	svn update tools/docutils
	cd tools/jinja2; hg update; cd ../..
	svn update tools/pygments
	svn update tools/converter

build: checkout
	mkdir -p build/$(BUILDER) build/doctrees
	$(PYTHON) tools/sphinx-build.py $(ALLSPHINXOPTS)
	@echo

build/userGuide.tex: userGuide.lyx
	/usr/local/bin/lyx --export pdflatex userGuide.lyx
	mv userGuide.tex build/

build/refManual.tex: refManual.lyx
	/usr/local/bin/lyx --export pdflatex refManual.lyx
	mv refManual.tex build/

rst: build/userGuide.tex build/refManual.tex
	$(PYTHON) tools/convert.py build/userGuide.tex build/userGuide.rst
	$(PYTHON) tools/convert.py build/refManual.tex build/refManual.rst

html: BUILDER = html
html: rst build
	@echo "Build finished. The HTML pages are in build/html."

htmlhelp: BUILDER = htmlhelp
htmlhelp: build
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      "build/htmlhelp/pydoc.hhp project file."

latex: BUILDER = latex
latex: build
	@echo "Build finished; the LaTeX files are in build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

text: BUILDER = text
text: build
	@echo "Build finished; the text files are in build/text."

changes: BUILDER = changes
changes: build
	@echo "The overview file is in build/changes."

linkcheck: BUILDER = linkcheck
linkcheck: build
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/$(BUILDER)/output.txt"

coverage: BUILDER = coverage
coverage: build
	@echo "Coverage finished; see c.txt and python.txt in build/coverage"

doctest: BUILDER = doctest
doctest: build
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in build/doctest/output.txt"

pydoc-topics: BUILDER = pydoc-topics
pydoc-topics: build
	@echo "Building finished; now copy build/pydoc-topics/pydoc_topics.py " \
	      "into the Lib/ directory"

htmlview: html
	 $(PYTHON) -c "import webbrowser; webbrowser.open('build/html/index.html')"

clean:
	-rm -rf build/*
	-rm -rf tools/sphinx

src_doc:
	cd ..; \
	SIMUPOP_VER=`cut -d'"' -f2 simuPOP_version.py | head -1`; export SIMUPOP_VER; \
	SIMUPOP_REV=`cut -d'"' -f2 simuPOP_version.py | tail -1`; export SIMUPOP_REV; \
	doxygen Doxy_web
	rsync -v -rlt -z ../doxygen_doc/html/ simupop,simupop@web.sourceforge.net:htdocs/src_doc

dist: html
	# remove old directory
	rm -rf build/empty
	mkdir build/empty
	cp build/example_*.py build/html/build
	rsync --delete --archive build/empty/ simupop,simupop@web.sourceforge.net:htdocs/manual/build
	rsync -v -rlt -z build/html/ simupop,simupop@web.sourceforge.net:htdocs/manual

pdf: userGuide.lyx refManual.lyx
	/usr/local/bin/lyx --export pdf2 userGuide.lyx
	/usr/local/bin/lyx --export pdf2 refManual.lyx