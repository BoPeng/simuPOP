#
# Makefile for processing simuPOP
#

rebuild:
	bcpp src/*.h src/*.cpp
	doxygen
	tools/doxy2swig.py
	scons install -j4 std
	tools/doxy2swig.py
	scons install -j4
	cd doc
	#make pdf
	python runSampleCode.py userGuide.py
	make dist_svn

2to3:
	2to3-3.1 -w setup.py
	2to3-3.1 -w src/sampling.py
	2to3-3.1 -w src/plotter.py
	2to3-3.1 -w src/utils.py
	2to3-3.1 -w tools/build.py
	2to3-3.1 -w tools/doxy2swig.py
	2to3-3.1 -w test/*.py

clean:
	@find . -name '*~' -delete
	@find . -name '*.bak' -delete
	@find . -name '*.py?' -delete
	@rm -f MANIFEST doc/*.aux doc/*.log doc/*.idx

