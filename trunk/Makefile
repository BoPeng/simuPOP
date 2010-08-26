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

clean:
	@find . -name '*~' -depth 1 -delete
	@find . -name '*.bak' -depth 1 -delete
	@find . -name '*.py?' -depth 1 -delete
	@rm -f MANIFEST doc/*.aux doc/*.log doc/*.idx

