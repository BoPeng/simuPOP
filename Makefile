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
	2to3-3.1 -w setup.py simuOpt.py
	2to3-3.1 -w src/sampling.py src/plotter.py src/utils.py
	2to3-3.1 -w tools/build.py tools/doxy2swig.py
	2to3-3.1 -w test/*.py
	2to3-3.1 -w doc/userGuide.py doc/runSampleCode.py

clean:
	@find . -name '*~' -delete
	@find . -name '*.bak' -delete
	@find . -name '*.py?' -delete
	@rm -f MANIFEST doc/*.aux doc/*.log doc/*.idx

