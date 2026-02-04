#
# Makefile for processing simuPOP
#


rebuild:
	bcpp src/simuPOP/*.h src/simuPOP/*.cpp
	doxygen
	tools/doxy2swig.py
	scons install -j4 std
	tools/doxy2swig.py
	scons install -j4
	cd doc; python runSampleCode.py; make dist_svn

rebuild3:
	pip install .

clean:
	@rm -rf build _skbuild
	@rm -f src/simuPOP/*wrap* src/simuPOP/simuPOP_*.py
	@find . -name '*~' -exec rm -f {} \;
	@find . -name '*.bak' -exec rm -f {} \;
	@find . -name '*.py?' -exec rm -f {} \;
	@rm -f MANIFEST doc/*.aux doc/*.log doc/*.idx

