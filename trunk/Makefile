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
	cd doc; python runSampleCode.py; make dist_svn

rebuild3:
	python3 setup.py install

clean:
	@rm -f src/*wrap* src/simuPOP_*.py
	@find . -name '*~' -delete
	@find . -name '*.bak' -delete
	@find . -name '*.py?' -delete
	@rm -f MANIFEST doc/*.aux doc/*.log doc/*.idx

