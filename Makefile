#
# Makefile for processing simuPOP
#


rebuild:
	bcpp simuPOP/*.h simuPOP/*.cpp
	doxygen
	tools/doxy2swig.py
	scons install -j4 std
	tools/doxy2swig.py
	scons install -j4
	cd doc; python runSampleCode.py; make dist_svn

rebuild3:
	python3 setup.py install

clean:
	@rm -rf build
	@rm -f simuPOP/*wrap* simuPOP/simuPOP_*.py
	@find . -name '*~' -exec rm -f {} \;
	@find . -name '*.bak' -exec rm -f {} \;
	@find . -name '*.py?' -exec rm -f {} \;
	@rm -f MANIFEST doc/*.aux doc/*.log doc/*.idx

