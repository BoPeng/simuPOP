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
	2to3 -w simuOpt.py  # 2to3-3.1 does not work for this file under mac
	2to3-3.1 -w src/__init__.py src/sampling.py src/plotter.py src/utils.py
	2to3-3.1 -w tools/build.py tools/doxy2swig.py
	2to3-3.1 -w test\run_tests.py               test\test_05_matings.py         test\test_11_terminator.py      test\test_17_utils.py$
	2to3-3.1 -w test\test_00_genoStru.py        test\test_06_initialization.py  test\test_12_migration.py       test\test_18_plotter.py$
	2to3-3.1 -w test\test_01_individual.py      test\test_07_tagging.py         test\test_13_mutation.py        test\test_19_performance.py$
	2to3-3.1 -w test\test_02_population.py      test\test_08_stat.py            test\test_14_transmitter.py$
	2to3-3.1 -w test\test_03_operator.py        test\test_09_selection.py       test\test_15_penetrance.py$
	2to3-3.1 -w test\test_04_simulator.py       test\test_10_qtrait.py          test\test_16_sampling.py$
	2to3-3.1 -w doc/userGuide.py doc/runSampleCode.py

build3: 2to3
	python3 setup.py install

clean:
	@rm -f src/*wrap* src/simuPOP_*.py
	@find . -name '*~' -delete
	@find . -name '*.bak' -delete
	@find . -name '*.py?' -delete
	@rm -f MANIFEST doc/*.aux doc/*.log doc/*.idx

