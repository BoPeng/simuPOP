2to3.py -w setup.py
2to3.py -w simuOpt.py 
2to3.py -w src\__init__.py src\sampling.py src\plotter.py src\utils.py
2to3.py -w tools\build.py tools\doxy2swig.py
2to3.py -w test\run_tests.py               test\test_05_matings.py         test\test_11_terminator.py      test\test_17_utils.py
2to3.py -w test\test_00_genoStru.py        test\test_06_initialization.py  test\test_12_migration.py       test\test_18_plotter.py
2to3.py -w test\test_01_individual.py      test\test_07_tagging.py         test\test_13_mutation.py        test\test_19_performance.py
2to3.py -w test\test_02_population.py      test\test_08_stat.py            test\test_14_transmitter.py
2to3.py -w test\test_03_operator.py        test\test_09_selection.py       test\test_15_penetrance.py
2to3.py -w test\test_04_simulator.py       test\test_10_qtrait.py          test\test_16_sampling.py
2to3.py -w doc\userGuide.py doc\runSampleCode.py
