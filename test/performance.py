#!/usr/bin/env python
#
# Purpose:
#    testing of performance of various evolutionary scenarios
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 475 $
# $LastChangedDate: 2006-10-04 01:00:38 -0500 (Wed, 04 Oct 2006) $
#
# The test results are recorded so we will be able to compare
# performance for different revisions.
#
# Usage:
#
#     performance.py [--alleleType=binary|short|long] [testname1] [testname2] ...
#
# where
#     testname can be any string within a test
#
#
#
import os, sys, time
from itertools import product

alleleType = None
# allele type can be specified by --alleleType=long/short/binary
if True in [x.startswith('--alleleType=') for x in sys.argv]:
    idx = [x.startswith('--alleleType=') for x in sys.argv].index(True)
    alleleType = sys.argv[idx][13:]
    sys.argv.pop(idx)

import simuOpt
simuOpt.setOptions(alleleType=alleleType, quiet=True, optimized=True)
from simuPOP import *

class PerformanceTest:
    def __init__(self, desc=''):
        self.description = desc
        self.expected_time = 0

    def describe(self):
        # machine name, os, time
        desc = 'HOST: %s (%s)\n' % (os.uname()[1], sys.platform)
        desc += 'TIME: %s\n' % time.asctime()
        info = moduleInfo()
        opt = ''
        if info['optimized']:
            opt = ', optimized'
        desc += 'SIMUPOP: %s (rev %d, module %s%s) on Python %s (%s)\n' % \
            (info['version'], info['revision'], info['alleleType'], opt, info['python'], info['compiler'])
        desc += 'TEST: %s' % self.description
        return desc

    def serialRun(self, **kwargs):
        '''
        Call self._run with a list of arguments, for example
            self.serialRun(a=[v1, v2], b=[v1, v2])
        will result in
            self._run(a=v1, b=v1)
            self._run(a=v2, b=v2)
        '''
        keys = kwargs.keys()
        for arg in zip(*[kwargs[x] for x in keys]):
            kwarg = dict(zip(keys, arg))
            case_desc = ', '.join(['%s=%s' % (x, kwarg[x]) for x in keys])
            try:
                res = self._run(**kwarg)
                print('%s: %s' % (case_desc, res))
            except Exception,e:
                print e
                print('%s: failed' % case_desc)
                pass

    def productRun(self, **kwargs):
        '''
        Call self._run with a combination of args from input arguments. For example
            self.productRun(a=[v1, v2], b=[v1, v2, v3])
        will result in
            self._run(a=v1, b=v1)
            self._run(a=v1, b=v2)
            self._run(a=v1, b=v3)
            self._run(a=v2, b=v1)
            self._run(a=v2, b=v2)
            self._run(a=v2, b=v3)
        '''
        keys = kwargs.keys()
        for arg in product(*[kwargs[x] for x in keys]):
            kwarg = dict(zip(keys, arg))
            case_desc = ', '.join(['%s=%s' % (x, kwarg[x]) for x in keys])
            try:
                res = self._run(**kwarg)
                print('%s: %s' % (case_desc, res))
            except Exception,e:
                print e
                print('%s: failed' % case_desc)
                pass

    def run(self):
        print('Please define your own run function')
        return 0

    def _run(self, *args, **kwargs):
        print('Please define your own _run function')
        return 0



class TestRandomMating(PerformanceTest):
    def __init__(self, time = 60):
        PerformanceTest.__init__(self, 'Standard random mating scheme, results are number of generations in %d seconds.' % time)
        self.time = time

    def run(self):
        # overall running case
        self.productRun(size=[1000, 10000, 100000],
            loci=[10, 100, 10000, 100000])
        return True

    def _run(self, size, loci):
        # single test case
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomMating(),
        )
        return gens



class TestRandomMatingWithSelecttion(PerformanceTest):
    def __init__(self, time = 60):
        PerformanceTest.__init__(self, 'Random mating with selection, results are number of generations in %d seconds.' % time)
        self.time = time

    def run(self):
        # overall running case
        self.productRun(size=[1000, 10000, 100000],
            loci=[10, 100, 10000, 100000])
        return True

    def _run(self, size, loci):
        # single test case
        pop = Population(size=size, loci=loci, infoFields='fitness')
        gens = pop.evolve(
            initOps=[
                InitSex(),
                InitGenotype(freq=[0.5, 0.5])
            ],
            preOps=[
                TicToc(output='', stopAfter=self.time),
                MaSelector(loci=0, fitness=[1, 0.999, 0.998]),
            ],
            matingScheme=RandomMating(),
        )
        return gens




if __name__ == '__main__':
    if len(sys.argv) > 1:
        # selected test to run
        tests = []
        for test in sys.argv[1:]:
            tests.extend([func for func in dir() if func.startswith('Test') and test in func])
    else:
        tests = [func for func in dir() if func.startswith('Test')]
    for test in tests:
        testObj = eval(test + '()')
        print testObj.describe()
        testObj.run()
        print 'End of test %s\n' % test
