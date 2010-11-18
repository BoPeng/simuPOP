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
#     testname can be any string within a test. If --alleleType is unspecified
#     this script will be run for all three allele types.
#
# All the test results will be written to a file 'performance.log'
#
#
import os, sys, time, platform, logging, subprocess
from itertools import product

alleleType = 'all'
# allele type can be specified by --alleleType=long/short/binary
if True in [x.startswith('--alleleType=') for x in sys.argv]:
    idx = [x.startswith('--alleleType=') for x in sys.argv].index(True)
    alleleType = sys.argv[idx][13:]
    if not alleleType in ['short', 'long', 'binary']:
        raise ValueError('Incorrect allele type. Please use --alleleType=long, short or binary')
    sys.argv.pop(idx)

if alleleType == 'all':
    for t in ['short', 'long', 'binary']:
        ret = subprocess.call(['python', sys.argv[0], '--alleleType=%s' % t] + sys.argv[1:])
        if ret != 0:  # if crash or killed
            sys.exit(ret)
    sys.exit(0)
    
import simuOpt
simuOpt.setOptions(alleleType=alleleType, quiet=True, optimized=True)
from simuPOP import *

class PerformanceTest:
    def __init__(self, desc, logger):
        self.description = desc
        self.expected_time = 0
        self.logger = logger
        # machine name, os, CPU, time
        uname = platform.uname()
        self.logger.debug('Host: %s (%s, %s)' % (uname[1], uname[0], uname[5]))
        self.logger.debug('Time: %s' % time.asctime())
        info = moduleInfo()
        opt = ''
        if info['optimized']:
            opt = ', optimized'
        self.logger.debug('Python: %s %s' % (info['python'], info['compiler']))
        self.logger.debug('simuPOP: %s (rev %d, module %s%s)' % \
            (info['version'], info['revision'], info['alleleType'], opt))
        self.logger.debug('Test: %s' % self.description)
        return desc

    def runTests(self):
        '''Run all tests and record the result to a log file.'''
        results = self.run()
        self.logger.info(', '.join(['%s' % x for x in results]))

    def serialRun(self, **kwargs):
        '''
        Call self._run with a list of arguments, for example
            self.serialRun(a=[v1, v2], b=[v1, v2])
        will result in
            self._run(a=v1, b=v1)
            self._run(a=v2, b=v2)

        This function should return the return value of each test
        '''
        results = []
        keys = kwargs.keys()
        for arg in zip(*[kwargs[x] for x in keys]):
            kwarg = dict(zip(keys, arg))
            case_desc = ', '.join(['%s=%s' % (x, kwarg[x]) for x in keys])
            try:
                res = self._run(**kwarg)
                self.logger.debug('%s: %s' % (case_desc, res))
                results.append(res)
            except Exception,e:
                self.logger.debug('%s: failed' % case_desc)
                results.append(0)
                pass
        return results

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
        results = []
        keys = kwargs.keys()
        # this is to make sure the results are outputted in the same order for all platofmrs
        keys.sort()
        for arg in product(*[kwargs[x] for x in keys]):
            kwarg = dict(zip(keys, arg))
            case_desc = ', '.join(['%s=%s' % (x, kwarg[x]) for x in keys])
            try:
                res = self._run(**kwarg)
                self.logger.debug('%s: %s' % (case_desc, res))
                results.append(res)
            except Exception,e:
                self.logger.debug('%s: failed' % case_desc)
                results.append(0)
                pass
        return results

    def run(self):
        '''Run tests and return results'''
        print('Please define your own run function')
        return []

    def _run(self, *args, **kwargs):
        print('Please define your own _run function')
        return 0



class TestBasicRandomMating(PerformanceTest):
    def __init__(self, logger, time=60):
        PerformanceTest.__init__(self, 'Standard random mating scheme, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[1000, 10000, 100000], loci=[10, 100, 10000, 100000])

    def _run(self, size, loci):
        # single test case
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomMating(),
        )
        return gens


class TestRandomMatingWithSelection(PerformanceTest):
    def __init__(self, logger, time=60):
        PerformanceTest.__init__(self, 'Random mating with selection, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[1000, 10000, 100000], loci=[10, 100, 10000, 100000])

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
    # 
    # Figure out tests to run
    #
    if len(sys.argv) > 1:
        # selected test to run
        tests = []
        for test in sys.argv[1:]:
            tests.extend([func for func in dir() if func.startswith('Test') and test in func])
    else:
        tests = [func for func in dir() if func.startswith('Test')]
    #
    #
    logging.basicConfig(level=logging.DEBUG, format='%(name)s: %(message)s')
    for test in tests:
        logger = logging.getLogger(test)
        # A log file with debug information
        logFile = logging.FileHandler('performance.log')
        logFile.setLevel(logging.DEBUG)
        logFile.setFormatter(logging.Formatter('%(name)s : %(message)s'))
        # A summary file with overall results
        summaryFile = logging.FileHandler('performance.csv')
        summaryFile.setLevel(logging.INFO)
        uname = platform.uname()
        info = moduleInfo()
        summaryFile.setFormatter(logging.Formatter('%%(name)s, %%(asctime)s, %s, %s, python%s, simuPOP-%s, rev%d, %s, %%(message)s' % \
            (uname[1], uname[4], info['python'], info['version'], info['revision'], info['alleleType'])))
        #
        logger.addHandler(logFile)
        logger.addHandler(summaryFile)
        #
        testObj = eval(test + '(logger)')
        testObj.runTests()
        logger.debug('End of test %s' % test)
