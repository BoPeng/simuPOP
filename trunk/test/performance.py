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
#     performance.py [-a] [-b|s|l] [-j #] [testname1] [testname2] ...
#
# where
#     testname can be any string within a test. If -b (binary) -s (short) or -l (long)
#     allele type is unspecified this script will be run for all three allele types.
#     -j # specifies number of threads, all threads will be used by default.
#
#     -a analysis mode. list performance test results in performance.csv in a
#        more readable format.
# All the test results will be written to a file 'performance.log'
#
#
import os, sys, time, platform, logging, subprocess, timeit, random, csv

try:
    from itertools import product
except:
    # python 2.4 does not have this function
    def product(*args, **kwds):
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)

if sys.version_info[0] >= 3:
    def callable(obj):
        return hasattr(obj, '__call__')

class PerformanceTest:
    def __init__(self, desc, logger):
        self.description = desc
        self.expected_time = 0
        self.logger = logger
        # machine name, os, CPU, time
        uname = platform.uname()
        self.logger.debug('Host: %s (%s, %s, %s)' % (uname[1], uname[0], uname[4], uname[5]))
        self.logger.debug('Time: %s' % time.asctime())
        info = moduleInfo()
        opt = ''
        if info['optimized']:
            opt = ', optimized'
        self.logger.debug('Python: %s %s' % (info['python'], info['compiler']))
        self.logger.debug('simuPOP: %s (rev %d, module %s%s, %d threads)' % \
            (info['version'], info['revision'], info['alleleType'], opt, info['threads']))
        self.logger.debug('Test: %s' % self.description)
        return desc

    def runTests(self):
        '''Run all tests and record the result to a log file.'''
        results = self.run()
        self.logger.info(', '.join(['%s' % x for x in results]))

    def sequentialRun(self, **kwargs):
        '''
        Call self._run with a list of arguments, for example
            self.sequentialRun(a=[v1, v2], b=[v1, v2])
        will result in
            self._run(a=v1, b=v1)
            self._run(a=v2, b=v2)

        This function should return the return value of each test
        '''
        results = []
        keys = kwargs.keys()
        for arg in zip(*[kwargs[x] for x in keys]):
            kwarg = dict(zip(keys, arg))
            case_desc = ', '.join(['%s=%s' % (x, kwarg[x].__name__ if callable(kwarg[x]) else kwarg[x]) for x in keys])
            try:
                res = self._run(**kwarg)
                self.logger.debug('%s: %s' % (case_desc, res))
                results.append(res)
            except Exception,e:
                print e
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
        keys = list(kwargs.keys())
        # this is to make sure the results are outputted in the same order for all platofmrs
        keys.sort()
        for arg in product(*[kwargs[x] for x in keys]):
            kwarg = dict(zip(keys, arg))
            case_desc = ', '.join(['%s=%s' % (x, kwarg[x].__name__ if callable(kwarg[x]) else kwarg[x]) for x in keys])
            try:
                res = self._run(**kwarg)
                self.logger.debug('%s: %s' % (case_desc, res))
                results.append(res)
            except:
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
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomMating(ops=[
                MendelianGenoTransmitter(),
                # in some cases, mating takes so much time so we have to stop in the middle
                TicToc(output='', stopAfter=self.time)
            ]),
        )
        return gens

class TestRandomSelection(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'Random selection scheme, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomSelection(ops=CloneGenoTransmitter()),
        )
        return gens

class TestPolygamousMating(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'Polygamous mating scheme, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time
 
    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci, infoFields=['father_idx','mother_idx'])
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=PolygamousMating(polySex=MALE, polyNum=2, ops=[MendelianGenoTransmitter(), ParentsTagger()]),
        )
        return gens

class TestPyParentChooser(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'HomoMating scheme with PyParentChooser , results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time
 
    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])


    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        def retInd(pop, subPop):
            while True:
                yield pop.individual(random.randint(0, pop.subPopSize(subPop) - 1)), pop.individual(random.randint(0, pop.subPopSize(subPop) - 1))
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=HomoMating(
                chooser=PyParentsChooser(retInd),
                generator=OffspringGenerator(
                    ops=[CloneGenoTransmitter()],
                    numOffspring=1,
                    sexMode=RANDOM_SEX),
            ),
        )
        return gens


class TestIdTagger(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'Test idTagger, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci, infoFields='ind_id')
        gens = pop.evolve(
            initOps=[
                    InitSex(),
                    IdTagger(),
            ],
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomMating(ops=IdTagger()),
        )
        return gens

class TestPedigreeTagger(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'Test PedigreeTagger, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci, infoFields=['ind_id', 'father_id', 'ped_id'], ancGen=1)
        tagID(pop, reset=True)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomSelection(ops=[ 
                CloneGenoTransmitter(), IdTagger(),
                PedigreeTagger(infoFields='father_id')]),
        )
        return gens

class TestInheritTagger(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'Test InheritTagger, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci, infoFields='x')
        for sp in range(pop.numSubPop()):
            pop.individual(0,sp).x = 1
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomMating(ops=[ 
                MendelianGenoTransmitter(),
                InheritTagger(mode=MAXIMUM, infoFields='x')]),
        )
        return gens



class TestSelfMating(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'Self mating scheme, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=SelfMating(ops=SelfingGenoTransmitter()),
        )
        return gens


class TestHaplodiploidMating(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'Haplodiploid mating scheme, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=HaplodiploidMating(ops=HaplodiploidGenoTransmitter()),
        )
        return gens


class TestMitochondrialGenoTransmitter(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'MitochondrialGenoTransmitter, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomMating(ops=MitochondrialGenoTransmitter()),
        )
        return gens

class TestRecombinator(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'Recombinator, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=RandomMating(ops=Recombinator(rates = 0.4, loci=[1,3,8])),
        )
        return gens

class TestCombinedParentsChooser(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'CombinedParentsChooser, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=HomoMating(
                chooser=CombinedParentsChooser(RandomParentChooser(MALE_ONLY),RandomParentChooser(FEMALE_ONLY)),
                generator=OffspringGenerator(
                    ops=MendelianGenoTransmitter(),
                    numOffspring=2,
                    sexMode=RANDOM_SEX),
                ),
        )
        return gens

class TestHeteroMating(PerformanceTest):
    def __init__(self, logger, time=30):
        PerformanceTest.__init__(self, 'HeteroMating, results are number of generations in %d seconds.' % int(time),
            logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000], loci=[10, 100, 10000])

    def _run(self, size, loci):
        # single test case
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
        pop = Population(size=size, loci=loci, infoFields=['father_idx','mother_idx'])
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=HeteroMating(
                [RandomMating(numOffspring=2, subPops=0, ops=[MendelianGenoTransmitter(), ParentsTagger()]),
                RandomMating(numOffspring=4, subPops=1, ops=[MendelianGenoTransmitter(), ParentsTagger()])]),
        )
        return gens

def createPop(size, loci=100, aff=False, vsp=False):
    pop = Population(size=size, loci=loci)
    initSex(pop)
    initGenotype(pop, freq=[0.5,0.5])
    if aff:
        maPenetrance(pop, loci=0, penetrance=[0.1, 0.3, 0.4])
    if vsp:
        pop.setVirtualSplitter(SexSplitter())
    return pop

class TestStatHaploFreq(PerformanceTest):
    def __init__(self, logger, repeats=200):
        PerformanceTest.__init__(self, 'Stat HaploFreq, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000, [10000]*10], loci=[100, 1000])

    def _run(self, size, loci):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import createPop, stat\n'
                'pop = createPop(size=%s, loci=%s)' % (size, loci),
            stmt = "stat(pop, haploFreq=[[0,1], [2,3], [3,4], [4,5], [5,6]], vars=['haploFreq', 'haploFreq_sp'])\n"
                "pop.vars().clear()")
        return t.timeit(number=self.repeats)

class TestStatHaploHomoFreq(PerformanceTest):
    def __init__(self, logger, repeats=200):
        PerformanceTest.__init__(self, 'Stat HaploHomoFreq, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000, [10000]*10], loci=[100, 1000])

    def _run(self, size, loci):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import createPop, stat\n'
                'pop = createPop(size=%s, loci=%s)' % (size, loci),
            stmt = "stat(pop, haploHomoFreq=[[0,1], [2,3], [3,4], [4,5], [5,6]], vars=['haploHomoFreq', 'haploHomoFreq_sp'])\n"
                "pop.vars().clear()")
        return t.timeit(number=self.repeats)

class TestStatAlleleFreq(PerformanceTest):
    def __init__(self, logger, repeats=30):
        PerformanceTest.__init__(self, 'Stat AlleleFreq, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000, [10000]*10], loci=[100, 1000])

    def _run(self, size, loci):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import createPop, stat, ALL_AVAIL\n'
                'pop = createPop(size=%s, loci=%s)' % (size, loci),
            stmt = "stat(pop, alleleFreq=ALL_AVAIL, vars=['allelFreq', 'alleleFreq_sp'])\n"
                "pop.vars().clear()")
        return t.timeit(number=self.repeats)

class TestStatGenoFreq(PerformanceTest):
    def __init__(self, logger, repeats=500):
        PerformanceTest.__init__(self, 'Stat GenoFreq, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000, [10000]*10], loci=[100, 1000])

    def _run(self, size, loci):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import createPop, stat\n'
                'pop = createPop(size=%s, loci=%s)' % (size, loci),
            stmt = "stat(pop, genoFreq=[0,1,2,3], vars=['genoFreq', 'genoFreq_sp'])\n"
                "pop.vars().clear()")
        return t.timeit(number=self.repeats)

class TestStatHeteroFreq(PerformanceTest):
    def __init__(self, logger, repeats=20):
        PerformanceTest.__init__(self, 'Stat HeteroFreq, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.productRun(size=[10000, 100000, [10000]*10], loci=[100, 1000])

    def _run(self, size, loci):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import createPop, stat, ALL_AVAIL\n'
                'pop = createPop(size=%s, loci=%s)' % (size, loci),
            stmt = "stat(pop, heteroFreq=ALL_AVAIL, vars=['heteroFreq', 'heteroFreq_sp'])\n"
                "pop.vars().clear()")
        return t.timeit(number=self.repeats)

class TestStatNeutrality(PerformanceTest):
    def __init__(self, logger, repeats=10):
        PerformanceTest.__init__(self, 'Stat Neutrality, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.productRun(size=[10000, [1000]*10], loci=[1,10,100])

    def _run(self, size, loci):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import createPop, stat, ALL_AVAIL\n'
                'pop = createPop(size=%s, loci=%s)' % (size, loci),
            stmt = "stat(pop, neutrality=ALL_AVAIL, vars=['neutrality', 'neutrality_sp'])\n"
                "pop.vars().clear()")
        return t.timeit(number=self.repeats)

class TestStatNumOfMales(PerformanceTest):
    def __init__(self, logger, repeats=5000):
        PerformanceTest.__init__(self, 'Stat NumOfMales, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(size=[100000, [10000]*10])

    def _run(self, size):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import createPop, stat, ALL_AVAIL\n'
                'pop = createPop(size=%s,vsp=True)' % (size),
            stmt = "stat(pop, numOfMales=True, subPops=[(ALL_AVAIL,ALL_AVAIL), (ALL_AVAIL,0), (ALL_AVAIL, 1)], vars=['numOfMales', 'numOfMales_sp'])\n"
                "pop.vars().clear()")
        return t.timeit(number=self.repeats)

class TestStatNumOfAffected(PerformanceTest):
    def __init__(self, logger, repeats=5000):
        PerformanceTest.__init__(self, 'Stat NumOfAffected, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(size=[100000, [10000]*10])

    def _run(self, size):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import createPop, stat, ALL_AVAIL\n'
                'pop = createPop(size=%s, aff=True, vsp=True)' % size,
            stmt = "stat(pop, numOfAffected=True, subPops=[(ALL_AVAIL,ALL_AVAIL), (ALL_AVAIL,0), (ALL_AVAIL, 1)], vars=['numOfAffected', 'numOfAffected_sp'])\n"
                "pop.vars().clear()")
        return t.timeit(number=self.repeats)

class TestInitGenotypeCase1(PerformanceTest):
    def __init__(self, logger, repeats=10):
        PerformanceTest.__init__(self, 'InitGenotype case 1, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(size=[100000, [10000]*10])

    def _run(self, size):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import Population, initGenotype\n'
                'pop = Population(size=%s,loci=1000)' % (size),
            stmt = 'initGenotype(pop,genotype=[1,2,3,4,5,6])') 
        return t.timeit(number=self.repeats)

class TestInitGenotypeCase2(PerformanceTest):
    def __init__(self, logger, repeats=10):
        PerformanceTest.__init__(self, 'InitGenotype case 2, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(size=[100000, [10000]*10])

    def _run(self, size):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import Population, initGenotype\n'
                'pop = Population(size=%s,loci=1000)' % (size),
            stmt = 'initGenotype(pop,prop=[.2,.8])') 
        return t.timeit(number=self.repeats)

class TestInitGenotypeCase3(PerformanceTest):
    def __init__(self, logger, repeats=10):
        PerformanceTest.__init__(self, 'InitGenotype case 3, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(size=[100000, [10000]*10])

    def _run(self, size):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import Population, initGenotype\n'
                'pop = Population(size=%s,loci=1000)' % (size),
            stmt = 'initGenotype(pop,haplotypes=[[0,0,1,1],[1,1,0,0]], prop=[.2,.8])')
        return t.timeit(number=self.repeats)

class TestInitGenotypeCase4(PerformanceTest):
    def __init__(self, logger, repeats=10):
        PerformanceTest.__init__(self, 'InitGenotype case 4, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(size=[100000, [10000]*10])

    def _run(self, size):
        # single test case
        t = timeit.Timer(
            setup = 'from __main__ import Population, initGenotype\n'
                'pop = Population(size=%s,loci=1000)' % (size),
            stmt = 'initGenotype(pop,freq=[0.5,0.5])')
        return t.timeit(number=self.repeats)

class TestInitSex(PerformanceTest):
    def __init__(self, logger, repeats=1000):
        PerformanceTest.__init__(self, 'InitSex, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(size=[1000000, [100000]*10])

    def _run(self, size):
        # single test case
        pop = Population(size = size)
        t = timeit.Timer(
            setup = 'from __main__ import Population,initSex, MALE,FEMALE\n'
                'pop = Population(size=%s,loci=1000)' % (size),
            stmt = 'initSex(pop)')
        return t.timeit(number=self.repeats)
    
class TestInitInfo(PerformanceTest):
    
    def __init__(self, logger, repeats=1000):
        PerformanceTest.__init__(self, 'InitInfo, results are time (not processor time) to apply operator for %d times.' % int(repeats),
            logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(size=[1000000, [100000]*10])

    def _run(self, size):
        # single test case
        pop = Population(size = size)
        t = timeit.Timer(
               setup = 'from __main__ import Population, initInfo\n'
        "pop = Population(size=%s,loci=100, infoFields=['a','b'])" % (size),
        stmt = "initInfo(pop,[1,2,3,4,5,6,7],infoFields=['a','b'])")
        return t.timeit(number=self.repeats)


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
        if size * loci * moduleInfo()['alleleBits'] / 8 > 1e9:
            return 0
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
            matingScheme=RandomMating(ops=[
                MendelianGenoTransmitter(),
                # in some cases, mating takes so much time so we have to stop in the middle
                TicToc(output='', stopAfter=self.time)
            ]),
        )
        return gens


class TestPyOperatorFuncCall(PerformanceTest):
    def __init__(self, logger, time=60):
        PerformanceTest.__init__(self, 'Test the performance of function call for operator PyOperator.'
            'The results are number of generations in %d seconds' % int(time), logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[1000, 10000, 100000],
                func=[self._popFunc, self._popFuncIterSex, self._indFunc, self._indFuncCallSex])

    def _popFunc(self, pop):
        return True

    def _popFuncIterSex(self, pop):
        for ind in pop.individuals():
            ind.sex()
        return True

    def _indFunc(self, ind):
        return True

    def _indFuncCallSex(self, ind):
        ind.sex()
        return True

    def _run(self, size, func):
        # single test case
        pop = Population(size=size, loci=1)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=[
                TicToc(output='', stopAfter=self.time),
                PyOperator(func=func)
            ],
            matingScheme=CloneMating(ops=[]),
        )
        return gens

class TestDuringMatingPyOperator(PerformanceTest):
    def __init__(self, logger, time=60):
        PerformanceTest.__init__(self, 'Test the performance of during-mating function call for operator PyOperator.'
            'The results are number of generations in %d seconds' % int(time), logger)
        self.time = time

    def run(self):
        # overall running case
        return self.productRun(size=[1000, 10000, 100000],
                func=[self._func1Param, self._func1ParamCallSex, self._func4Param])

    def _func1Param(self, ind):
        return True

    def _func1ParamCallSex(self, ind):
        ind.sex()
        return True

    def _func4Param(self, ind, pop, dad, mom):
        return True

    def _run(self, size, func):
        # single test case
        pop = Population(size=size, loci=1)
        gens = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(output='', stopAfter=self.time),
            matingScheme=CloneMating(ops=[
                PyOperator(func=func)]),
        )
        return gens


def getIteratablePop(vsp):
    pop = Population(size=100000, loci=2, infoFields=['x', 'y'])
    initSex(pop)
    initInfo(pop, lambda:random.randint(0, 10), infoFields=['x', 'y'])
    initGenotype(pop, freq=[0.5, 0.5])
    maPenetrance(pop, loci=0, penetrance=[0.1, 0.2, 0.3])
    if vsp == 'all ind':
        pass
    elif vsp == 'all males':
        pop.setVirtualSplitter(SexSplitter())
    elif vsp == 'all unaffected':
        pop.setVirtualSplitter(AffectionSplitter())
    elif vsp == 'first 50%':
        pop.setVirtualSplitter(ProportionSplitter([0.5, 0.5]))
    elif vsp == '20000 to 50000':
        pop.setVirtualSplitter(RangeSplitter([20000, 50000]))
    elif vsp == 'x < 5':
        pop.setVirtualSplitter(InfoSplitter(field='x', cutoff=[5]))
    elif vsp == 'y < 2':
        pop.setVirtualSplitter(InfoSplitter(field='y', cutoff=[2]))
    elif vsp == 'y == 5':
        pop.setVirtualSplitter(InfoSplitter(field='y', values=[5]))
    elif vsp == 'x in [3, 5]':
        pop.setVirtualSplitter(CombinedSplitter([
            InfoSplitter(field='x', values=[3, 5])], vspMap=[(0,1)]))
    elif vsp == 'x < 5 and y < 5':
        pop.setVirtualSplitter(ProductSplitter([
            InfoSplitter(field='x', cutoff=[5]),
            InfoSplitter(field='y', cutoff=[5])]))
    elif vsp == 'g1 == [0,0]':
        pop.setVirtualSplitter(GenotypeSplitter(loci=0, alleles=[0,0]))
    elif vsp == 'g1 == [0,0], g2 == [0,0]':
        pop.setVirtualSplitter(GenotypeSplitter(loci=[0,1], alleles=[0,0,0,0]))
    else:
        raise ValueError('Incorrect virtual subpopulation name')
    return pop


class TestIteratingVSPs(PerformanceTest):
    def __init__(self, logger, repeats=500):
        PerformanceTest.__init__(self, 'Test the performance of iterating through virtual subpopulations. '
            'The results are time of iterating for %d times' % repeats, logger)
        self.repeats = repeats

    def run(self):
        # overall running case
        return self.sequentialRun(vsp=[
            'all ind',
            'all males',
            'all unaffected',
            'first 50%',
            '20000 to 50000',
            'x < 5',
            'y < 2',
            'y == 5',
            'x in [3, 5]',
            'x < 5 and y < 5',
            'g1 == [0,0]',
            'g1 == [0,0], g2 == [0,0]',
        ])


    def _run(self, vsp):
        if vsp == 'all ind':
            s = 'for ind in pop.individuals():\n    pass\n'
        else:
            s = 'for ind in pop.individuals([0,0]):\n    pass\n'
        t = timeit.Timer(stmt=s, setup='from __main__ import getIteratablePop\n'
            'pop = getIteratablePop("%s")' % vsp)
        return t.timeit(number=self.repeats)


def analyze(test):
    '''Output performance statistics for a test
    '''
    # read csv file
    reader = csv.reader(open('performance.csv', 'r'), delimiter=',',
        skipinitialspace=True) 
    # record for specified test
    # name, date, sec, machine, platform-threads, python, ver, rev, type, rec...
    uname = platform.uname()
    records = [rec for rec in reader if rec[0] == test]
    # get different platforms and number of threads
    pfs = set([(rec[3], rec[4]) for rec in records if rec[3] != uname[1]])
    pfs = list(pfs)
    pfs.sort()
    # put results for local machine the last so that we can locate them easily
    pfs_local = set([(rec[3], rec[4]) for rec in records if rec[3] == uname[1]])
    pfs_local = list(pfs_local)
    pfs_local.sort()
    #
    pfs.extend(pfs_local)
    for pf in pfs:
        pfRecords = [rec for rec in records if (rec[3], rec[4]) == pf]
        revs = list(set([rec[7] for rec in pfRecords]))
        revs.sort()
        print test, pf[0], pf[1], revs
        numStat = len(pfRecords) - 9
        for stat in range(9, len(pfRecords[0])):
            print 'STAT %2d' % (stat - 8),
            # each stat, for different revisions, and type
            for type in ['short', 'long', 'binary']:
                print '|' if len(revs) < 4 else '\n %3s|' % type[0],
                for rev in revs:
                    one = [rec for rec in pfRecords if rec[8] == type and rec[7] == rev]
                    if len(one) == 0:
                        print '     ??',
                    elif len(one) == 1:
                        if '.' in one[0][stat]:
                            print '%7.2f' % float(one[0][stat]),
                        else:
                            print '%7s' % (one[0][stat]),
                    else:
                        if '.' in one[0][stat]:
                            print '%17s' % ('(%7.2f - %7.2f)' % (min([float(x[stat]) for x in one]), max([float(x[stat]) for x in one]))),
                        else:
                            print '%17s' % ('(%d - %d)' % (min([int(x[stat]) for x in one]), max([int(x[stat]) for x in one]))),
            print
        print
      
if __name__ == '__main__':
    # 
    # Figure out tests to run
    #
    if sum([not x.startswith('-') for x in sys.argv[1:]]) > 0:
        # selected test to run
        tests = []
        for test in sys.argv[1:]:
            tests.extend([func for func in dir() if func.startswith('Test') and test in func])
    else:
        tests = [func for func in dir() if func.startswith('Test')]
    #
    # if in analysis mode
    if '-a' in sys.argv:
        for test in tests:
            analyze(test)
        sys.exit(0)
    # 
    # Perform tests
    #
    alleleType = 'all'
    # allele type can be specified by --alleleType=long/short/binary
    if '-b' in sys.argv:
        alleleType = 'binary'
        sys.argv.remove('-b')
    elif '-s' in sys.argv:
        alleleType = 'short'
        sys.argv.remove('-s')
    elif '-l' in sys.argv:
        alleleType = 'long'
        sys.argv.remove('-l')
    #
    if alleleType == 'all':
        for t in ['s', 'l', 'b']:
            ret = subprocess.call([sys.executable, sys.argv[0], '-%s' % t] + sys.argv[1:])
            if ret != 0:  # if crash or killed
                print 'Error: A non-zero return value is returned for module %s' % t
                sys.exit(ret)
        sys.exit(0)
    #    
    numThreads = None
    if '-j' in sys.argv:
        idx = sys.argv.index('-j')
        numThreads = int(sys.argv[idx + 1])
        sys.argv.pop(idx+1)
        sys.argv.pop(idx)
    elif True in [x.startswith('-j') for x in sys.argv]:
        idx = [x.startswith('-j') for x in sys.argv].index(True)
        numThreads = int(sys.argv[idx][2:])
        sys.argv.pop(idx)
    #
    import simuOpt
    simuOpt.setOptions(alleleType=alleleType, quiet=True, optimized=True, numThreads=numThreads)
    from simuPOP import *
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
        compiler = 'gcc'
        if info['compiler'].find("Intel") > -1: 
            compiler = 'intel'
        elif info['compiler'].find("MSC") > -1:
            compiler = 'vc'
        summaryFile.setFormatter(logging.Formatter('%%(name)s, %%(asctime)s, %s, %s-%dthreads, python%s, simuPOP-%s, rev%d-%s, %s, %%(message)s' % \
            (uname[1], uname[4], info['threads'], info['python'], info['version'], info['revision'], compiler, info['alleleType'])))
        #
        logger.addHandler(logFile)
        logger.addHandler(summaryFile)
        #
        testObj = eval(test + '(logger)')
        testObj.runTests()
        logger.debug('End of test %s' % test)
