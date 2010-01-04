#!/usr/bin/env python
#
# Purpose:
#     Testing selection.
#
# Author:
#     Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

class TestSelector(unittest.TestCase):

    def testMapSelector(self):
        'Testing a map selector'
        pop = Population(size=1000, loci=[1], infoFields=['a', 'fitness', 'b'])
        initGenotype(pop, freq=[.2, .8])
        MapSelector(loci=[0], fitness={(0,0):1, (0,1):0.5, (1,1):0.25}).apply(pop)
        for ind in pop.individuals():
            if ind.genotype() == (0,0):
                self.assertEqual(ind.fitness, 1.0)
            elif ind.genotype() == (0,1):
                self.assertEqual(ind.fitness, 0.5)
            elif ind.genotype() == (1,0):
                self.assertEqual(ind.fitness, 0.5)
            elif ind.genotype() == (1,1):
                self.assertEqual(ind.fitness, 0.25)
        #
        if moduleInfo()['alleleType'] == 'binary':
            initGenotype(pop, freq=[.2, .8])
            MapSelector(loci=[0], fitness={(0,0):1, (0,1):0.5, (1,0):0.3, (1,1):0.25}).apply(pop)
            for ind in pop.individuals():
                if ind.genotype() == (0,0):
                    self.assertEqual(ind.fitness, 1.0)
                elif ind.genotype() == (0,1):
                    self.assertEqual(ind.fitness, 0.5)
                elif ind.genotype() == (1,0):
                    self.assertEqual(ind.fitness, 0.3)
                elif ind.genotype() == (1,1):
                    self.assertEqual(ind.fitness, 0.25)
            #
            pop = Population(size=1000, loci=[1], infoFields=['a', 'fitness', 'b'],
                ploidy=1)
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0], fitness={(0,):1, (1,):0.9}).apply(pop)
            for ind in pop.individuals():
                if ind.genotype() == (0,):
                    self.assertEqual(ind.fitness, 1.0)
                elif ind.genotype() == (2,):
                    self.assertEqual(ind.fitness, 0.9)
            # haplodiploid population
            pop = Population(size=1000, loci=[1], infoFields=['a', 'fitness', 'b'],
                ploidy=HAPLODIPLOID)
            initSex(pop)
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0], fitness={(0,0):1, (0,1):0.5, (1,0):0.3, (1,1):0.25,
                (0,):0.8, (1,):0.9}).apply(pop)
            for ind in pop.individuals():
                if ind.genotype() == (0,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif ind.genotype() == (0,1) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.5)
                elif ind.genotype() == (1,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.3)
                elif ind.genotype() == (1,1) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.25)
                elif ind.genotype() == (0,0) and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.8)
                elif ind.genotype() == (0,1) and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.8)
                elif ind.genotype()[0] == 1 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.9)
            #
            # sex chromosome
            pop = Population(size=1000, loci=[1, 2], infoFields=['a', 'fitness', 'b'],
                chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
            initSex(pop)
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0], fitness={(0,):0.9, (1,):0.8, (0,0):1, (0,1):0.5, (1,0):0.3, (1,1):0.25}).apply(pop)
            for ind in pop.individuals():
                if (ind.allele(0,0),ind.allele(0,1)) == (0,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif (ind.allele(0,0),ind.allele(0,1)) == (0,1) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.5)
                elif (ind.allele(0,0),ind.allele(0,1)) == (1,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.3)
                elif (ind.allele(0,0),ind.allele(0,1)) == (1,1) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.25)
                elif ind.allele(0,0) == 0 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.9)
                elif ind.allele(0,0) == 2 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.8)
            # 
            MapSelector(loci=[2], fitness={(0,):0.9, (1,):0.8, ():1}).apply(pop)
            for ind in pop.individuals():
                if ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif ind.allele(2,1) == 0 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.9)
                elif ind.allele(2,1) == 2 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.8)
            # multiple loci
            pop = Population(size=1000, loci=[1, 2], infoFields=['a', 'fitness', 'b'],
                chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
            initSex(pop)
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0,1], fitness={(0,0):1, (0,1):0.5, (1,0):0.3, (1,1):0.25}).apply(pop)
            for ind in pop.individuals():
                if (ind.allele(0,0),ind.allele(0,1)) == (0,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif (ind.allele(0,0),ind.allele(0,1)) == (0,1) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.5)
                elif (ind.allele(0,0),ind.allele(0,1)) == (1,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.3)
                elif (ind.allele(0,0),ind.allele(0,1)) == (1,1) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.25)
                elif (ind.allele(0,0),ind.allele(1,1)) == (0,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif (ind.allele(0,0),ind.allele(1,1)) == (0,1) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.5)
                elif (ind.allele(0,0),ind.allele(1,1)) == (1,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.3)
                elif (ind.allele(0,0),ind.allele(1,1)) == (1,1) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.25)
        else:
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0], fitness={(0,0):1, (0,2):0.5, (2,0):0.3, (2,2):0.25}).apply(pop)
            for ind in pop.individuals():
                if ind.genotype() == (0,0):
                    self.assertEqual(ind.fitness, 1.0)
                elif ind.genotype() == (0,2):
                    self.assertEqual(ind.fitness, 0.5)
                elif ind.genotype() == (2,0):
                    self.assertEqual(ind.fitness, 0.3)
                elif ind.genotype() == (2,2):
                    self.assertEqual(ind.fitness, 0.25)
            #
            pop = Population(size=1000, loci=[1], infoFields=['a', 'fitness', 'b'],
                ploidy=1)
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0], fitness={(0,):1, (2,):0.9}).apply(pop)
            for ind in pop.individuals():
                if ind.genotype() == (0,):
                    self.assertEqual(ind.fitness, 1.0)
                elif ind.genotype() == (2,):
                    self.assertEqual(ind.fitness, 0.9)
            # haplodiploid population
            pop = Population(size=1000, loci=[1], infoFields=['a', 'fitness', 'b'],
                ploidy=HAPLODIPLOID)
            initSex(pop)
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0], fitness={(0,0):1, (0,2):0.5, (2,0):0.3, (2,2):0.25,
                (0,):0.8, (2,):0.9}).apply(pop)
            for ind in pop.individuals():
                if ind.genotype() == (0,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif ind.genotype() == (0,2) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.5)
                elif ind.genotype() == (2,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.3)
                elif ind.genotype() == (2,2) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.25)
                elif ind.genotype() == (0,0) and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.8)
                elif ind.genotype() == (0,2) and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.8)
                elif ind.genotype()[0] == 1 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.9)
            #
            # sex chromosome
            pop = Population(size=1000, loci=[1, 2], infoFields=['a', 'fitness', 'b'],
                chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
            initSex(pop)
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0], fitness={(0,):0.9, (2,):0.8, (0,0):1, (0,2):0.5, (2,0):0.3, (2,2):0.25}).apply(pop)
            for ind in pop.individuals():
                if (ind.allele(0,0),ind.allele(0,1)) == (0,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif (ind.allele(0,0),ind.allele(0,1)) == (0,2) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.5)
                elif (ind.allele(0,0),ind.allele(0,1)) == (2,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.3)
                elif (ind.allele(0,0),ind.allele(0,1)) == (2,2) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.25)
                elif ind.allele(0,0) == 0 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.9)
                elif ind.allele(0,0) == 2 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.8)
            # 
            MapSelector(loci=[2], fitness={(0,):0.9, (2,):0.8, ():1}).apply(pop)
            for ind in pop.individuals():
                if ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif ind.allele(2,1) == 0 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.9)
                elif ind.allele(2,1) == 2 and ind.sex() == MALE:
                    self.assertEqual(ind.fitness, 0.8)
            # multiple loci
            pop = Population(size=1000, loci=[1, 2], infoFields=['a', 'fitness', 'b'],
                chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
            initSex(pop)
            initGenotype(pop, freq=[.2, 0, .8])
            MapSelector(loci=[0,1], fitness={(0,0):1, (0,2):0.5, (2,0):0.3, (2,2):0.25}).apply(pop)
            for ind in pop.individuals():
                if (ind.allele(0,0),ind.allele(0,1)) == (0,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif (ind.allele(0,0),ind.allele(0,1)) == (0,2) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.5)
                elif (ind.allele(0,0),ind.allele(0,1)) == (2,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.3)
                elif (ind.allele(0,0),ind.allele(0,1)) == (2,2) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.25)
                elif (ind.allele(0,0),ind.allele(1,1)) == (0,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 1.0)
                elif (ind.allele(0,0),ind.allele(1,1)) == (0,2) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.5)
                elif (ind.allele(0,0),ind.allele(1,1)) == (2,0) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.3)
                elif (ind.allele(0,0),ind.allele(1,1)) == (2,2) and ind.sex() == FEMALE:
                    self.assertEqual(ind.fitness, 0.25)


    def testMaSelector(self):
        'Testing multi-allele selector'
        pop = Population(size=10, loci=[1], infoFields=['a', 'fitness', 'b'])
        initGenotype(pop, freq=[.2, .8])
        MaSelector(loci=[0], fitness=[1, 0.5, 0.25], wildtype=[0]).apply(pop)
        for ind in pop.individuals():
            if ind.genotype() == (0,0):
                assert ind.info('fitness') == 1
            elif ind.genotype() == (0,1):
                assert ind.info('fitness') == 0.5
            elif ind.genotype() == (1,0):
                assert ind.info('fitness') == 0.5
            elif ind.genotype() == (1,1):
                assert ind.info('fitness') == 0.25
        #
        # selector on a population with selection on is not allowed
        # explicitly walk around this.
        initGenotype(pop, freq=[.2, 0, .3, .4, .1])
        # other than 1 alleles
        MaSelector(loci=[0], fitness=[1, 0.5, 0.25], wildtype=[0]).apply(pop)
        for ind in pop.individuals():
            #print ind.genotype(), ind.info('fitness')
            if 0 in ind.genotype():
                assert ind.info('fitness') > 0.25
        # selector on a population with selection on is not allowed
        # explicitly walk around this.
        initGenotype(pop, freq=[.2, 0, .3, .4, .1])
        # more than one wild type
        MaSelector(loci=[0], fitness=[1, 0.5, 0.25], wildtype=[0, 2]).apply(pop)
        for ind in pop.individuals():
            # print ind.genotype(), ind.info('fitness')
            if 0 in ind.genotype() or 2 in ind.genotype():
                assert ind.info('fitness') > 0.25

    def TestAgeOfDistinction(self):
        '''Testing selection in a long time, this is a long test, and
        will not be performed automatically '''
        import rpy
        N = 5000
        p = 0.2
        m = 1000
        s = 0.10
        if os.path.isfile('dist.py'):
            res = {}
            execfile('dist.py', res, res)
            lost = res['lost']
        else:
            fixed = 0
            length_lost = []
            length_fixed = []
            for i in range(m):
                traj = FreqTrajectoryStoch(curGen=200000, fitness=[1, 1+s/2, 1+s], \
                    freq=p, N=N, restartIfFail=True, allowFixation=True)
                if traj[0] > 0.5:
                    fixed += 1
                    length_fixed.append(len(traj))
                else:
                    length_lost.append(len(traj))
                print len(traj)
            out = open('dist.py', 'w')
            fixed = rpy.r.quantile(length_fixed, [0.05, 0.25, 0.5, 0.75, 0.95])
            lost = rpy.r.quantile(length_lost,  [0.05, 0.25, 0.5, 0.75, 0.95])
            if len(length_fixed) > 0:
                print >> out, 'fixed = %s' % str(fixed)
            if len(length_lost) > 0:
                print >> out, 'lost = %s' % str(lost)
            out.close()
            print open('dist.py').read()
        # real case
        simulated = []
        sel = MlSelector(
                        [MaSelector(loci=[0], wildtype=[0], fitness=[1, 1-s/2, 1-s])],
                        mode=HETEROGENEITY
                        )
        # sel = MaSelector(loci=[0], wildtype=[0], fitness=[1, 1-s/2, 1-s])]
        for i in range(100):
            pop = Population(size=[200,N], loci=[1], infoFields=['fitness'])
            simu = Simulator(pop)
            simu.evolve(
                initOps = [InitSex(), InitGenotype(freq=[1-p]+[p/10.]*10)],
                matingScheme = RandomMating(),
                preOps = sel,
                postOps = [
                    Stat(alleleFreq=[0]),
                    TerminateIf('subPop[1]["alleleNum"][0][0] == %d*2' % N),
                    #PyEval(r'"%d\n"%alleleNum[0][0]', step=100)
                ]
            )
            print i, simu.dvars(0).gen
            simulated.append(simu.dvars(0).gen)
            #if simu.dvars(0).gen < lost[0] or simu.dvars(0).gen > lost[1]:
            #    print "Warning: something may be wrong %d outside: [%f %f]. " % (simu.dvars(0).gen, lost[0], lost[1])
        print rpy.r.quantile(simulated,  [0.05, 0.25, 0.5, 0.75, 0.95])


    def testSubPopDirSelection(self):
        'Testing directional selection in s subpopulation using a map selector'
        simu = Simulator(
            Population(size=[200, 1000], ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        # 1. directional selection
        #     w11 > w12 > w22
        #    p -> 1
        simu.evolve(
            initOps=[ InitSex(), InitGenotype(freq=[.5, .5])],
                matingScheme = RandomMating(),
            preOps = MapSelector(loci=0,
                    fitness={(0,0):1, (0,1):0.9, (1,1):.8}),
            postOps = [
                Stat( alleleFreq=[0], genoFreq=[0], vars='alleleFreq_sp'),
                TerminateIf('subPop[1]["alleleFreq"][0][0] < 0.4'),
                TerminateIf('subPop[1]["alleleFreq"][0][0] < 0.8', begin=50)
            ],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)


    def testMapSelectorDirSelection(self):
        'Testing directional selection using a map selector'
        # specify relative fitness: w11, w12/w21, w22
        # NOTE: use spare here to make sure that slector will work
        # when more than one information fields are available.
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        # 1. directional selection
        #     w11 > w12 > w22
        #    p -> 1
        simu.evolve(
            initOps =[ InitSex(), InitGenotype(freq=[.5,.5])],
                matingScheme = RandomMating(),
            preOps = MapSelector(loci=0,
                    fitness={(0,0):1, (0,1):0.9, (1,1):.8}),
            postOps = [
                Stat( alleleFreq=[0], genoFreq=[0]),
                TerminateIf('alleleFreq[0][0] < 0.4'),
                TerminateIf('alleleFreq[0][0] < 0.8', begin=50)
            ],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)


    def testMaSelectorDirSelection(self):
        'Testing directional selection using a multi-allele selector'
        # specify relative fitness: w11, w12/w21, w22
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        # 1. directional selection
        #     w11 > w12 > w22
        #    p -> 1
        simu.evolve(
            initOps=[ InitSex(), InitGenotype(freq=[.5,.5])],
                matingScheme = RandomMating(),
            preOps = MaSelector(loci=0, wildtype=[0],
                    fitness=[1, 0.9, .8]),
            postOps = [
                Stat( alleleFreq=[0], genoFreq=[0]),
                TerminateIf('alleleFreq[0][0] < 0.4'),
                TerminateIf('alleleFreq[0][0] < 0.8', begin=50)
            ],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)

    def testMapSelectorHeteroAdv(self):
        'Testing heterozygous advantage using map selector'
        # specify relative fitness: w11, w12/w21, w22
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        s1 = .1
        s2 = .2
        p = .2/ (.1+.2)
        # 2. heterozygote superiority
        #     w11 < w12, w12 > w22
        #    stable.
        # let
        #        s1 = w12 -    w11
        #        s2 = w12 - w22
        #    p_ = s2/ (s1+s2)
        simu.evolve(
            preOps = MapSelector(loci=0,
                    fitness={(0,0):1-s1, (0,1):1, (1,1):1-s2}),
                matingScheme = RandomMating(),
            postOps = [
                Stat( alleleFreq=[0], genoFreq=[0]),
                TerminateIf('alleleFreq[0][0] < 0.5', begin=50),
                TerminateIf('alleleFreq[0][0] > 0.9', begin=50)
            ],
            initOps=[ InitSex(), InitGenotype(freq=[.5,.5])],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)

    def testMaSelectorHeteroAdv(self):
        'Testing heterozygous advantage using a multi-allele selector'
        s1 = .1
        s2 = .2
        p = .2/ (.1+.2)
        # specify relative fitness: w11, w12/w21, w22
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        # 2. heterozygote superiority
        #     w11 < w12, w12 > w22
        #    stable.
        # let
        #        s1 = w12-    w11
        #        s2 = w12 - w22
        #    p_ = s2/ (s1+s2)
        #
        simu.evolve(
            preOps = MaSelector(loci=0, wildtype=0,
                    fitness=[1-s1, 1, 1-s2]),
                matingScheme = RandomMating(),
            postOps = [
                Stat( alleleFreq=[0], genoFreq=[0]),
                TerminateIf('alleleFreq[0][0] < 0.5', begin=50),
                TerminateIf('alleleFreq[0][0] > 0.9', begin=50)
            ],
            initOps=[ InitSex(), InitGenotype(freq=[.5,.5])],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)

    def testMapSelectorHeteroDisadv(self):
        'Testing heterozygous disadvantage using map selector'
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        # 2. heterozygote inferiority
        #     w11 > w12, w12 < w22
        #    p unstable, become fix
        simu.evolve(
            preOps = MapSelector(loci=0,
                    fitness={(0,0):1, (0,1):0.8, (1,1):1}),
                matingScheme = RandomMating(),
            postOps = [
                Stat( alleleFreq=[0], genoFreq=[0]),
                # PyEval('alleleFreq[0][0]'),
                TerminateIf('alleleFreq[0][0] > 0.4 and    alleleFreq[0][0]    < 0.6',
                    begin=50),
            ],
            initOps=[ InitSex(), InitGenotype(freq=[.5,.5])],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)

    def testMaSelectorHeteroDisadv(self):
        'Testing heterozygous advantage using a multi-allele selector'
        s1 = .1
        s2 = .2
        p = .2/ (.1+.2)
        # specify relative fitness: w11, w12/w21, w22
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        # 2. heterozygote inferiority
        #     w11 > w12, w12 < w22
        #    p unstable, become fix
        simu.evolve(
            preOps = MaSelector(loci=0, wildtype=0,
                    fitness=[1, 0.7, 1]),
                matingScheme = RandomMating(),
            postOps = [
                Stat( alleleFreq=[0], genoFreq=[0]),
                #PyEval('alleleFreq[0][0]'),
                TerminateIf('alleleFreq[0][0] > 0.3 and    alleleFreq[0][0]    < 0.7',
                    begin=50)
            ],
            initOps=[ InitSex(), InitGenotype(freq=[.5,.5])],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)

    def testMultiLocusMaSelector(self):
        'Testing the multi-locus version of the MaSelector'
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[3,6],
            infoFields=['fitness', 'spare']))
        simu.evolve(
            preOps = MaSelector(loci=[3,6], wildtype=0,
                    fitness=[1, 0.7, 1, 0.99, 0.98, 0.97, 1, 1, 0.5]),
                matingScheme = RandomMating(),
            initOps=[InitSex(), InitGenotype(freq=[.5,.5])],
            gen=100
        )

    def testMultiLocusMapSelector(self):
        'Testing multiple loci map selector.'
        pop = Population(10, loci=[2],
            infoFields=['fitness'])
        initGenotype(pop, genotype=[[0,0],[1,1]], proportions=[0.5,0.5])
        #MapSelector(pop, loci=[0,1],
        #    fitness={(0,0,0,0):0, (1,1,1,1):0.25, (0,1,0,1):0.5, (1,0,1,0):0.75})
        # there is only one field, so fitness is continuous
##         ft = pop.arrIndInfo()
##         for ind in range(pop.popSize()):
##             gt = pop.individual(ind).genotype()
##             if gt == [0,0,1,1]:
##                 self.assertEqual( ft[ind], 0.5)
##             # note that 10 => 01
##             elif gt == [1,1,0,0]:
##                 self.assertEqual( ft[ind], 0.5)
##             elif gt == [0,0,0,0]:
##                 self.assertEqual( ft[ind], 0)
##             elif gt == [1,1,1,1]:
##                 self.assertEqual( ft[ind], 0.25)
        # selector on a population with selection on is not allowed
        # explicitly walk around this.
        # test phase
        #MapSelector(pop, loci=[0,1], phase=True,
        #    fitness={(0,0,0,0):0, (1,1,1,1):0.25, (0,1,0,1):0.5, (1,0,1,0):0.75})
##         ft = pop.arrIndInfo()
##         for ind in range(pop.popSize()):
##             gt = pop.individual(ind).genotype()
##             if gt == [0,0,1,1]:
##                 self.assertEqual( ft[ind], 0.5)
##             # note that 10 != 01
##             elif gt == [1,1,0,0]:
##                 self.assertEqual( ft[ind], 0.75)
##             elif gt == [0,0,0,0]:
##                 self.assertEqual( ft[ind], 0)
##             elif gt == [1,1,1,1]:
##                 self.assertEqual( ft[ind], 0.25)

    def testPySelector(self):
        'Testing heterozygous advantage  using PySelector'
        s1 = .1
        s2 = .2
        p = .2/ (.1+.2)
        # gen may not be used.
        def sel(geno):
            if geno == (0, 0):
                return 1 - s1
            elif geno == (0, 1):
                return 1
            elif geno == (1, 0):
                return 1
            else:
                return 1 - s2
        #
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        # 2. heterozygote superiority
        #     w11 < w12, w12 > w22
        #    stable.
        # let
        #        s1 = w12 -    w11
        #        s2 = w12 - w22
        #    p_ = s2/ (s1+s2)
        simu.evolve(
            initOps = [
                InitSex(),
                InitGenotype(freq=[.5,.5])],
                matingScheme = RandomMating(),
            preOps = PySelector(loci=0, func=sel),
            postOps = [
                Stat( alleleFreq=[0], genoFreq=[0]),
                TerminateIf('alleleFreq[0][0] < 0.5', begin=50),
                TerminateIf('alleleFreq[0][0] > 0.9', begin=50)
            ],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)

    def testPySelectorWithGen(self):
        'Testing varying selection pressure using PySelector'
        s1 = .1
        s2 = .2
        p = .2/ (.1+.2)
        # gen may not be used.
        def sel(geno, gen):
            if gen > 50:
                if geno == (0, 0):
                    return 1 - s1
                elif geno == (0, 1):
                    return 1
                elif geno == (1, 0):
                    return 1
                else:
                    return 1 - s2
            else:
                if geno == (0, 0):
                    return 1 - s1/2.
                elif geno == (0, 1):
                    return 1
                elif geno == (1, 0):
                    return 1
                else:
                    return 1 - s2/2.
        #
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[1],
            infoFields=['fitness', 'spare']))
        # 2. heterozygote superiority
        #     w11 < w12, w12 > w22
        #    stable.
        # let
        #        s1 = w12 -    w11
        #        s2 = w12 - w22
        #    p_ = s2/ (s1+s2)
        simu.evolve(
            initOps = [
                InitSex(),
                InitGenotype(freq=[.5,.5])],
            preOps = PySelector(loci=[0], func=sel),
                matingScheme = RandomMating(),
            postOps = [
                Stat(alleleFreq=[0]),
                TerminateIf('alleleFreq[0][0] < 0.5', begin=50),
                TerminateIf('alleleFreq[0][0] > 0.9', begin=50)
            ],
            gen=100
        )
        # simulation did not terminate unexpectedly
        self.assertEqual(simu.dvars(0).gen, 100)

    def testMlSelector(self):
        'Testing multi-locus selector'
        simu = Simulator(
            Population(size=1000, ploidy=2, loci=[2],
            infoFields=['fitness', 'spare']))
        sel = MlSelector(
            [
                    MapSelector(loci=0, fitness={(0,0):1,(0,1):1,(1,1):.8}),
                    MapSelector(loci=1, fitness={(0,0):1,(0,1):1,(1,1):.8}),
                ], mode=ADDITIVE),
        simu.evolve(
            preOps = sel,
                matingScheme = RandomMating(),
            initOps = [ InitSex(), InitGenotype(freq=[.2,.8])],
            gen=100
        )
        #
        simu.evolve(preOps = [
            MlSelector([
                    MapSelector(loci=0, fitness={(0,0):1,(0,1):1,(1,1):.8}),
                    MaSelector(loci = 1, wildtype=[1], fitness=[1,1,.8])
                ], mode=MULTIPLICATIVE),
            ],
                matingScheme = RandomMating(),
            initOps=[ InitSex(), InitGenotype(freq=[.2,.8])],
            gen=100
        )

    def testSubPops(self):
        'Testing the subPops parameter of selector'
        simu = Simulator(
            Population(size=[20, 30, 40], loci=[2],
                infoFields='fitness'))
        def testFitness(pop, param):
            for sp in param[0]:
                for ind in pop.individuals(sp):
                    self.assertEqual(ind.info('fitness'), 0.)
        simu.evolve(
            initOps = [InitSex(), InitGenotype(freq=[.4, .6])],
            preOps = [
                MapSelector(loci = 1, fitness = {(0,0):1.,(0,1):1.,(1,1):.8}, subPops=[1]),
                PyOperator(func=testFitness, param=([0, 2],)),
                ],
                matingScheme = RandomMating(),
            gen = 5
        )
        # subPop is also allowed
        simu = Simulator(
            Population(size=[20, 30, 40], loci=[2],
                infoFields=['fitness']))
        simu.evolve(
            initOps = [InitSex(), InitGenotype(freq=[.4, .6])],
            preOps = [
                MaSelector(loci=0, wildtype=[0], fitness = [0.5, 0.4, 0.6],
                    subPops=1),
                MaSelector(loci=0, wildtype=[0], fitness = [0.6, 0.4, 0.6],
                    subPops=2),
                PyOperator(func=testFitness, param=([0],)),
                ],
                matingScheme = RandomMating(),
            gen = 5
        )



if __name__ == '__main__':
    unittest.main()
