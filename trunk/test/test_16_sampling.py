#!/usr/bin/env python
#
# Purpose:
#
#     Unittest testcases for ascertainment operators
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
from simuPOP.utils import *
from simuPOP.sampling import *
import unittest, os, sys, exceptions

class TestSampling(unittest.TestCase):

    def setUp(self):
        simu = simulator(
            population(size=[1000,2000], ploidy=2, loci=[5,10],
                ancGen=2,
                infoFields=['fitness', 'father_idx', 'mother_idx', 'migrate_to', 'oldindex']),
            randomMating(numOffspring=2))
        simu.evolve(
            preOps = migrator(rate=[[0.1,0.1], [0.1,0.1]]),
            duringOps = parentsTagger(),
            postOps = [
                stat( alleleFreq=[0,1], genoFreq=[0,1]),
                mapPenetrance(loci=0,
                    penetrance={(0,0):0,(0,1):.7,(1,1):1}),
            ],
            initOps = [
                 initSex(),
                 initByFreq(alleleFreq=[.2, .8], loci=[0]),
                 initByFreq(alleleFreq=[.2]*5, loci=range(1, simu.population(0).totNumLoci()))
            ],
            gen = 4
        )
        self.pop = simu.extract(0)
        for gen in range(self.pop.ancestralGens(), -1, -1):
            self.pop.useAncestralGen(gen)
            self.pop.setIndInfo(range(self.pop.popSize()), 'oldindex')
        # more complicated one
        simu1 = simulator(
            population(size=[5000,20000], ploidy=2, loci=[5,10],
                ancGen=2, infoFields=['fitness', 'father_idx', 'mother_idx', 'migrate_to', 'oldindex']),
            randomMating(numOffspring=(UniformDistribution, 2, 5)))
        simu1.evolve(
            preOps = migrator(rate=[[0.1,0.1],[0.1,0.1]]),
            duringOps = parentsTagger(),
            postOps = [
                stat( alleleFreq=[0,1], genoFreq=[0,1]),
                mapPenetrance(loci=0,
                    penetrance={(0,0):0,(0,1):.7,(1,1):1}),
            ],
            initOps=[
                 initSex(),
                 initByFreq(alleleFreq=[.2, .8], loci=[0]),
                 initByFreq(alleleFreq=[.2]*5, loci=range(1, simu1.population(0).totNumLoci()))
            ],
            gen = 10
        )
        self.largepop = simu1.extract(0)
        for gen in range(self.largepop.ancestralGens(), -1, -1):
            self.largepop.useAncestralGen(gen)
            self.largepop.setIndInfo(range(self.largepop.popSize()), 'oldindex')

    def testRandomSample(self):
        'Testing random sampling (imcomplete)'
        #
        s = DrawRandomSample(self.pop, 10)
        self.assertEqual(s.popSize(), 10)
        for ind in s.individuals():
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        # 
        s = DrawRandomSample(self.pop, [2, 8])
        self.assertEqual(s.subPopSize(0), 2)
        self.assertEqual(s.subPopSize(1), 8)
        #
        for ind in s.individuals():
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        #
        self.pop.setVirtualSplitter(sexSplitter())
        s = DrawRandomSample(self.pop, 10, subPops=[(0,0), (1,0)])
        # all samples should be Male
        self.assertEqual(s.popSize(), 10)
        for ind in s.individuals():
            self.assertEqual(ind.sex(), Male)
        #
        s = DrawRandomSample(self.pop, [2, 8], subPops=[(0,0), (1,0)])
        # all samples should be Male
        self.assertEqual(s.subPopSizes(), (2, 8))
        for ind in s.individuals():
            self.assertEqual(ind.sex(), Male)
        #
        samples = DrawRandomSamples(self.pop, [2, 8], subPops=[(0,0), (1,0)], times=10)
        self.assertEqual(len(samples), 10)
        for s in samples:
            # all samples should be Male
            self.assertEqual(s.subPopSizes(), (2, 8))
            for ind in s.individuals():
                self.assertEqual(ind.sex(), Male)

    def testCaseControlSample(self):
        'Testing case control sampling (imcomplete)'
        # case control sampling.
        (s,) = CaseControlSample(self.pop, 10, 10)
        self.assertEqual(s.subPopSize(0), 10)
        self.assertEqual(s.subPopSize(1), 10)
        #
        (s,) = CaseControlSample(self.pop, cases=[1,2], controls=[5,4])
        self.assertEqual(s.subPopSize(0), 3)
        self.assertEqual(s.subPopSize(1), 9)
        # # old index
        self.assertEqual(s.hasInfoField('oldindex'), True)
        #
        for ind in s.individuals(0):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        for ind in s.individuals(1):
            self.assertEqual(ind.affected(), False)
            #old index?
            inpop = self.pop.individual(int(ind.oldindex))


    def testAffectedSibpairSample(self):
        'Testing affected sibpair sampling (imcomplete)'
		# FIXME: testing sharing of parents
		# (father_idx and mother_idx of original and sample population,
		# and if the parents are the same.)
		#
        # find sibpairs
        (s,) = AffectedSibpairSample(self.pop, [2, 3])
        assert s.subPopSize(0) <= 4
        assert s.subPopSize(1) <= 6
        for ind in s.individuals(0):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        for ind in s.individuals(1):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        #
        (s,) = AffectedSibpairSample(self.pop, 2)
        s.useAncestralGen(1)
        self.assertEqual(s.subPopSizes(), (2, 2))
        s.useAncestralGen(0)
        self.assertEqual(s.subPopSizes(), (2, 2))
        for ind in s.individuals():
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)


    def TestLargePedigreeSample(self):
        'Testing large pedigree sampling (FIXME)'
        (s,) = LargePedigreeSample(self.largepop, minTotalSize=20, maxOffspring=5,
            minPedSize=3, minAffected=0)
        assert s.ancGen() == 2
        for ind in s.individuals():
            #old index?
            inpop = self.largepop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        (s,) = LargePedigreeSample(self.largepop, 50, minTotalSize=5, maxOffspring=5,
            minPedSize=3, minAffected=0)
        for ind in s.individuals():
            #old index?
            inpop = self.largepop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        #


    def TestNuclearFamilySample(self):
        'Testing nuclear family sampling (imcomplete)'
        (s,) = NuclearFamilySample(self.largepop, 50, minTotalSize=50, maxOffspring=5,
            minPedSize=5, minAffected=0)
        #print s.subPopSizes()
        assert s.subPopSize(0) <= 5
        assert s.subPopSize(1) <= 5
        for ind in s.individuals():
            #old index?
            inpop = self.largepop.individual(ind.oldindex)
            self.assertEqual(ind, inpop)
        #
        (s,) = NuclearFamilySample(self.largepop, size=[2, 3], maxOffspring=5)
        assert s.subPopSize(0) <= 5
        for ind in s.individuals():
            #old index?
            inpop = self.largepop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)



if __name__ == '__main__':
    unittest.main()
