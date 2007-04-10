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
import unittest, os, sys, exceptions

class TestAscertainment(unittest.TestCase):

    def setUp(self):
        if alleleType() == 'binary':
            simu = simulator(
                population(subPop=[1000,2000], ploidy=2, loci=[5,10],
                    ancestralDepth=2, maxAllele=1, 
                    infoFields=['fitness', 'father_idx', 'mother_idx']),
                randomMating(numOffspring=2))
        else:
            simu = simulator(
                population(subPop=[1000,2000], ploidy=2, loci=[5,10],
                    ancestralDepth=2, maxAllele=9, 
                    infoFields=['fitness', 'father_idx', 'mother_idx']),
                randomMating(numOffspring=2))
        simu.evolve(
            [
                stat( alleleFreq=[0,1], genoFreq=[0,1]),
                migrator(rate=[[0.1,0.1],[0.1,0.1]]),
                mapPenetrance(locus=0,
                    penetrance={'0-0':0,'0-1':.7,'1-1':1}),
                parentsTagger(),
            ],
            preOps=[    
                 initByFreq(alleleFreq=[.2, .8], atLoci=[0]),
                 initByFreq(alleleFreq=[.2]*5, atLoci=range(1, simu.totNumLoci()))     
            ],
            end=4
        )
        self.pop = simu.getPopulation(0)
    
    def testRandomSample(self):
        'Testing random sampling (imcomplete)'
        (s,) = RandomSample(self.pop, 10)
        self.assertEqual(s.popSize(), 10)
        #
        (s,) = RandomSample(self.pop,[2,8])
        self.assertEqual(s.subPopSize(0), 2)
        self.assertEqual(s.subPopSize(1), 8)
        for ind in s.individuals():
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        
    def testCaseControlSample(self):
        'Testing case control sampling (imcomplete)'
        # case control sampling.
        (s,) = CaseControlSample(self.pop, 10, 10)
        self.assertEqual(s.subPopSize(0), 10)
        self.assertEqual(s.subPopSize(1), 10)
        #
        (s,) = CaseControlSample(self.pop, [1,2],[5,4])
        self.assertEqual(s.subPopSize(0), 3)
        self.assertEqual(s.subPopSize(1), 9)
        # old index
        self.assertEqual(s.hasInfoField('oldindex'), True)
        #
        for ind in s.individuals(0):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        for ind in s.individuals(1):
            self.assertEqual(ind.affected(), False)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))

    def testAffectedSibpairSample(self):
        'Testing affected sibpair sampling (imcomplete)'
		# FIXME: testing sharing of parents 
		# (father_idx and mother_idx of original and sample population,
		# and if the parents are the same.)
		#
        # find sibpairs
        (s,) = AffectedSibpairSample(self.pop, [2,3])
        assert s.subPopSize(0) <= 4
        assert s.subPopSize(1) <= 6
        for ind in s.individuals(0):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        for ind in s.individuals(1):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        #
        (s,) = AffectedSibpairSample(self.pop, 2)
        assert s.subPopSize(0) <= 4
        for ind in s.individuals():
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)


    def testLargePedigreeSample(self):
        'Testing large pedigree sampling (FIXME)'
        (s,) = LargePedigreeSample(self.pop, 5, minTotalSize=10, maxOffspring=5, 
            minPedSize=3, minAffected=0)
        print s.subPopSizes()
        return
        assert s.subPopSize(0) <= 4
        assert s.subPopSize(1) <= 6
        for ind in s.individuals(0):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        for ind in s.individuals(1):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        #
        (s,) = AffectedSibpairSample(self.pop, 2)
        assert s.subPopSize(0) <= 4
        for ind in s.individuals():
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        pass


    def testNuclearFamilySample(self):
        'Testing nuclear family sampling'
        (s,) = NuclearFamilySample(self.pop, 5, )
        return
        assert s.subPopSize(0) <= 4
        assert s.subPopSize(1) <= 6
        for ind in s.individuals(0):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        for ind in s.individuals(1):
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)
        #
        (s,) = AffectedSibpairSample(self.pop, 2)
        assert s.subPopSize(0) <= 4
        for ind in s.individuals():
            self.assertEqual(ind.affected(), True)
            #old index?
            inpop = self.pop.individual(int(ind.info('oldindex')))
            self.assertEqual(ind, inpop)

        

if __name__ == '__main__':
    unittest.main()
