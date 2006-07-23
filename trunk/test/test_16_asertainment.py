#!/usr/bin/env python
#
# Purpose:
#
#   Unittest testcases for ascertainment operators
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
    simu = simulator(
      population(subPop=[100,200], ploidy=2, loci=[5,10],
        ancestralDepth=1, maxAllele=MaxAllele),
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
         initByFreq(alleleFreq=[.2,.8], atLoci=[0]),
         initByFreq(alleleFreq=[.2]*5, atLoci=range(1, simu.totNumLoci()))   ],
       end=4
    )
    self.pop = simu.getPopulation(0)
  
  def testRandomSamle(self):
    'Testing random sampling (imcomplete)'
    (s,) = RandomSample(self.pop, 10)
    self.assertEqual(s.popSize(), 10)
    #
    (s,) = RandomSample(self.pop,[2,8])
    self.assertEqual(s.subPopSize(0), 2)
    self.assertEqual(s.subPopSize(1), 8)
    
  def testCaseControlSample(self):
    'Testing case control sampling (imcomplete)'
    # case control sampling.
    (s,) = CaseControlSample(self.pop, 10, 10)
    self.assertEqual(s.subPopSize(0), 10)
    self.assertEqual(s.subPopSize(1), 10)
    #
    (s,) =  CaseControlSample(self.pop, [1,2],[5,4])
    self.assertEqual(s.subPopSize(0), 3)
    self.assertEqual(s.subPopSize(1), 9)

  def testAffectedSibpairSample(self):
    'Testing affected sibpair sampling (imcomplete)'
    # find sibpairs
    (s,) = AffectedSibpairSample(self.pop, [2,3])
    assert s.subPopSize(0) <= 4
    assert s.subPopSize(1) <= 6
    # 
    (s,) = AffectedSibpairSample(self.pop, 2)
    assert s.subPopSize(0) <= 4
    

if __name__ == '__main__':
  unittest.main()
