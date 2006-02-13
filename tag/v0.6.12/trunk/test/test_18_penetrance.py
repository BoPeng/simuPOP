#!/usr/bin/env python
#
# Purpose:
#   Testing penetrance.
# 
# Author:
#   Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
# 

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

class TestPenetrance(unittest.TestCase):
  
  def setUp(self):
    self.pop = population(subPop=[500,100,1000], 
      ploidy=2, loci = [1])
    if alleleType() == 'binary':
      InitByValue(self.pop, 
        value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
        indRange = [[0,125], [125,375],[375,500],[500,550],
          [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
      self.a1, self.a2 = 0, 1
      self.key11, self.key12, self.key22 = '0-0', '0-1', '1-1'
      self.wildtype = 0
    else:
      InitByValue(self.pop, 
        value = [[1,1],[1,2],[2,2],[1,1],[1,2],[2,2],[1,2],[1,2],[2,2]],
        indRange = [[0,125], [125,375],[375,500],[500,550],
          [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
      self.a1, self.a2 = 1, 2
      self.key11, self.key12, self.key22 = '1-1', '1-2', '2-2'
      self.wildtype = 1
 
  def testMapPenetrance(self):
    'Testing map penetrance'
    MapPenetrance(self.pop, locus=0, 
      penetrance={self.key11:0, self.key12:1, self.key22:1})
    Stat(self.pop, numOfAffected=1)
    self.assertEqual(self.pop.dvars().numOfAffected, 1425)
    self.assertEqual(self.pop.dvars(0).numOfAffected, 375)
    self.assertEqual(self.pop.dvars(1).numOfAffected, 50)
    self.assertEqual(self.pop.dvars(2).numOfAffected, 1000)
    #
    # imcomlete penetrance
    MapPenetrance(self.pop, locus=0, 
      penetrance={self.key11:0, self.key12:.3, self.key22:.5})
    Stat(self.pop, numOfAffected=1)
    assert abs(self.pop.dvars().numOfAffected -  880*0.3 - 545*0.5) < 100
    assert abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5) < 30
    assert abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5) < 15
    assert abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5) < 50
    
  def testMaPenetrance(self):
    'Testing multi-allele penetrance'
    MaPenetrance(self.pop, locus=0, wildtype=self.wildtype,
      penetrance=[0, 1, 1])
    Stat(self.pop, numOfAffected=1)
    self.assertEqual(self.pop.dvars().numOfAffected, 1425)
    self.assertEqual(self.pop.dvars(0).numOfAffected, 375)
    self.assertEqual(self.pop.dvars(1).numOfAffected, 50)
    self.assertEqual(self.pop.dvars(2).numOfAffected, 1000)
    #
    # imcomlete penetrance
    self.pop.dvars().clear()
    MaPenetrance(self.pop, locus=0,  wildtype=self.wildtype,
      penetrance=[0, .3, .5])
    Stat(self.pop, numOfAffected=1)
    assert abs(self.pop.dvars().numOfAffected -  880*0.3 - 545*0.5) < 100
    assert abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5) < 30
    assert abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5) < 15
    assert abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5) < 50
    
  def testMultiLocusMaPenetrance(self):
    'Test the multi-locus version of maPenetrance'
    pop = population(1000, loci=[3,5])
    InitByFreq(pop, [.3, .7])
    #
    MaPenetrance(pop, loci=[3,5], wildtype=self.wildtype,
      penetrance=[0, .3, .5, 0.3, 0.6, 0.8, 0.1, 1, 0.8])

    
  def testMlPenetrance(self):
    'Testing multi-locus penetrance'
    pop = population(1000, loci=[3,5])
    InitByFreq(pop, [.3, .7])
    #
    MlPenetrance(pop, [
      maPenetrance(locus=0,  wildtype=self.wildtype,
        penetrance=[0, .3, .5]),
      mapPenetrance(locus=1, 
        penetrance={self.key11:0, self.key12:1, self.key22:1})
      ],
      mode=PEN_Additive
    )
    #
    MlPenetrance(pop, [
      maPenetrance(locus=2,  wildtype=self.wildtype,
        penetrance=[0, .3, .5]),
      mapPenetrance(locus=4, 
        penetrance={self.key11:0, self.key12:1, self.key22:1})
      ],
      mode=PEN_Multiplicative
    )

  def testPyPenetrance(self):
    def pen(geno):
      if geno == [self.a1, self.a1]:
        return 0
      elif geno == [self.a1, self.a2]:
        return 0.5
      elif geno == [self.a2, self.a1]:
        return 0.5
      else:
        return 1
    PyPenetrance(self.pop, locus=0,  func=pen)
    Stat(self.pop, numOfAffected=1)
    assert abs(self.pop.dvars().numOfAffected -  880*0.5 - 545) < 100
    assert abs(self.pop.dvars(0).numOfAffected - 250*0.5 - 125) < 30
    assert abs(self.pop.dvars(1).numOfAffected - 30*0.5 - 20) < 15
    assert abs(self.pop.dvars(2).numOfAffected - 600*0.5 - 400) < 50
    
if __name__ == '__main__':
  unittest.main()    
