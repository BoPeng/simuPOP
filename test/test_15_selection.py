#!/usr/bin/env python
#
# Purpose:
#   Testing selection.
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

class TestSelector(unittest.TestCase):
  def setUp(self):
    if alleleType() == 'binary':
      self.a1, self.a2 = 0, 1
      self.key11, self.key12, self.key22 = '0-0', '0-1', '1-1'
      self.freqStr = 'alleleFreq[0][0]'
      self.wildtype = 0
    else:
      self.a1, self.a2 = 1, 2
      self.key11, self.key12, self.key22 = '1-1', '1-2', '2-2'
      self.freqStr = 'alleleFreq[0][1]'
      self.wildtype = 1

  def testMapSelectorDirSelection(self):
    'Testing directional selection using a map selector'
    # specify relative fitness: w11, w12/w21, w22
    simu = simulator(
      population(size=1000, ploidy=2, loci=[1]),
      randomMating() )
    # 1. directional selection
    #   w11 > w12 > w22
    #  p -> 1
    simu.evolve(
      [
        stat( alleleFreq=[0], genoFreq=[0]),
        mapSelector(locus=0, 
          fitness={self.key11:1, self.key12:0.9, self.key22:.8}),
        terminateIf(self.freqStr + '< 0.4'),
        terminateIf(self.freqStr + '< 0.8', begin=50)
      ],
      preOps=[ initByFreq(alleleFreq=[.5,.5])],
      end=100
    )
    # simulation did not terminate unexpectedly
    self.assertEqual(simu.gen(), 101)
    
  def testMaSelectorDirSelection(self):
    'Testing directional selection using a multi-allele selector'
    # specify relative fitness: w11, w12/w21, w22
    simu = simulator(
      population(size=1000, ploidy=2, loci=[1]),
      randomMating() )
    # 1. directional selection
    #   w11 > w12 > w22
    #  p -> 1
    simu.evolve(
      [
        stat( alleleFreq=[0], genoFreq=[0]),
        maSelector(locus=0, wildtype=[self.wildtype], 
          fitness=[1, 0.9, .8]),
        terminateIf(self.freqStr + '< 0.4'),
        terminateIf(self.freqStr + '< 0.8', begin=50)
      ],
      preOps=[ initByFreq(alleleFreq=[.5,.5])],
      end=100
    )
    # simulation did not terminate unexpectedly
    self.assertEqual(simu.gen(), 101)
    
  def testMapSelectorHeteroAdv(self):
    'Testing heterozygous advantage using map selector'
    # specify relative fitness: w11, w12/w21, w22
    simu = simulator(
      population(size=1000, ploidy=2, loci=[1]),
      randomMating() )
    s1 = .1
    s2 = .2
    p = .2/ (.1+.2)
    # 2. heterozygote superiority
    #   w11 < w12, w12 > w22
    #  stable. 
    # let
    #    s1 = w12 -  w11
    #    s2 = w12 - w22
    #  p_ = s2/ (s1+s2)
    simu.evolve(
      [
        stat( alleleFreq=[0], genoFreq=[0]),
        mapSelector(locus=0, 
          fitness={self.key11:1-s1, self.key12:1, self.key22:1-s2}),
        terminateIf(self.freqStr + '< 0.5', begin=50),
        terminateIf(self.freqStr + '> 0.9', begin=50)
      ],
      preOps=[ initByFreq(alleleFreq=[.5,.5])],
      end=100
    )
    # simulation did not terminate unexpectedly
    self.assertEqual(simu.gen(), 101)
    
  def testMaSelectorHeteroAdv(self):
    'Testing heterozygous advantage using a multi-allele selector'
    s1 = .1
    s2 = .2
    p = .2/ (.1+.2)
    # specify relative fitness: w11, w12/w21, w22
    simu = simulator(
      population(size=1000, ploidy=2, loci=[1]),
      randomMating() )
    # 2. heterozygote superiority
    #   w11 < w12, w12 > w22
    #  stable. 
    # let
    #    s1 = w12-  w11
    #    s2 = w12 - w22
    #  p_ = s2/ (s1+s2)
    #
    simu.evolve(
      [
        stat( alleleFreq=[0], genoFreq=[0]),
        maSelector(locus=0, wildtype=self.wildtype, 
          fitness=[1-s1, 1, 1-s2]),
        terminateIf(self.freqStr + '< 0.5', begin=50),
        terminateIf(self.freqStr + '> 0.9', begin=50)
      ],
      preOps=[ initByFreq(alleleFreq=[.5,.5])],
      end=100
    )
    # simulation did not terminate unexpectedly
    self.assertEqual(simu.gen(), 101)

  def testMapSelectorHeteroDisadv(self):
    'Testing heterozygous disadvantage using map selector'
    simu = simulator(
      population(size=1000, ploidy=2, loci=[1]),
      randomMating() )
    # 2. heterozygote inferiority
    #   w11 > w12, w12 < w22
    #  p unstable, become fix
    simu.evolve(
      [
        stat( alleleFreq=[0], genoFreq=[0]),
        mapSelector(locus=0, 
          fitness={self.key11:1, self.key12:0.8, self.key22:1}),
        # pyEval(self.freqStr),
        terminateIf(self.freqStr + '> 0.4 and ' + self.freqStr + ' < 0.6', 
          begin=50),
      ],
      preOps=[ initByFreq(alleleFreq=[.5,.5])],
      end=100
    )
    # simulation did not terminate unexpectedly
    self.assertEqual(simu.gen(), 101)
    
  def testMaSelectorHeteroDisadv(self):
    'Testing heterozygous advantage using a multi-allele selector'
    s1 = .1
    s2 = .2
    p = .2/ (.1+.2)
    # specify relative fitness: w11, w12/w21, w22
    simu = simulator(
      population(size=1000, ploidy=2, loci=[1]),
      randomMating() )
    # 2. heterozygote inferiority
    #   w11 > w12, w12 < w22
    #  p unstable, become fix
    simu.evolve(
      [
        stat( alleleFreq=[0], genoFreq=[0]),
        maSelector(locus=0, wildtype=self.wildtype, 
          fitness=[1, 0.7, 1]), 
        #pyEval(self.freqStr),
        terminateIf(self.freqStr + '> 0.3 and ' + self.freqStr + ' < 0.7', 
          begin=50)
      ],
      preOps=[ initByFreq(alleleFreq=[.5,.5])],
      end=100
    )
    # simulation did not terminate unexpectedly
    self.assertEqual(simu.gen(), 101)

  def testMultiLocusMaSelector(self):
    'Testing the multi-locus version of the maSelector'
    simu = simulator(
      population(size=1000, ploidy=2, loci=[3,6]),
      randomMating() )
    simu.evolve(
      [
        maSelector(loci=[3,6], wildtype=self.wildtype, 
          fitness=[1, 0.7, 1, 0.99, 0.98, 0.97, 1, 1, 0.5]), 
      ],
      preOps=[ initByFreq(alleleFreq=[.5,.5])],
      end=100
    )
  
  def testMultiLocusMapSelector(self):
    'Testing basic parameters of selector'
    pop = population(10, loci=[2])
    InitByValue(pop, value=[[0,0],[1,1]], proportions=[0.5,0.5])
    mapSelector(loci=[0,1], 
      fitness={'0-0|0-0':0, '1-1|1-1':0.25,
               '0-1|0-1':0.5, '1-0|1-0':0.75}).apply(pop)
    ft = pop.arrFitness()
    for ind in range(pop.popSize()):
      gt = pop.individual(ind).arrGenotype()
      if gt == [0,0,1,1]:
        self.assertEqual( ft[ind], 0.5)
      # note that 10 => 01
      elif gt == [1,1,0,0]:
        self.assertEqual( ft[ind], 0.5)
      elif gt == [0,0,0,0]:
        self.assertEqual( ft[ind], 0)
      elif gt == [1,1,1,1]:
        self.assertEqual( ft[ind], 0.25)
    # test phase
    mapSelector(loci=[0,1], phase=True,
      fitness={'0-0|0-0':0, '1-1|1-1':0.25,
               '0-1|0-1':0.5, '1-0|1-0':0.75}).apply(pop)
    ft = pop.arrFitness()
    for ind in range(pop.popSize()):
      gt = pop.individual(ind).arrGenotype()
      if gt == [0,0,1,1]:
        self.assertEqual( ft[ind], 0.5)
      # note that 10 != 01
      elif gt == [1,1,0,0]:
        self.assertEqual( ft[ind], 0.75)
      elif gt == [0,0,0,0]:
        self.assertEqual( ft[ind], 0)
      elif gt == [1,1,1,1]:
        self.assertEqual( ft[ind], 0.25)
    
  def testPySelector(self):
    'Testing heterozygous advantage  using pySelector'
    s1 = .1
    s2 = .2
    p = .2/ (.1+.2)
    def sel(arr):
      if arr == [self.a1, self.a1]:
        return 1 - s1
      elif arr == [self.a1, self.a2]:
        return 1
      elif arr == [self.a2, self.a1]:
        return 1
      else:
        return 1 - s2
    #
    simu = simulator(
      population(size=1000, ploidy=2, loci=[1]),
      randomMating() )
    # 2. heterozygote superiority
    #   w11 < w12, w12 > w22
    #  stable. 
    # let
    #    s1 = w12 -  w11
    #    s2 = w12 - w22
    #  p_ = s2/ (s1+s2)
    simu.evolve(
      [
        stat( alleleFreq=[0], genoFreq=[0]),
        pySelector(loci=[0], func=sel),
        terminateIf(self.freqStr + '< 0.5', begin=50),
        terminateIf(self.freqStr + '> 0.9', begin=50)
      ],
      preOps=[ initByFreq(alleleFreq=[.5,.5])],
      end=100
    )
    # simulation did not terminate unexpectedly
    self.assertEqual(simu.gen(), 101)
    

  def testMlSelector(self):
    'Testing multi-locus selector'
    simu = simulator(
      population(size=1000, ploidy=2, loci=[2]),
      randomMating())
    simu.evolve(
      [
        mlSelector([
          mapSelector(locus=0, fitness={self.key11:1,self.key12:1,self.key22:.8}),
          mapSelector(locus=1, fitness={self.key11:1,self.key12:1,self.key22:.8}),
        ], mode=SEL_Additive),
      ],
      preOps=[ initByFreq(alleleFreq=[.2,.8])],
      end=100
    )
    # 
    simu.setGen(0)
    simu.evolve([
      mlSelector(
        [
          mapSelector(locus=0, fitness={self.key11:1,self.key12:1,self.key22:.8}),
          maSelector(locus=1, wildtype=[1], fitness=[1,1,.8])
        ], mode=SEL_Multiplicative),
      ],
      preOps=[ initByFreq(alleleFreq=[.2,.8])],
      end=100
    )
    
if __name__ == '__main__':
  unittest.main()    
