#!/usr/bin/env python
#
#  This is a unittest file for operator recombinator
#
#  Bo Peng (bpeng@rice.edu)
#
# Last Modified: Sep, 2005
#

from simuPOP import *
from simuUtil import *

import unittest

class TestRecombinator(unittest.TestCase):
  def testNoNewAllele(self):
    ' cross same haplotype should not result in new alleles at loci'
    geno = [1,2,3,4,5,6,7]
    pop = population(1000, loci=[2,3,2])
    InitByValue(pop, value=geno)
    simu = simulator(pop, randomMating())
    simu.evolve( [ recombinator(rate = 0.4) ], end=100)
    Stat(simu.population(0), alleleFreq=range(0,7))
    for i in range(7):
      self.assertEqual(simu.dvars(0).alleleFreq[i][i+1], 1.)
    
  def testCrossBetweenChrom(self):
    ' see if chromsomes are crossed by default'
    h1 = [1,1,1,1,1,1,1]
    h2 = [2,2,2,2,2,2,2]
    pop = population(10000, loci=[2,3,2])
    InitByValue(pop, value=[h1,h2], proportions=[0.5,.5])
    simu = simulator(pop, randomMating())
    # no recombination within chromosome
    TurnOnDebug(DBG_RECOMBINATOR)
    simu.step( [ 
      stat( haploFreq = [[0,1], [2,3], [2,4], [5,6], [0,2], [0,6], [3,6]],
        stage=PrePostMating), 
      pyEval('haploFreq', stage=PrePostMating),
      recombinator(rate = 0) ] )
#    Stat(simu.population(0), 
      # first four: within chrom. Only original haplotypes
      # last three: new haplotypes
#      haploFreq = [[0,1], [2,3], [2,4], [5,6], [0,2], [0,6], [3,6]])
    print simu.dvars(0).haploFreq
    self.assertEqual(simu.dvars(0).haploFreq['0-1'].setdefault('1-2', 0), 0.)
    self.assertEqual(simu.dvars(0).haploFreq['2-3'].setdefault('1-2', 0), 0.)
    self.assertEqual(simu.dvars(0).haploFreq['2-4'].setdefault('1-2', 0), 0.)
    self.assertEqual(simu.dvars(0).haploFreq['5-6'].setdefault('1-2', 0), 0.)
    assert abs(simu.dvars(0).haploFreq['0-2'].setdefault('1-2', 0) - 0.25) < 0.1
    assert abs(simu.dvars(0).haploFreq['0-6'].setdefault('1-2', 0) - 0.25) < 0.1
    assert abs(simu.dvars(0).haploFreq['3-6'].setdefault('1-2', 0) - 0.25) < 0.1
   
    
  def testRecombine(self):
    ' verify table 4 of H&C 3nd edition P49 '
    N = 10000
    r = 0.1
    
    genoDad = [[1,1],[1,1],[1,1],[1,1],[1,2],[1,2],[1,2],[2,1],[2,1],[2,2]]
    genoMom = [[1,1],[1,2],[2,1],[2,2],[1,2],[2,1],[2,2],[2,1],[2,2],[2,2]]
    prop = [ [1, 0, 0, 0], [.5, .5, 0, 0], [.5, .5, 0, 0],
      [0.5-r/2, r/2, r/2, 0.5-r/2], [0, 1, 0, 0], [r/2, .5-r/2, .5-r/2, r/2], 
      [0, .5, 0, .5], [0, 1, 0, 0], [0, 0, .5, .5], [0, 0, 0, 1] ]
    
    for i in range(0, len(genoDad)):
      pop = population(size=N, loci=[2])
      InitByValue(pop, value=[genoDad[i], genoMom[i]], \
        proportions=[0.5,0.5])
      simu = simulator(pop, randomMating())
      simu.step(
          [ recombinator(rate=r),
            stat(haploFreq=[0,1]) 
          ])
      hf = simu.dvars(0).haploFreq['0-1']
      assert not ((hf.setdefault('1-1',0) - prop[i][0]) > 0.1 or \
        (hf.setdefault('1-2',0) - prop[i][1]) > 0.1 or \
        (hf.setdefault('2-1',0) - prop[i][2]) > 0.1 or \
        (hf.setdefault('2-2',0) - prop[i][3]) > 0.1), \
        "Recombination results in potentially wrong proportions."
        
      
    # table is correct.
    # 1/2*(1-r)=0.5*0.9=0.45
    
    
#  def testLDDecay(self):
#    ' verify Dn=(1-r)^n D0 '
#    i=1
#    pop = population(size=N, loci=[5,7,8], sexChrom=True)
#    InitByValue(pop, value=[genoDad[i], genoMom[i]], proportions=[0.5,0.5])
#    simu = simulator(pop, randomMating())
#    simu.evolve(
#        [ recombinator(rate=r, maleRate=r/2., maleAfterLoci=[2,6,7]),
#          stat(haploFreq=[[0,1]]) ],
#        end=100)


if __name__ == '__main__':
  unittest.main()
   


