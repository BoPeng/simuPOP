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
  def testRecRate(self):
    ' see if we actually recombine at this rate '
    pop = population(10000, loci=[2,3,2])
    InitByValue(pop, value=[1]*7+[2]*7)
    simu = simulator(pop, randomMating())
    simu.step( [ 
      stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]), 
      recombinator(rate = 0.1) ] )
    # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
    #print simu.dvars(0).haploFreq
    assert (simu.dvars(0).haploFreq['0-1']['1-2'] - 0.05) < 0.01
    assert (simu.dvars(0).haploFreq['2-3']['1-2'] - 0.05) < 0.01
    assert (simu.dvars(0).haploFreq['3-4']['1-2'] - 0.05) < 0.01
    assert (simu.dvars(0).haploFreq['4-5']['1-2'] - 0.25) < 0.01
    assert (simu.dvars(0).haploFreq['5-6']['1-2'] - 0.05) < 0.01
    
  def testRecRates(self):
    ' see if we actually recombine at this rate '
    pop = population(10000, loci=[2,3,2])
    InitByValue(pop, value=[1]*7+[2]*7)
    simu = simulator(pop, randomMating())
    simu.step( [ 
      stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]), 
      recombinator(rate = [0.1,0.15,0.2,0.3], afterLoci=[0,2,3,5] ) ] )
    # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
    #print simu.dvars(0).haploFreq
    assert (simu.dvars(0).haploFreq['0-1']['1-2'] - 0.05) < 0.01
    assert (simu.dvars(0).haploFreq['2-3']['1-2'] - 0.075)< 0.01
    assert (simu.dvars(0).haploFreq['3-4']['1-2'] - 0.1)  < 0.01
    assert (simu.dvars(0).haploFreq['4-5']['1-2'] - 0.25) < 0.01
    assert (simu.dvars(0).haploFreq['5-6']['1-2'] - 0.15) < 0.01
          
  def testRecIntensity(self):
    ' see if we actually recombine at this rate '
    pop = population(10000, loci=[2,3,2], lociDist=[0,1,0,2,4,0,4] )
    InitByValue(pop, value=[1]*7+[2]*7)
    simu = simulator(pop, randomMating())
    simu.step( [ 
      stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]), 
      recombinator(intensity = 0.1) ] )
    # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
    #print simu.dvars(0).haploFreq
    assert (simu.dvars(0).haploFreq['0-1']['1-2'] - 0.05) < 0.01
    assert (simu.dvars(0).haploFreq['2-3']['1-2'] - 0.1)  < 0.01
    assert (simu.dvars(0).haploFreq['3-4']['1-2'] - 0.1)  < 0.01
    assert (simu.dvars(0).haploFreq['4-5']['1-2'] - 0.25) < 0.01
    assert (simu.dvars(0).haploFreq['5-6']['1-2'] - 0.2) < 0.01
  
  def testRecIntensityAfterLoci(self):
    ' see if we actually recombine at this rate '
    pop = population(10000, loci=[2,3,2], lociDist=[0,1,0,2,4,0,4] )
    InitByValue(pop, value=[1]*7+[2]*7)
    simu = simulator(pop, randomMating())
    simu.step( [ 
      stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]), 
      recombinator(intensity = 0.1, afterLoci=[2,5]) ] )
    # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
    #print simu.dvars(0).haploFreq
    assert simu.dvars(0).haploFreq['0-1'].setdefault('1-2',0) == 0
    assert (simu.dvars(0).haploFreq['2-3']['1-2'] - 0.1)  < 0.01
    assert simu.dvars(0).haploFreq['3-4'].setdefault('1-2',0) == 0
    assert (simu.dvars(0).haploFreq['4-5']['1-2'] - 0.25) < 0.01
    assert (simu.dvars(0).haploFreq['5-6']['1-2'] - 0.2) < 0.01    
    
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
    pop = population(10000, loci=[2,3,2])
    InitByValue(pop, value=[1]*7+[2]*7)
    simu = simulator(pop, randomMating())
    #TurnOnDebug(DBG_RECOMBINATOR)
    simu.step( [ 
      stat( haploFreq = [[0,1], [2,3], [2,4], [5,6], [0,2], [0,6], [3,6]],
        stage=PrePostMating), 
      ## for debug purpose, output haploFreq
      #pyEval('haploFreq', stage=PrePostMating),
      recombinator(rate = 0) ] )
    self.assertEqual(simu.dvars(0).haploFreq['0-1'].setdefault('1-2', 0), 0.)
    self.assertEqual(simu.dvars(0).haploFreq['2-3'].setdefault('1-2', 0), 0.)
    self.assertEqual(simu.dvars(0).haploFreq['2-4'].setdefault('1-2', 0), 0.)
    self.assertEqual(simu.dvars(0).haploFreq['5-6'].setdefault('1-2', 0), 0.)
    assert abs(simu.dvars(0).haploFreq['0-2'].setdefault('1-2', 0) - 0.25) < 0.01
    assert abs(simu.dvars(0).haploFreq['0-6'].setdefault('1-2', 0) - 0.25) < 0.01
    assert abs(simu.dvars(0).haploFreq['3-6'].setdefault('1-2', 0) - 0.25) < 0.01
    
  def testRecProportion(self):
    ' verify table 4 of H&C 3nd edition P49 '
    N = 10000
    r = 0.1
    
    genoDad = [[1,1],[1,1],[1,1],[1,1],[1,2],[1,2],[1,2],[2,1],[2,1],[2,2]]
    genoMom = [[1,1],[1,2],[2,1],[2,2],[1,2],[2,1],[2,2],[2,1],[2,2],[2,2]]
    prop = [ [1, 0, 0, 0], [.5, .5, 0, 0], [.5, 0, 0.5, 0],
      [0.5-r/2, r/2, r/2, 0.5-r/2], [0, 1, 0, 0], [r/2, .5-r/2, .5-r/2, r/2], 
      [0, .5, 0, .5], [0, 0, 1, 0], [0, 0, .5, .5], [0, 0, 0, 1] ]
    
    for i in range(0, len(genoDad)):
      pop = population(size=N, loci=[2])
      InitByValue(pop, value=genoDad[i]+genoMom[i])
      simu = simulator(pop, randomMating())
      simu.step(
          [ recombinator(rate=r),
            stat(haploFreq=[0,1]) 
          ])
      hf = simu.dvars(0).haploFreq['0-1']
      assert not ((hf.setdefault('1-1',0) - prop[i][0]) > 0.01 or \
        (hf.setdefault('1-2',0) - prop[i][1]) > 0.01 or \
        (hf.setdefault('2-1',0) - prop[i][2]) > 0.01 or \
        (hf.setdefault('2-2',0) - prop[i][3]) > 0.01), \
        "Recombination results in potentially wrong proportions." +  \
        str( genoDad[i]) + ' crossing ' + str(genoMom[i])
      
  def testLDDecay(self):
    ' verify Dn=(1-r)^n D0 '
    r = 0.1
    N = 10000
    pop = population(size=N, loci=[2])
    # genotype 11/22, with LD=1
    InitByValue(pop, value=[1,1,2,2])
    simu = simulator(pop, randomMating())
    simu.evolve( 
      ops  = [ recombinator(rate=r) ],
      postOps = [  stat(LD=[0,1]) ], 
      end=9)
    # check the change of LD, hopefully, the variation is not too high.
    assert abs(simu.dvars(0).LD[0][1] - 0.25*(1-r)**10) < 0.02, \
      "Decay of LD is not as expected: " + str(simu.dvars(0).LD[0][1]) + " vs expected " \
      + str( 0.25*(1-r)**10 )

  def testNoMaleRec(self):
    ' male chromosome is not supposed to recombine '
    # create such an situation where female has 11111/22222, male has 1111/33333
    # we will test if 3333 is untouched under recombination
    r = 0.1
    N = 100
    pop = population(size=N, loci=[2,5], sexChrom=True)
    InitByValue(pop, indRange=[0,N/2-1], sex=[Male]*(N/2), atPloidy=0, value=[1]*7)
    InitByValue(pop, indRange=[0,N/2-1], sex=[Male]*(N/2), atPloidy=1, value=[3]*7)
    InitByValue(pop, indRange=[N/2,N-1], sex=[Female]*(N/2), value=[1]*7+[2]*7)
    # now let us recombine
    simu = simulator(pop, randomMating())
    simu.evolve( [ recombinator(rate=r) ], end=100)
    pop = simu.population(0)
    # we can see the things are going well
    #Dump(pop)
    # but we need to check it automatically
    for i in range( pop.popSize() ):
      ind = pop.individual(i)
      if ind.sex() == Male:
        # check the second chromosome
        # arrGenotype(ploidy, chrom), no dict parameter for efficiency purpose
        assert ind.arrGenotype(1, 1) == carray('B', [3]*5)
      else: 
        # there is no allele 3 anywhere
        assert ind.arrGenotype().count(3) == 0
       


if __name__ == '__main__':
  unittest.main()
   


