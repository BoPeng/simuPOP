#!/usr/bin/env python
#
# unittests for mating schemes
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
import unittest, os, sys, random, math

def setGen(pop, off, dad, mom):
  off.setAllele(pop.gen(), 0)
  return True
    
class TestMatingSchemes(unittest.TestCase):

  def testNoMating(self):
    'Testing noMating mating scheme'
    simu = simulator(population(10, loci=[1], ploidy=1),
      noMating())
    # during mating operator will be applied
    simu.step( ops=[ pyOperator(func=setGen, stage=DuringMating) ])
    self.assertEqual( simu.population(0).arrGenotype(),
      [0]*10)
    simu.step( ops=[ pyOperator(func=setGen, stage=DuringMating) ])
    self.assertEqual( simu.population(0).arrGenotype(),
      [1]*10)

  def testBinormialSelection(self):
    'Testing binomialSelection mating scheme (FIXME: imcomplete)'
    simu = simulator(population(10, loci=[1], ploidy=1),
      binomialSelection())
    
  def testSelection(self):
    'Testing selections (FIXME: imcomplete)'
    pass
   
  def testPopSizeChange(self):
    'Testing means to change population size (FIXME: imcomplete)'
    pass
    
  def getFamSize(self, mate, endGen=0, size=1000 ):
    TurnOnDebug(DBG_MATING)
    simu = simulator(population(size, loci=[1]), mate)
    simu.evolve(ops=[], end=endGen)
    return simu.population(0).dvars().famSizes
    
  def testNumOffspring(self):
    'Testing means to control number of offspring (FIXME: check distribution)'
    self.assertEqual( 
      self.getFamSize( binomialSelection(numOffspring=2) ),
      [2]*500)
    # numOffspringFunc
    def nos(gen):
      return gen%2+1
    self.assertEqual( 
      self.getFamSize( binomialSelection(numOffspringFunc=nos), endGen=1),
      [2]*500)
    self.assertEqual( 
      self.getFamSize( binomialSelection(numOffspringFunc=nos), endGen=2),
      [1]*1000)
    # what if each family have different number of offspring?
    # MATE_NumOffspringEachFamily
    def nos(gen):
      return random.randrange(1,4)
    # default mode: the first guess used by all
    cnt = self.getFamSize( randomMating(numOffspringFunc=nos) )
    assert cnt[0]==cnt[1] and cnt[2]==cnt[3]
    #
    cnt = self.getFamSize( randomMating(numOffspringFunc=nos, 
      mode= MATE_NumOffspringEachFamily))
    self.assertEqual( sum(cnt), 1000)
    num = [ cnt.count(i) for i in range(1,4) ]
    # test for uniform?
    mean = sum(num)/3.
    for i in range(3):
      assert num[i] < mean + 50 and num[i] > mean - 50
    #
    # MATE_GeometricDistribution
    cnt = self.getFamSize( randomMating(numOffspring=.3, 
        mode=MATE_GeometricDistribution))
    #print cnt  
    # MATE_BinomialDistribution
    cnt = self.getFamSize( randomMating(numOffspring=.3, 
      maxNumOffspring=5, mode=MATE_BinomialDistribution))
    #print cnt
    # MATE_PoissonDistribution
    cnt = self.getFamSize( randomMating(numOffspring=.3, 
      mode=MATE_PoissonDistribution))
    #print cnt
    
  def testTrajectory(self):
    'Testing trajectory prediction functions'
    sel, Ne, freq, h, selection = 0.5, 100, 0.50, 2, 1
    FreqTrajectorySelSim(sel, Ne, freq, h, selection)
    # the second method, forward, with population expansion
    low, high = 0.5, 0.55
    mutage, grate, N0, sco = 840, 0.01, 1000000, 0.0
    FreqTrajectoryForward(low, high, mutage, grate, N0, sco)

  def testTrajectoryStoch(self):
    'Testing the trajectory obtained from backward binomial sampling'
    TurnOnDebug(DBG_MATING)
    #TurnOnDebug(DBG_DEVEL)
    # fitness
    #   AA     Aa      aa
    #    1     1+s1    1+s2
    # constant population size
    def NtFunc(gen):
      return 10000
    # neutral
    T, freq, s1, s2 = 100000, 0.3, 0., 0.
    path = FreqTrajectoryStoch(NtFunc, T, freq, s1, s2)
    # advantageous allele, s2>s1>0 
    T, freq, s1, s2 = 10000, 0.3, 0., 0.01
    path = FreqTrajectoryStoch(NtFunc, T, freq, s1, s2)
    # overdominance, s1 > s2 > 0
    T, freq, s1, s2 = 10000, 0.3, 0.02, 0.
    path = FreqTrajectoryStoch(NtFunc, T, freq, s1, s2)
    # with week purifying selection (additive)
    T, freq, s1, s2 = 10000, 0.3, -0.0001, -0.0002
    path = FreqTrajectoryStoch(NtFunc, T, freq, s1, s2)
    #
    # population growth
    def NtFunc(gen):
      return 1000+10000*math.exp(-0.001*(10000-gen))
    # neutral
    T, freq, s1, s2 = 100000, 0.3, 0., 0.
    path = FreqTrajectoryStoch(NtFunc, T, freq, s1, s2)
    # advantageous allele, s2>s1>0 
    T, freq, s1, s2 = 10000, 0.3, 0., 0.01
    path = FreqTrajectoryStoch(NtFunc, T, freq, s1, s2)
    # overdominance, s1 > s2 > 0
    T, freq, s1, s2 = 10000, 0.3, 0.02, 0.
    path = FreqTrajectoryStoch(NtFunc, T, freq, s1, s2)
    # with week purifying selection (additive)
    T, freq, s1, s2 = 10000, 0.3, -0.0001, -0.0002
    path = FreqTrajectoryStoch(NtFunc, T, freq, s1, s2)

  def checkRoot(self):
    'Testing the algorithm'
    for s1 in [ x/200. - 0.9 for x in range(400) ]:
      for s2 in [ x/200. - 0.9 for x in range(400) ]:
        for x in [ x/100. for x in range(0,100)]:
          a = s2*x-2*s1*x-s2+s1
          if a == 0:
            continue
          b = 2*s1*x - 1 - s1
          c = x
          b4ac = b*b-4*a*c
          if b4ac < 0:
            raise exceptions.ValueError("b2ac<0")
          y1 = (-b+math.sqrt(b4ac))/(2*a)
          y2 = (-b-math.sqrt(b4ac))/(2*a)
          if (y1 < 0 or y1 > 1) and (y2 < 0 or y2 > 1):
            print s1,s2,x,a,b,c,y1,y2
            #raise exceptions.ValueError("no valid solution")
          if (y1 >= 0 and y1 <= 1) and (y2 >= 0 and y2 <= 1):
            print "over", s1,s2,x,a,b,c,y1,y2
            #raise exceptions.ValueError("no valid solution")

  def testControlledMating(self):
    'Testing controlled mating'
    # parameters
    sel, Ne, freq, h, selection = 0.5, 100, 0.50, 2, 1
    # planned trajectory
    freq = FreqTrajectoryStochSel(sel, Ne, freq, h, selection)
    # staring from when?
    burnin = 100
    mutAge = len(freq)
    # trajectory function
    # 0 ...., 100, 101, .... 100+mutAge
    #              x         freq
    def freqRange(gen):
      if gen <= burnin:
        # whatever
        return [0,1]
      expected = freq[gen-1-burnin]
      return [expected, expected + 0.05]
    #
    # turn On debug
    #TurnOnDebug(DBG_MATING)
    simu = simulator( population(Ne, loci=[1], ploidy=2), 
      controlledMating( matingScheme=randomMating(), 
        locus=0, allele=StartingAllele+1, freqFunc=freqRange ) 
      )
    print "Simulator created"
    simu.evolve( 
      preOps=[
        initByValue([StartingAllele])
        ],
      ops=[
        pointMutator(atLoci=[0], 
          toAllele=StartingAllele+1, 
          inds = [0],
          at = [burnin+1],
          stage = PreMating),
        stat(alleleFreq=[0]),
        pyEval(r'"%%f\n"%%alleleFreq[0][%d]'%StartingAllele)
      ], 
      end=burnin+mutAge
    )
      
    

if __name__ == '__main__':
  unittest.main()
  sys.exit(0)


