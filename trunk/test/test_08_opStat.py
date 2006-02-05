#!/usr/bin/env python
#
# testing statistics calculation
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

class TestStat(unittest.TestCase):

  def testPopSize(self):
    'Testing calculation of population (subPopulation) size'
    pop = population(subPop=[200,800])
    Stat(pop, popSize=1)
    self.assertEqual(pop.dvars().numSubPop, 2)
    self.assertEqual(pop.dvars().popSize, 1000)
    self.assertEqual(pop.dvars(0).popSize, 200)
    self.assertEqual(pop.dvars(1).popSize, 800)
    self.assertEqual(pop.dvars().subPopSize, [200,800])

  def testNumOfMale(self):
    'Tsting counting number of male'
    pop = population(subPop=[200, 800])
    for i in range(100):
      pop.individual(i,0).setSex(Male)
      pop.individual(i,1).setSex(Male)
    for i in range(100,200):
      pop.individual(i,0).setSex(Female)
    for i in range(100,800):
      pop.individual(i,1).setSex(Female)
    Stat(pop, numOfMale=1)
    self.assertEqual(pop.dvars().numOfMale, 200)
    self.assertEqual(pop.dvars().numOfFemale, 800)
    self.assertEqual(pop.dvars(0).numOfMale, 100)
    self.assertEqual(pop.dvars(0).numOfFemale, 100)
    self.assertEqual(pop.dvars(1).numOfMale, 100)
    self.assertEqual(pop.dvars(1).numOfFemale, 700)
      
  def testNumOfAffected(self):
    'Tsting counting number of affected individuals'
    pop = population(subPop=[200, 800])
    for i in range(100):
      pop.individual(i,0).setAffected(True)
      pop.individual(i,1).setAffected(True)
    for i in range(100,200):
      pop.individual(i,0).setAffected(False)
    for i in range(100,800):
      pop.individual(i,1).setAffected(False)
    Stat(pop, numOfAffected=1)
    self.assertEqual(pop.dvars().numOfAffected, 200)
    self.assertEqual(pop.dvars().numOfUnaffected, 800)
    self.assertEqual(pop.dvars(0).numOfAffected, 100)
    self.assertEqual(pop.dvars(0).numOfUnaffected, 100)
    self.assertEqual(pop.dvars(1).numOfAffected, 100)
    self.assertEqual(pop.dvars(1).numOfUnaffected, 700)
    self.assertEqual(pop.dvars().propOfAffected, 0.2)
    self.assertEqual(pop.dvars().propOfUnaffected, 0.8)
    self.assertEqual(pop.dvars(0).propOfAffected, 0.5)
    self.assertEqual(pop.dvars(0).propOfUnaffected, 0.5)
    self.assertEqual(pop.dvars(1).propOfAffected, 1/8.)
    self.assertEqual(pop.dvars(1).propOfUnaffected, 7/8.)
    
  def testAlleleFreq(self):
    'Testing calculation of allele frequency and number of alleles'
    pop = population(subPop=[500,100,1000], 
      ploidy=2, loci = [1])
    InitByValue(pop, 
      value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
      indRange = [[0,125], [125,375],[375,500],[500,550],
        [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
    Stat(pop, alleleFreq=[0], numOfAlleles=[0])
    self.assertEqual(pop.dvars().alleleFreq[0], [1230./3200, 1970./3200])
    self.assertEqual(pop.dvars().alleleNum[0], [1230, 1970])
    self.assertEqual(pop.dvars(0).alleleFreq[0], [.5, .5])
    self.assertEqual(pop.dvars(0).alleleNum[0], [500, 500])
    self.assertEqual(pop.dvars(1).alleleFreq[0], [13./20, 7./20])
    self.assertEqual(pop.dvars(1).alleleNum[0], [130, 70])
    self.assertEqual(pop.dvars(2).alleleFreq[0], [0.3, 0.7])
    self.assertEqual(pop.dvars(2).alleleNum[0], [600, 1400])
    #
    self.assertEqual(pop.dvars().numOfAlleles[0], 2)
    self.assertEqual(pop.dvars(0).numOfAlleles[0], 2)
    self.assertEqual(pop.dvars(1).numOfAlleles[0], 2)

  def testNumOfAlleles(self):
    'Testing calculation of number of alleles'
    pop = population(subPop=[5000,15000], loci=[1])
    InitByFreq(pop, [.2, 0, .5, .3])
    Stat(pop, numOfAlleles=[0])
    if alleleType() == 'binary':
      self.assertEqual(pop.dvars().numOfAlleles[0], 2)
      self.assertEqual(pop.dvars(0).numOfAlleles[0], 2)
      self.assertEqual(pop.dvars(1).numOfAlleles[0], 2)
    else:
      self.assertEqual(pop.dvars().numOfAlleles[0], 3)
      self.assertEqual(pop.dvars(0).numOfAlleles[0], 3)
      self.assertEqual(pop.dvars(1).numOfAlleles[0], 3)
    
  def testHeteroFreq(self):
    'Testing counting of heterozygote frequency'
    pop = population(subPop=[500,100,1000], 
      ploidy=2, loci = [1])
    if alleleType() == 'binary':
      InitByValue(pop, 
        value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
        indRange = [[0,125], [125,375],[375,500],[500,550],
          [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
      Stat(pop, heteroFreq=[0])
      self.assertEqual(pop.dvars().heteroNum[0][0], 880)
      self.assertEqual(pop.dvars().heteroNum[0][1], 880)
      self.assertEqual(pop.dvars().heteroFreq[0][0], 880./1600)
      self.assertEqual(pop.dvars().heteroFreq[0][1], 880./1600)
      self.assertEqual(pop.dvars(0).heteroNum[0][0], 250)
      self.assertEqual(pop.dvars(0).heteroNum[0][1], 250)
      self.assertEqual(pop.dvars(0).heteroFreq[0][0], .5)
      self.assertEqual(pop.dvars(0).heteroFreq[0][1], .5)
      self.assertEqual(pop.dvars(1).heteroNum[0][0], 30)
      self.assertEqual(pop.dvars(1).heteroNum[0][1], 30)
      self.assertEqual(pop.dvars(1).heteroFreq[0][0], 0.3)
      self.assertEqual(pop.dvars(1).heteroFreq[0][1], 0.3)
      self.assertEqual(pop.dvars(2).heteroNum[0][0], 600)
      self.assertEqual(pop.dvars(2).heteroNum[0][1], 600)
      self.assertEqual(pop.dvars(2).heteroFreq[0][0], 0.6)
      self.assertEqual(pop.dvars(2).heteroFreq[0][1], 0.6)
    else:
      InitByValue(pop, 
        value = [[1,1],[1,2],[2,3],[1,1],[3,2],[2,2],[1,2],[3,2],[2,2]],
        indRange = [[0,125], [125,375],[375,500],[500,550],
          [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
      Stat(pop, heteroFreq=[0])
      self.assertEqual(pop.dvars().heteroNum[0][0], 1005)
      self.assertEqual(pop.dvars().heteroNum[0][1], 350)
      self.assertEqual(pop.dvars().heteroNum[0][2], 1005)
      self.assertEqual(pop.dvars().heteroNum[0][3], 655)
      self.assertEqual(pop.dvars().heteroFreq[0][0], 1005./1600)
      self.assertEqual(pop.dvars().heteroFreq[0][1], 350./1600)
      self.assertEqual(pop.dvars().heteroFreq[0][2], 1005./1600)
      self.assertEqual(pop.dvars().heteroFreq[0][3], 655./1600)
      self.assertEqual(pop.dvars(0).heteroNum[0][0], 375)
      self.assertEqual(pop.dvars(0).heteroNum[0][1], 250)
      self.assertEqual(pop.dvars(0).heteroNum[0][2], 375)
      self.assertEqual(pop.dvars(0).heteroNum[0][3], 125)
      self.assertEqual(pop.dvars(0).heteroFreq[0][0], 375/500.)
      self.assertEqual(pop.dvars(0).heteroFreq[0][1], .5)
      self.assertEqual(pop.dvars(0).heteroFreq[0][2], 375/500.)
      self.assertEqual(pop.dvars(0).heteroFreq[0][3], 125/500.)
      self.assertEqual(pop.dvars(1).heteroNum[0][0], 30)
      self.assertEqual(pop.dvars(1).heteroNum[0][1], 0)
      self.assertEqual(pop.dvars(1).heteroNum[0][2], 30)
      self.assertEqual(pop.dvars(1).heteroNum[0][3], 30)
      self.assertEqual(pop.dvars(1).heteroFreq[0][0], 0.3)
      self.assertEqual(pop.dvars(1).heteroFreq[0][1], 0)
      self.assertEqual(pop.dvars(1).heteroFreq[0][2], 0.3)
      self.assertEqual(pop.dvars(1).heteroFreq[0][3], 0.3)
      self.assertEqual(pop.dvars(2).heteroNum[0][0], 600)
      self.assertEqual(pop.dvars(2).heteroNum[0][1], 100)
      self.assertEqual(pop.dvars(2).heteroNum[0][2], 600)
      self.assertEqual(pop.dvars(2).heteroNum[0][3], 500)
      self.assertEqual(pop.dvars(2).heteroFreq[0][0], 0.6)
      self.assertEqual(pop.dvars(2).heteroFreq[0][1], 0.1)
      self.assertEqual(pop.dvars(2).heteroFreq[0][2], 0.6)
      self.assertEqual(pop.dvars(2).heteroFreq[0][3], 0.5)
    
  def testExpHetero(self):
    'Testing expected heterozygosity 1-sum p_i^2'
    pop = population(subPop=[500,100,1000], 
      ploidy=2, loci = [1])
    InitByValue(pop, 
      value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
      indRange = [[0,125], [125,375],[375,500],[500,550],
        [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
    Stat(pop, expHetero=[0])
    #
    assert abs(pop.dvars().expHetero[0] - (1-(123./320)**2-(197./320)**2)) < 0.00001
    assert abs(pop.dvars(0).expHetero[0] - (1-0.5**2-0.5**2)) < 0.00001
    assert abs(pop.dvars(1).expHetero[0] - (1-(13./20)**2-(7./20)**2)) < 0.00001
    assert abs(pop.dvars(2).expHetero[0] - (1-0.3**2-0.7**2)) < 0.00001

  def testGenoFreq(self):
    'Testing the counting of genotype frequency'
    pop = population(subPop=[500,100,1000], 
      ploidy=2, loci = [1])
    InitByValue(pop, 
      value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
      indRange = [[0,125], [125,375],[375,500],[500,550],
        [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
    Stat(pop, genoFreq=[0])
    self.assertEqual(pop.dvars().genoNum[0][0][0], 175)
    self.assertEqual(pop.dvars().genoNum[0][0][1], 880)
    self.assertEqual(pop.dvars().genoNum[0][1][1], 545)
    self.assertEqual(pop.dvars().genoFreq[0][0][0], 175./1600)
    self.assertEqual(pop.dvars().genoFreq[0][0][1], 880./1600)
    self.assertEqual(pop.dvars().genoFreq[0][1][1], 545./1600)
    #
    self.assertEqual(pop.dvars(0).genoNum[0][0][0], 125)
    self.assertEqual(pop.dvars(0).genoNum[0][0][1], 250)
    self.assertEqual(pop.dvars(0).genoNum[0][1][1], 125)
    self.assertEqual(pop.dvars(0).genoFreq[0][0][0], .25)
    self.assertEqual(pop.dvars(0).genoFreq[0][0][1], .5)
    self.assertEqual(pop.dvars(0).genoFreq[0][1][1], .25)
    #
    self.assertEqual(pop.dvars(1).genoNum[0][0][0], 50)
    self.assertEqual(pop.dvars(1).genoNum[0][0][1], 30)
    self.assertEqual(pop.dvars(1).genoNum[0][1][1], 20)
    self.assertEqual(pop.dvars(1).genoFreq[0][0][0], 0.5)
    self.assertEqual(pop.dvars(1).genoFreq[0][0][1], 0.3)
    self.assertEqual(pop.dvars(1).genoFreq[0][1][1], 0.2)
    # key does not exist
    self.assertEqual(pop.dvars(2).genoNum[0][0].setdefault(0,0), 0)
    self.assertEqual(pop.dvars(2).genoNum[0][0][1], 600)
    self.assertEqual(pop.dvars(2).genoNum[0][1][1], 400)
    # key does not exist
    self.assertEqual(pop.dvars(2).genoFreq[0][0].setdefault(0,0), 0)
    self.assertEqual(pop.dvars(2).genoFreq[0][0][1], 0.6)
    self.assertEqual(pop.dvars(2).genoFreq[0][1][1], 0.4)

  def testFst(self):
    'Testing calculation of Fst value'  
    pop = population(subPop=[500,100,1000], 
      ploidy=2, loci = [1])
    if alleleType() == 'binary':
      InitByValue(pop, 
        value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
        indRange = [[0,125], [125,375],[375,500],[500,550],
          [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
    else:
      InitByValue(pop, 
        value = [[1,1],[1,2],[2,2],[1,1],[1,2],[2,2],[1,2],[1,2],[2,2]],
        indRange = [[0,125], [125,375],[375,500],[500,550],
          [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
    #SaveFstat(simu.population(0), "p1.dat", maxAllele=2)
    # Fst is compaared with result from Fstat.
    #
    #For locus : loc_0_
    #Allele  Capf	 Theta	Smallf	 Relat	Relatc	 Sig_a	 Sig_b	 Sig_w
    #    1	-0.103	 0.102	-0.229	 0.228
    #    2	-0.103	 0.102	-0.229	 0.228
    #  All	-0.103	 0.102	-0.229	 0.228	 0.373	 0.051	-0.102	 0.550
    Stat(pop, Fst=[0])
    assert abs(pop.dvars().Fis[0] - (-0.229)) < 0.01
    assert abs(pop.dvars().Fit[0] - (-0.103)) < 0.01
    assert abs(pop.dvars().Fst[0] - (0.102)) < 0.01

  def testHaploFreq(self):
    'Testing calculation of haplotype frequency'
    # test haplotype frequency
    pop = population(subPop=[5000,1000], ploidy=2, loci = [10])
    if alleleType() == 'binary':
      InitByValue(pop, value=[[0]*10,[1]*10], proportions=[.3,.7])
      Stat(pop, haploFreq=[[0,1,5],[2,5]])
      assert abs(pop.dvars().haploFreq['0-1-5']['0-0-0'] - 0.3) < 0.05
      assert abs(pop.dvars().haploFreq['0-1-5']['1-1-1'] - 0.7) < 0.05
      assert abs(pop.dvars().haploFreq['2-5']['0-0'] - 0.3) < 0.05
      assert abs(pop.dvars().haploFreq['2-5']['1-1'] - 0.7) < 0.05
    else:
      InitByValue(pop, value=[[1]*10,[2]*10,[3]*10], 
        proportions=[.2,.3,.5])
      Stat(pop, haploFreq=[[0,1,5],[2,5]])
      assert abs(pop.dvars().haploFreq['0-1-5']['1-1-1'] - 0.2) < 0.05
      assert abs(pop.dvars().haploFreq['0-1-5']['2-2-2'] - 0.3) < 0.05
      assert abs(pop.dvars().haploFreq['0-1-5']['3-3-3'] - 0.5) < 0.05
      assert abs(pop.dvars().haploFreq['2-5']['1-1'] - 0.2) < 0.05
      assert abs(pop.dvars().haploFreq['2-5']['2-2'] - 0.3) < 0.05
      assert abs(pop.dvars().haploFreq['2-5']['3-3'] - 0.5) < 0.05

  def testLD(self):
    'Testing calculation of LD (imcomplete)'
    pop = population(subPop=[500,100,1000], 
      ploidy=2, loci = [4,5])
    InitByFreq(pop, [.2, .3, .5] )
    if alleleType() == 'binary':
      # [4,5,0,1] should not make any difference
      Stat(pop, LD=[[1,2],[4,5,0,1],[2,5]])
    else:
      Stat(pop, LD=[[1,2],[4,5,1,2],[2,5]])
    #
    #assert pop.dvars().LD

if __name__ == '__main__':
  unittest.main()
