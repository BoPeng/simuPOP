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

  def testAlleleFreq(self):
    'Testing allele frequency (imcomplete)'
    simu = simulator(population(subPop=[500,100,1000], 
      ploidy=2, loci = [1]), randomMating()       )
    init = initByValue(
      value = [[1,1],[1,2],[2,2],[1,1],[1,2],[2,2],[1,2],[1,2],[2,2]],
      indRange = [[0,125], [125,375],[375,500],[500,550],
        [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
    count = stat(
      popSize = 1,
      numOfMale = 1,
      numOfAffected = 1,
      alleleFreq = [0],
      heteroFreq = [0],
      homoFreq = [0],
      expHetero = [0],
      genoFreq = [0],
      Fst = [ 0 ]
    )
    pene = mapPenetrance(locus=0, stage=PreMating,
      penetrance={'1-1':0,'1-2':1,'2-2':0} )
    #simu.apply([init, pene])
    #
    #SaveFstat(simu.population(0), "p1.dat", maxAllele=2)
    # Fst is compaared with result from Fstat.
    #
    #For locus : loc_0_
    #Allele  Capf	 Theta	Smallf	 Relat	Relatc	 Sig_a	 Sig_b	 Sig_w
    #    1	-0.103	 0.102	-0.229	 0.228
    #    2	-0.103	 0.102	-0.229	 0.228
    #  All	-0.103	 0.102	-0.229	 0.228	 0.373	 0.051	-0.102	 0.550
    #
    simu.apply([count])
    # test number of alleles
    pop = population(subPop=[500,1000,1000], ploidy=2, loci = [20])
    InitByFreq(pop, [.2]*5)
    Stat(pop, LD=[[1,2],[2,3,4,5]])


  def testHaplotypeFreq(self):
    'Testing calculation of haplotype frequency (imcomplete)'
    # test haplotype frequency
    pop = population(subPop=[5,10], ploidy=2, loci = [10])
    InitByValue(pop, value=[[1]*10,[2]*10,[3]*10], proportions=[.2,.3,.5])
    Stat(pop, haploFreq=[[0,1],[0,5]], homoFreq=[0,3], heteroFreq=[0,3])

  def testLD(self):
    'Testing calculation of LD (imcomplete)'
    pop = population(subPop=[5,10], ploidy=2, loci = [20]*20)
    InitByFreq(pop, [.2,.3,.5])
    Stat(pop, LD=[[4,5],[100,120]])


if __name__ == '__main__':
  unittest.main()
