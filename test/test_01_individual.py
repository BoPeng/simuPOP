#!/usr/bin/env python
#
# Purpose:
#
# This is a unittest file for individual object
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

class TestIndividual(unittest.TestCase):

  def testIndProperties(self):
    'Testing individual properties'
    if alleleType() != 'binary':
      pop = population(size=100, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        maxAllele=4, alleleNames=['_','A','C','T','G']) 
    else:
      pop = population(size=100, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        alleleNames=['1','2']) 
    #
    pop.individual(99)
    pop.individual(19, 0)
    pop.individual(79, 1)
    self.assertRaises(exceptions.IndexError, pop.individual, 100)
    # this one is allowed.
    pop.individual(20, 0)
    self.assertRaises(exceptions.IndexError, pop.individual, 80, 1)
    self.assertRaises(exceptions.IndexError, pop.individual, 0, 2)
    # 
    ind = pop.individual(10)
    #
    self.assertEqual(ind.ploidy(), 2)
    self.assertEqual(ind.ploidyName(), 'diploid')
    #
    self.assertEqual(ind.numChrom(), 2)
    #
    self.assertEqual(ind.numLoci(0), 5)
    self.assertEqual(ind.numLoci(1), 7)
    self.assertRaises(exceptions.IndexError, ind.numLoci, 2 )
    #
    self.assertEqual(ind.locusPos(10), 12)
    self.assertRaises(exceptions.IndexError, ind.locusPos, 20 )
    self.assertRaises( (exceptions.TypeError, exceptions.OverflowError), ind.locusPos, -1 )
    # more use of arr.. please see test_00_carray
    self.assertEqual(len(ind.arrLociPos()), 12)
    self.assertEqual(ind.arrLociPos().tolist(), [2,3,4,5,6,2,4,6,8,10,12,14])
    #
    self.assertEqual(ind.chromBegin(0), 0)
    self.assertEqual(ind.chromBegin(1), 5)
    self.assertEqual(ind.chromEnd(0), 5)
    self.assertEqual(ind.chromEnd(1), 12)
    self.assertRaises(exceptions.IndexError, ind.chromBegin, 2 )
    self.assertRaises(exceptions.IndexError, ind.chromEnd, 2 )
    #
    self.assertEqual(ind.absLocusIndex(1,5), 10)
    self.assertEqual(ind.locusPos(ind.absLocusIndex(1,2) ), 6)
    self.assertRaises(exceptions.IndexError, ind.absLocusIndex, 2, 5 )
    #
    self.assertEqual(ind.chromLocusPair(10), (1,5) )
    self.assertRaises(exceptions.IndexError, ind.chromLocusPair, 50 )
    #
    if alleleType() == 'binary':
      self.assertEqual(pop.alleleNames(), ('1','2') )
      self.assertEqual(pop.alleleName(0), '1')
      self.assertEqual(pop.alleleName(1), '2')
      # 5 is passed to be function as bool
      self.assertEqual(pop.alleleName(5), '2')
    else:
      self.assertEqual(pop.alleleName(0), '_')
      self.assertEqual(pop.alleleName(1), 'A')
      self.assertEqual(pop.alleleName(2), 'C')
      self.assertEqual(pop.alleleName(3), 'T')
      self.assertEqual(pop.alleleName(4), 'G')
      self.assertRaises(exceptions.IndexError, pop.alleleName, 5)
    # loci name, default, the name will be used by other programs
    # or file format, so we set it to be one based.
    self.assertEqual(pop.locusName(0), 'loc1-1')
    self.assertEqual(pop.locusName(1), 'loc1-2')
    self.assertEqual(pop.locusName(5), 'loc2-1')
    pop = population(loci=[1,2], lociNames=['la','lb','lc'])
    self.assertEqual(pop.locusName(0), 'la')
    self.assertEqual(pop.locusName(1), 'lb')
    self.assertEqual(pop.locusName(2), 'lc')
    self.assertRaises(exceptions.IndexError, pop.locusName, 5)

  def testIndGenotype(self):
    'Testing individual genotype manipulation function'
    if alleleType() != 'binary':
      pop = population(size=100, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        maxAllele=4, alleleNames=['_','A','C','T','G']) 
    else: # binary
      pop = population(size=100, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        alleleNames=['1','2']) 
    ind = pop.individual(0)
    gt = ind.arrGenotype()
    if alleleType() == 'binary':
      gt[:] = [0,1]*12
      self.assertEqual(gt, [0,1]*12)
    else:
      gt[:] = [2,3,4]*8
      self.assertEqual(gt, [2,3,4]*8)
    # ploidy 1
    gt = ind.arrGenotype(1)
    if alleleType() == 'binary':
      self.assertEqual(gt, [0,1]*6)
    else:
      self.assertEqual(gt, [2,3,4]*4)
    # ploidy 1, ch 1
    gt = ind.arrGenotype(1, 1)
    if alleleType() == 'binary':
      self.assertEqual(gt, [1,0,1,0,1,0,1])
    else:
      self.assertEqual(gt, [4,2,3,4,2,3,4])
    self.assertRaises(exceptions.IndexError, ind.arrGenotype, 2)
    self.assertRaises(exceptions.IndexError, ind.arrGenotype, 0, 2)
    if alleleType() == 'binary':
      # layout 0 1 0 1 0 | 1 0 1 0 1 0 1 || 0 1 0 1 0 | 1 0 1 0 1 0 1 
      self.assertEqual( ind.allele(0), 0)
      self.assertEqual( ind.allele(1), 1)
      self.assertEqual( ind.allele(5), 1)
      self.assertEqual( ind.allele(1, 1), 1)
      self.assertEqual( ind.allele(2, 1), 0)
      self.assertEqual( ind.allele(2, 1, 1), 1)
      self.assertEqual( ind.allele(2, 0, 1), 1)
      self.assertRaises(exceptions.IndexError, ind.allele, 24)
      self.assertRaises(exceptions.IndexError, ind.allele, 12, 0)
      self.assertRaises(exceptions.IndexError, ind.allele, 5, 0, 0)
      self.assertRaises(exceptions.IndexError, ind.allele, 0, 2)
      self.assertRaises(exceptions.IndexError, ind.allele, 0, 0, 2)
    else:
      # layout 2 3 4 2 3 | 4 2 3 4 2 3 4 || 2 3 4 2 3 | 4 2 3 4 2 3 4 
      self.assertEqual( ind.allele(0), 2)
      self.assertEqual( ind.allele(1), 3)
      self.assertEqual( ind.allele(5), 4)
      self.assertEqual( ind.allele(1, 1), 3)
      self.assertEqual( ind.allele(2, 1), 4)
      self.assertEqual( ind.allele(2, 1, 1), 3)
      self.assertEqual( ind.allele(2, 0, 1), 3)
      self.assertRaises(exceptions.IndexError, ind.allele, 24)
      self.assertRaises(exceptions.IndexError, ind.allele, 12, 0)
      self.assertRaises(exceptions.IndexError, ind.allele, 5, 0, 0)
      self.assertRaises(exceptions.IndexError, ind.allele, 0, 2)
      self.assertRaises(exceptions.IndexError, ind.allele, 0, 0, 2)
    if alleleType() == 'binary':
      # layout 0 1 0 1 0 | 1 0 1 0 1 0 1 || 0 1 0 1 0 | 1 0 1 0 1 0 1 
      self.assertEqual( ind.alleleChar(0), '1')
      self.assertEqual( ind.alleleChar(1), '2')
      self.assertEqual( ind.alleleChar(5), '2')
      self.assertEqual( ind.alleleChar(1, 1), '2')
      self.assertEqual( ind.alleleChar(2, 1), '1')
      self.assertEqual( ind.alleleChar(2, 1, 1), '2')
      self.assertEqual( ind.alleleChar(2, 0, 1), '2')
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 24)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 12, 0)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 5, 0, 0)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 0, 2)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 0, 0, 2)
    else:
      # layout 2 3 4 2 3 | 4 2 3 4 2 3 4 || 2 3 4 2 3 | 4 2 3 4 2 3 4 
      # ACTG
      self.assertEqual( ind.alleleChar(0), 'C')
      self.assertEqual( ind.alleleChar(1), 'T')
      self.assertEqual( ind.alleleChar(5), 'G')
      self.assertEqual( ind.alleleChar(1, 1), 'T')
      self.assertEqual( ind.alleleChar(2, 1), 'G')
      self.assertEqual( ind.alleleChar(2, 1, 1), 'T')
      self.assertEqual( ind.alleleChar(2, 0, 1), 'T')
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 24)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 12, 0)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 5, 0, 0)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 0, 2)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 0, 0, 2)
    if alleleType() == 'binary':
      # set allele
      # layout 0 1 0 1 0 | 1 0 1 0 1 0 1 || 0 1 0 1 0 | 1 0 1 0 1 0 1 
      ind.setAllele(0,1)
      self.assertEqual( ind.allele(0), 0)
      ind.setAllele(0,5)
      self.assertEqual( ind.allele(5), 0)
      ind.setAllele(0, 1, 1)
      self.assertEqual( ind.allele(1, 1), 0)
      ind.setAllele(1, 2, 1)
      self.assertEqual( ind.allele(2, 1), 1)
      ind.setAllele(0, 2, 1, 1)
      self.assertEqual( ind.allele(2, 1, 1), 0)
      ind.setAllele(0, 2, 0, 1)
      self.assertEqual( ind.allele(2, 0, 1), 0)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 24)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 12, 0)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 5, 0, 0)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 0, 2)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 0, 0, 2)
    else:
      # layout 2 3 4 2 3 | 4 2 3 4 2 3 4 || 2 3 4 2 3 | 4 2 3 4 2 3 4 
      ind.setAllele(1,1)
      self.assertEqual( ind.allele(1), 1)
      ind.setAllele(2,1)
      self.assertEqual( ind.allele(1), 2)
      ind.setAllele(3,5)
      self.assertEqual( ind.allele(5), 3)
      ind.setAllele(1,1,1)
      self.assertEqual( ind.allele(1, 1), 1)
      ind.setAllele(2,2,1)
      self.assertEqual( ind.allele(2, 1), 2)
      ind.setAllele(1,2,1,1)
      self.assertEqual( ind.allele(2, 1, 1), 1)
      ind.setAllele(1,2,0,1)
      self.assertEqual( ind.allele(2, 0, 1), 1)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 24)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 12, 0)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 5, 0, 0)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 0, 2)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 0, 0, 2)
    ind.setTag((1,2))
    self.assertEqual(ind.tag(), (1,2))
    # default to Male
    self.assertEqual(ind.sex(), Male)
    self.assertEqual(ind.sexChar(), 'M')
    ind.setSex(Female)
    self.assertEqual(ind.sex(), Female)
    self.assertEqual(ind.sexChar(), 'F')
    # affectedness
    self.assertEqual(ind.affected(), False)
    self.assertEqual(ind.unaffected(), True)
    ind.setAffected(True)
    self.assertEqual(ind.affected(), True)
    self.assertEqual(ind.unaffected(), False)
    # info
    ind.setInfo(20)
    self.assertEqual(ind.info(), 20)
    
  def testCompare(self):
    'Testing individual comparison'
    pop = population(10, loci=[2])
    self.assertEqual( pop.individual(0) == pop.individual(1), True)
    pop.individual(0).setAllele(1, 0)
    self.assertEqual( pop.individual(0) == pop.individual(1), False)
    pop1 = population(10, loci=[2])
    self.assertEqual( pop.individual(1) == pop1.individual(1), True)
    pop1 = population(10, loci=[3])
    self.assertEqual( pop.individual(1) == pop1.individual(1), False)


if __name__ == '__main__':
  unittest.main()
