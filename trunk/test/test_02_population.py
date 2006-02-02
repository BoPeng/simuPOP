#!/usr/bin/env python
#
# Purpose:
#
# This is a unittest file for population object
#
# Bo Peng (bpeng@rice.edu)
# 
# Last Modified: $Date$
# 
 
from simuPOP import *
import unittest, os, sys, exceptions

class TestPopulation(unittest.TestCase):

  def testNewPopulation(self):
    self.assertRaises(exceptions.ValueError,
      [1,2,4,5].index, 6)
    # no exception '
    if alleleType() == 'binary':
      population(size=100, ploidy=2, loci=[5, 7], subPop=[20, 80],
        lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        maxAllele=1, alleleNames=['_','A','C','T','G'])
    else:
      population(size=100, ploidy=2, loci=[5, 7], subPop=[20, 80],
        lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        maxAllele=4, alleleNames=['_','A','C','T','G'])
    # raise exception when max allele is wrong
    self.assertRaises(exceptions.ValueError,
      population, size=10, maxAllele=MaxAllele+1)
    # raise exception when subpop size does not match
    self.assertRaises(exceptions.ValueError, 
      population, size=1000, ploidy=4, loci=[5, 7]*20, \
      subPop=[20, 80]*20)
    # raise exception when ploidy value is wrong
    self.assertRaises(exceptions.ValueError, 
      population, subPop=[20,20], ploidy=0)
    # raise exceptions with negative values
    self.assertRaises(exceptions.TypeError,
      population, size=-10)
    self.assertRaises(exceptions.TypeError,
      population, subPop=[-10])
    self.assertRaises(exceptions.TypeError,
      population, ploidy=-1)
    # lociDist is depreciated
    self.assertRaises(exceptions.TypeError,
      population, loci=[2], lociDist=[1,2])
    # loci distance in order.
    self.assertRaises(exceptions.ValueError,
      population, loci=[2], lociPos=[3,2])
    # loci distance is in order.
    population(loci=[2,3], lociPos=[1,2,3,5,6])
    # loci distance can be given in another form
    population(loci=[2,3], lociPos=[[1,2],[3,5,6]])
    # but not if number mismatch
    self.assertRaises(exceptions.ValueError,
      population, loci=[2,3], lociPos=[[1,2,3], [8,9]])
    # size mismatch
    self.assertRaises(exceptions.ValueError,
      population, loci=[2], lociPos=[1,2,3])
    self.assertRaises(exceptions.ValueError,
      population, loci=[2], lociPos=[1])
    # loci names
    self.assertRaises(exceptions.ValueError,
      population, loci=[2], lociNames=['1', '2' , '3'])
    # allele names, should only gives warning
    population(10, alleleNames=['-'])
    if alleleType() == 'binary':
      population(10, maxAllele=1, alleleNames=['A', 'B'])
    else:
      population(10, maxAllele=2, alleleNames=['-', 'A', 'B'])
    # sex chrom, 2 does not trigger type error
    p = population(loci=[2], sexChrom=2)
    self.assertEqual( p.sexChrom(), True)
    p = population(loci=[2], sexChrom=False)
    self.assertEqual( p.sexChrom(), False)
    # default population size is 0
    p = population()
    self.assertEqual( p.popSize(), 0) 
    self.assertEqual( p.subPopSize(0), 0)
    self.assertEqual( p.sexChrom(), False)
    # ancestral depth
    p = population(ancestralDepth=-1)
    p = population(ancestralDepth=5)
    # max allele can not be zero
    self.assertRaises(exceptions.ValueError,
      population, maxAllele=0)

  def testPopProperties(self):
    if alleleType() != 'binary':
      pop = population(size=100, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        maxAllele=4, alleleNames=['_','A','C','T','G']) 
    else:
      pop = population(size=100, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        alleleNames=['1','2']) 
    #
    self.assertEqual(pop.popSize(), 100)
    #
    self.assertEqual(pop.subPopSize(0), 20)
    self.assertEqual(pop.subPopSize(1), 80)
    self.assertRaises(exceptions.IndexError, pop.subPopSize, 2 )
    #
    self.assertEqual(pop.subPopSizes(), (20,80) )
    #
    self.assertEqual(pop.subPopBegin(1), 20)
    self.assertRaises(exceptions.IndexError, pop.subPopBegin, 2 )
    self.assertEqual(pop.subPopEnd(0), 20)
    self.assertRaises(exceptions.IndexError, pop.subPopEnd, 2 )
    #
    self.assertEqual(pop.ploidy(), 2)
    self.assertEqual(pop.ploidyName(), 'diploid')
    #
    self.assertEqual(pop.numChrom(), 2)
    #
    self.assertEqual(pop.numLoci(0), 5)
    self.assertEqual(pop.numLoci(1), 7)
    self.assertRaises(exceptions.IndexError, pop.numLoci, 2 )
    #
    self.assertEqual(pop.locusPos(10), 12)
    self.assertRaises(exceptions.IndexError, pop.locusPos, 20 )
    self.assertRaises(exceptions.TypeError, pop.locusPos, -1 )
    # more use of arr.. please see test_00_carray
    self.assertEqual(len(pop.arrLociPos()), 12)
    self.assertEqual(pop.arrLociPos().tolist(), [2,3,4,5,6,2,4,6,8,10,12,14])
    #
    self.assertEqual(pop.chromBegin(0), 0)
    self.assertEqual(pop.chromBegin(1), 5)
    self.assertEqual(pop.chromEnd(0), 5)
    self.assertEqual(pop.chromEnd(1), 12)
    self.assertRaises(exceptions.IndexError, pop.chromBegin, 2 )
    self.assertRaises(exceptions.IndexError, pop.chromEnd, 2 )
    #
    self.assertEqual(pop.absLocusIndex(1,5), 10)
    self.assertEqual(pop.locusPos(pop.absLocusIndex(1,2) ), 6)
    self.assertRaises(exceptions.IndexError, pop.absLocusIndex, 2, 5 )
    #
    self.assertEqual(pop.chromLocusPair(10), (1,5) )
    self.assertRaises(exceptions.IndexError, pop.chromLocusPair, 50 )
    #
    self.assertEqual(pop.numSubPop(), 2)
    #
    self.assertEqual(pop.totNumLoci(), 12)
    #
    self.assertEqual(pop.genoSize(), pop.totNumLoci()*pop.ploidy() )
    # ind, subPop
    self.assertEqual(pop.absIndIndex(1,1), 21)
    self.assertRaises(exceptions.IndexError, pop.absIndIndex, 0, 2 )
    #
    self.assertEqual(pop.subPopIndPair(21), (1,1) )
    self.assertRaises(exceptions.IndexError, pop.subPopIndPair, 200 )
    #
    if alleleType() == 'binary':
      self.assertEqual(pop.alleleNames(), ('1','2') )
      self.assertEqual(pop.alleleName(0), '1')
      self.assertEqual(pop.alleleName(1), '2')
      # 5 is passed to be function as bool
      self.assertEqual(pop.alleleName(5), '2')
    else:
      self.assertEqual(pop.alleleName(0), '-')
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
    
  def testIndProperties(self):
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
    self.assertRaises(exceptions.IndexError, pop.individual, 101)
    self.assertRaises(exceptions.IndexError, pop.individual, 20, 0)
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
    self.assertRaises(exceptions.TypeError, ind.locusPos, -1 )
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
      self.assertEqual(pop.alleleName(0), '-')
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
    if alleleType() != 'binary':
      pop = population(size=1, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        maxAllele=4, alleleNames=['_','A','C','T','G']) 
    else:
      pop = population(size=1, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        alleleNames=['1','2']) 
    ind = pop.individual(0)
    gt = ind.arrGenotype()
    if alleleType() == 'binary':
      gt[:] = [0,1]*6
      self.assertEqual(list(gt), [0,1]*6)
    else:
      gt[:] = [2,3,4]*4
      self.assertEqual(list(gt), [2,3,4]*6)
    # ploidy 1
    gt = ind.arrGenotype(1)
    if alleleType() == 'binary':
      self.assertEqual(list(gt1), [0,1]*3)
    else:
      self.assertEqual(list(gt1), [2,3,4]*3)
    # ploidy 1, ch 1
    gt = ind.arrGenotype(p=1, ch=1)
    if alleleType() == 'binary':
      self.assertEqual(list(gt1), [1,0,1,0,1,0,1])
    else:
      self.assertEqual(list(gt1), [4,2,3,4,2,3,4])
    self.assertRaises(exceptions.IndexError, ind.arrGenotype, 2)
    self.assertRaises(exceptions.IndexError, ind.arrGenotype, 0, 2)
    if alleleType() == 'binary':
      # layout 0 1 0 1 0 | 1 0 1 0 1 0 1 || 0 1 0 1 0 | 1 0 1 0 1 0 1 
      self.assertEqual( ind.allele(0), 0)
      self.assertEqual( ind.allele(1), 1)
      self.assertEqual( ind.allele(5), 1)
      self.assertEqual( ind.allele(1, p=1), 1)
      self.assertEqual( ind.allele(2, p=1), 0)
      self.assertEqual( ind.allele(2, p=1, ch=1), 1)
      self.assertEqual( ind.allele(2, p=0, ch=1), 1)
      self.assertRaises(exceptions.IndexError, ind.allele, 24)
      self.assertRaises(exceptions.IndexError, ind.allele, 12, p=0)
      self.assertRaises(exceptions.IndexError, ind.allele, 5, p=0, ch=1)
      self.assertRaises(exceptions.IndexError, ind.allele, 0, p=2)
      self.assertRaises(exceptions.IndexError, ind.allele, 0, p=0, ch=2)
    else:
      # layout 2 3 4 2 3 | 4 2 3 4 2 3 4 || 2 3 4 2 3 | 4 2 3 4 2 3 4 
      self.assertEqual( ind.allele(0), 2)
      self.assertEqual( ind.allele(1), 3)
      self.assertEqual( ind.allele(5), 4)
      self.assertEqual( ind.allele(1, p=1), 2)
      self.assertEqual( ind.allele(2, p=1), 3)
      self.assertEqual( ind.allele(2, p=1, ch=1), 2)
      self.assertEqual( ind.allele(2, p=0, ch=1), 2)
      self.assertRaises(exceptions.IndexError, ind.allele, 24)
      self.assertRaises(exceptions.IndexError, ind.allele, 12, p=0)
      self.assertRaises(exceptions.IndexError, ind.allele, 5, p=0, ch=1)
      self.assertRaises(exceptions.IndexError, ind.allele, 0, p=2)
      self.assertRaises(exceptions.IndexError, ind.allele, 0, p=0, ch=2)
    if alleleType() == 'binary':
      # layout 0 1 0 1 0 | 1 0 1 0 1 0 1 || 0 1 0 1 0 | 1 0 1 0 1 0 1 
      self.assertEqual( ind.alleleChar(0), '1')
      self.assertEqual( ind.alleleChar(1), '2')
      self.assertEqual( ind.alleleChar(5), '2')
      self.assertEqual( ind.alleleChar(1, p=1), '2')
      self.assertEqual( ind.alleleChar(2, p=1), '1')
      self.assertEqual( ind.alleleChar(2, p=1, ch=1), '2')
      self.assertEqual( ind.alleleChar(2, p=0, ch=1), '2')
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 24)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 12, p=0)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 5, p=0, ch=1)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 0, p=2)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 0, p=0, ch=2)
    else:
      # layout 2 3 4 2 3 | 4 2 3 4 2 3 4 || 2 3 4 2 3 | 4 2 3 4 2 3 4 
      self.assertEqual( ind.alleleChar(0), 'C')
      self.assertEqual( ind.alleleChar(1), 'T')
      self.assertEqual( ind.alleleChar(5), 'G')
      self.assertEqual( ind.alleleChar(1, p=1), 'C')
      self.assertEqual( ind.alleleChar(2, p=1), 'T')
      self.assertEqual( ind.alleleChar(2, p=1, ch=1), 'C')
      self.assertEqual( ind.alleleChar(2, p=0, ch=1), 'C')
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 24)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 12, p=0)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 5, p=0, ch=1)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 0, p=2)
      self.assertRaises(exceptions.IndexError, ind.alleleChar, 0, p=0, ch=2)
    if alleleType() == 'binary':
      # set allele
      # layout 0 1 0 1 0 | 1 0 1 0 1 0 1 || 0 1 0 1 0 | 1 0 1 0 1 0 1 
      ind.setAllele(0,1)
      self.assertEqual( ind.allele(0), 1)
      ind.setAllele(0,5)
      self.assertEqual( ind.allele(5), 0)
      ind.setAllele(0, 1, p=1)
      self.assertEqual( ind.allele(1, p=1), 0)
      ind.setAllele(1, 2, p=1)
      self.assertEqual( ind.allele(2, p=1), 1)
      ind.setAllele(0, 2, p=1, ch=1)
      self.assertEqual( ind.allele(2, p=1, ch=1), 0)
      ind.setAllele(0, 2, p=0, ch=1)
      self.assertEqual( ind.allele(2, p=0, ch=1), 0)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 24)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 12, p=0)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 5, p=0, ch=1)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 0, p=2)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 0, p=0, ch=2)
    else:
      # layout 2 3 4 2 3 | 4 2 3 4 2 3 4 || 2 3 4 2 3 | 4 2 3 4 2 3 4 
      ind.setAllele(1,1)
      self.assertEqual( ind.allele(0), 1)
      ind.setAllele(2,1)
      self.assertEqual( ind.allele(1), 2)
      ind.setAllele(3,5)
      self.assertEqual( ind.allele(5), 3)
      ind.setAllele(1,1,p=1)
      self.assertEqual( ind.allele(1, p=1), 1)
      ind.setAllele(2,2,p=1)
      self.assertEqual( ind.allele(2, p=1), 2)
      ind.setAllele(1,2,p=1,ch=1)
      self.assertEqual( ind.allele(2, p=1, ch=1), 1)
      ind.setAllele(1,2,p=0,ch=1)
      self.assertEqual( ind.allele(2, p=0, ch=1), 1)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 24)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 12, p=0)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 5, p=0, ch=1)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 0, p=2)
      self.assertRaises(exceptions.IndexError, ind.setAllele, 0, 0, p=0, ch=2)
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
    self.assertEqual(ind.Unaffected(), True)
    ind.setAffected(True)
    self.assertEqual(ind.affected(), True)
    self.assertEqual(ind.Unaffected(), False)
    # info
    ind.setInfo(20)
    self.assertEqual(ind.info(), 20)
    
  def testPopManupulation(self):
    self.pop.setIndInfoWithSubPopID()
    
  def testIndividualInfo(self):
    self.assertEqual( self.pop.individual(20).info(), 0)
    self.pop.individual(20).setInfo(1)
    self.assertEqual( self.pop.individual(20).info(), 1)
    
  def testSaveLoadPopulation(self):
    for file in ['a.txt', 'a.xml', 'a.bin']:
      self.pop.savePopulation(file)
      assert os.path.isfile(file)
      pop1 = LoadPopulation(file)
      os.remove(file)
    
  def testArrGenotype(self):
    self.assertEqual( len( self.pop.arrGenotype() ),
      self.pop.genoSize() * self.pop.popSize() ) 
    self.assertEqual( len( self.pop.arrGenotype(0) ),
      self.pop.genoSize() * self.pop.subPopSize(0) ) 
    
  def testIndividual(self):
    self.pop.individual(0)
    
  def testSetIndInfo(self):
    self.pop.setIndInfo(range(0,100))
    
  def testGrp(self):
    self.assertEqual( self.pop.grp(), -1)
    
  def testRep(self):
    self.assertEqual( self.pop.rep(), -1)
    
  def testVars(self):
    self.pop.vars()

  def testSwap(self):
    pop1 = population()
    InitByFreq(self.pop, [.2,.3,.5])
    #Dump(self.pop)
    pop1.swap(self.pop)
    #Dump(pop1)
    #Dump(self.pop)
    #assert(pop1.
    
#  def testAncestry(self):
#    # test ancestry history features
#    pop = population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
#    InitByFreq(pop, [.2,.8])
#    Dump(pop, ancestralPops=True)
#    pop1 = population(subPop=[2,3], loci=[2,3], ancestralDepth=2)
#    InitByFreq(pop1, [.8,.2])
#    Dump(pop1)
#    
#    print popncestralDepth()
#    pop.pushAndDiscard(pop1)
#    print popncestralDepth()
#    Dump(pop1)
#    Dump(pop, ancestralPops=True)
#    
#    tmp = pop.clone()
#    pop.pushAndDiscard(tmp)
#    print popncestralDepth()
#    Dump(pop, ancestralPops=True)
#    
#    # test SavePopulations and LoadPopulations
#    pop = population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
#    InitByFreq(pop, [.2,.8])
#    pop1 = population(subPop=[2,3], loci=[2,3], ancestralDepth=2)
#    InitByFreq(pop1, [.8,.2])
#    SavePopulations([pop,pop1], 'pop01.bin')
#    # pop and pop1 is still alive :-)
#    Dump(pop)
#    Dump(pop1)
#    a = LoadPopulations('pop01.bin')
#    print len(a)
#    Dump(a[0])
#    Dump(a[1])
#    
#    # test split and merge subpopulations
#    pop = population(subPop=[5,6,7])
#    InitByFreq(pop,[.2,.8])
#    Dump(pop)
#    SplitSubPop(pop, 1, [2,4],subPopID=[4,1])
#    Dump(pop)
#    pop = population(subPop=[5,6,7])
#    InitByFreq(pop,[.2,.8])
#    Dump(pop)
#    SplitSubPop(pop, 2, proportions=[.5,.5])
#    Dump(pop)
#    MergeSubPops(pop)
#    Dump(pop)
#    SplitSubPop(pop, 0, proportions=[.2,.3,.5])
#    Dump(pop)
#    MergeSubPops(pop,[0,2])
#    Dump(pop)
#    SplitSubPop(pop, 0, proportions=[.2,.3,.5])
#    Dump(pop)
#    MergeSubPops(pop,[2,0])
#    Dump(pop)
#    SplitSubPop(pop, 3, proportions=[.5,.5], subPopID=[-1,0])
#    Dump(pop)
#    
#    #
#    #
#    ## # testing serialization of shared vars
#    ## pop = population(10)
#    ## InitByFreq(pop, [.2,.4,.4])
#    ## Stat(pop, alleleFreq=[0])
#    ## d = pop.dvars()
#    ## d.a = 1
#    ## d.b = 1.0
#    ## d.c = [1,2,3.5, "a"]
#    ## d.d = {'1':3,'4':6}
#    ## s =  pop.varsAsString()
#    ## print s
#    ## pop.varsFromString(s)
#    ## listVars(pop.dvars())
#    
#    
#    #
#    # test remove loci
#    pop = population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
#    pop1 = population(subPop=[2,3], loci=[2,3], ancestralDepth=2)
#    InitByFreq(pop, [.2,.8])
#    InitByFreq(pop1, [.5,.5])
#    pop.pushAndDiscard(pop1)
#    Dump(pop, ancestralPops=1)
#    pop.removeLoci(remove=[0,1])
#    #pop.removeLoci(keep=[3])
#    Dump(pop, ancestralPops=1)
#    
#    
#    #
#    # load population for testing.
#    #a = LoadPopulation("pop1.bin")
#    # save and load
#    #SaveFstat(a,'a.dat')
#    #p = LoadFstat('a.dat')
#    #Dump(p)
#    
#    turnOnDebug(DBG_UTILITY)
#    # copy of variables
#    pop =  population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
#    InitByFreq(pop, [.2,.8])
#    Stat(pop, LD=[[1,2]])
#    pop.savePopulation('a.bin')
#    pop1 = LoadPopulation('a.bin')
#    
#    pop1 = pop.clone()
#    pop1.dvars()
#    pop1.savePopulation('a.bin')
#    pop.pushAndDiscard(pop1)
#    pop.dvars().a=1
#    Dump(pop, ancestralPops=1)
#    pop.savePopulation('a.bin')
#    
#    pop2 = pop.clone()
#    Dump(pop2, ancestralPops=1)
#    simu = simulator(pop, randomMating())
#    SavePopulation(simu.getPopulation(0), 'a.bin')
#    simu.population(0).dvars()
#    p = simu.getPopulation(0)
#    p.savePopulation('a.bin')
#    p.dvars().clear()
#
#
if __name__ == '__main__':
  unittest.main()
