#!/usr/bin/env python
#
# Purpose:
#
# This is a unittest file for population object
#
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision$
# $LastChangedDate$
# 

import simuOpt
#simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

class TestPopulation(unittest.TestCase):

  def assertGenotype(self, pop, subPop, genotype):
    'Assert if the genotype of subPop of pop is genotype '
    gt = list(pop.arrGenotype(subPop))
    gt.sort()
    if alleleType() == 'binary':
      self.assertEqual(gt, [x>0 for x in genotype])
    else:
      self.assertEqual(gt, genotype)

  def testNewPopulation(self):
    'Testing the creation of populations'
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
    # MaxAllele may be too big.
    self.assertRaises( (exceptions.TypeError, exceptions.ValueError, exceptions.OverflowError),
      population, size=10, maxAllele=MaxAllele+1)
    # raise exception when subpop size does not match
    self.assertRaises(exceptions.ValueError, 
      population, size=1000, ploidy=4, loci=[5, 7]*20, \
      subPop=[20, 80]*20)
    # raise exception when ploidy value is wrong
    self.assertRaises(exceptions.ValueError, 
      population, subPop=[20,20], ploidy=0)
    # raise exceptions with negative values
    self.assertRaises((exceptions.TypeError, exceptions.OverflowError),
      population, size=-10)
    self.assertRaises((exceptions.TypeError, exceptions.OverflowError),
      population, subPop=[-10])
    self.assertRaises((exceptions.TypeError, exceptions.OverflowError), 
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
    population(10, alleleNames=['_'])
    population(10, maxAllele=1, alleleNames=['A', 'B'])
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

  def testGenoStructure(self):
    'Testing geno structure related functions'
    if alleleType() != 'binary':
      pop = population(size=100, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        maxAllele=3, alleleNames=['A','C','T','G']) 
    else:
      pop = population(size=100, ploidy=2, loci=[5, 7], 
        subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        alleleNames=['1','2']) 
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
    self.assertRaises((exceptions.TypeError, exceptions.OverflowError), pop.locusPos, -1 )
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
    self.assertEqual(pop.totNumLoci(), 12)
    #
    self.assertEqual(pop.genoSize(), pop.totNumLoci()*pop.ploidy() )
    #
    if alleleType() == 'binary':
      self.assertEqual(pop.alleleNames(), ('1','2') )
      self.assertEqual(pop.alleleName(0), '1')
      self.assertEqual(pop.alleleName(1), '2')
      # 5 is passed to be function as bool
      self.assertEqual(pop.alleleName(5), '2')
    else:
      self.assertEqual(pop.alleleName(0), 'A')
      self.assertEqual(pop.alleleName(1), 'C')
      self.assertEqual(pop.alleleName(2), 'T')
      self.assertEqual(pop.alleleName(3), 'G')
      self.assertRaises(exceptions.IndexError, pop.alleleName, 4)
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
    
  def testPopProperties(self):
    'Testing population properties'
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
    self.assertEqual(pop.numSubPop(), 2)
    #
    # ind, subPop
    self.assertEqual(pop.absIndIndex(1,1), 21)
    self.assertRaises(exceptions.IndexError, pop.absIndIndex, 0, 2 )
    #
    self.assertEqual(pop.subPopIndPair(21), (1,1) )
    self.assertRaises(exceptions.IndexError, pop.subPopIndPair, 200 )
   
  def testSetSubPopStru(self):
    'Testing function setSubPopStru'
    pop = population(size=1, loci=[1])
    # pop1 is only a reference to pop
    pop1 = pop
    pop.individual(0).setAllele(1,0)
    # the genotype of pop1 is also changed
    self.assertEqual( pop1.individual(0).allele(0), 1)
    pop2 = pop.clone()
    pop.individual(0).setAllele(0,0)
    # pop2 does not change with pop
    self.assertEqual( pop2.individual(0).allele(0), 1)
    #
    pop = population(size=10)
    self.assertEqual( pop.subPopSizes(), (10,) )
    pop.setSubPopStru(newSubPopSizes=[2,8], allowPopSizeChange=False)
    self.assertEqual( pop.subPopSizes(), (2,8) )
    # can set empty subpop
    pop.setSubPopStru(newSubPopSizes=[0,0,1,0,2,7,0], allowPopSizeChange=False)
    self.assertEqual( pop.subPopSizes(), (0,0,1,0,2,7,0) )
    # by default, can not change population size
    self.assertRaises(exceptions.ValueError, 
      pop.setSubPopStru, newSubPopSizes=[10,20], allowPopSizeChange=False)
    # can change population size if allow... is set to True
    pop.setSubPopStru(newSubPopSizes=[10, 20], allowPopSizeChange=True)
    self.assertEqual( pop.subPopSizes(), (10, 20) )

  def testPopSwap(self):
    'Testing population swap'
    pop = population(10, loci=[2])
    pop1 = population(5, loci=[3])
    InitByFreq(pop, [.2,.3,.5])
    InitByFreq(pop1, [.2,.3,.5])
    popa = pop.clone()
    pop1a = pop1.clone()
    pop1.swap(pop)
    self.assertEqual( pop, pop1a)
    self.assertEqual( pop1, popa)

  def testSplitSubPop(self):
    'Testing function splitSubPop'
    pop = population(subPop=[5,6,7], ploidy=1, loci=[1])
    arr = pop.arrGenotype()
    # mapArr is separate
    mapArr = list(arr)
    mapArr = range(pop.popSize())
    arr[:] = mapArr
    # split, the additional subpop will be put at the end
    # member function form
    pop.splitSubPop(1, [2,4])
    self.assertEqual(pop.subPopSizes(), (5,2,7,4) )
    # check if 7,8,9,10 is moved to subpopulation 3.
    # underlying genotype will *not* be sorted
    self.assertGenotype(pop, 0, [0,1,2,3,4])
    self.assertGenotype(pop, 1, [5,6])
    self.assertGenotype(pop, 2, [11,12,13,14,15,16,17])
    self.assertGenotype(pop, 3, [7,8,9,10])
    #, 5,6, 11,12,13,14,15,16,17, 7,8,9,10])
    #
    # recover population
    pop.setSubPopStru([5,6,7], False)
    pop.arrGenotype()[:] = range(pop.popSize())
    # function form
    SplitSubPop(pop, 1, [2,4], subPopID=[4,1], randomize=False)
    self.assertEqual(pop.subPopSizes(), (5,4,7,0,2))
    self.assertGenotype(pop, 0, [0,1,2,3,4])
    self.assertGenotype(pop, 1, [7,8,9,10])
    self.assertGenotype(pop, 2, [11,12,13,14,15,16,17])
    self.assertGenotype(pop, 3, [])
    self.assertGenotype(pop, 4, [5,6])
    # given wrong split size?
    self.assertRaises(exceptions.ValueError,
      SplitSubPop, pop, 1, [2,4], subPopID=[4,1])
    # if given subPopID is already used?
    print "A warning should be issued"
    SplitSubPop(pop, 0, [2,3], subPopID=[2,3])
    
  def testSplitSubPopByProportion(self):
    'Testing function splitSubPopByProportion'
    # split by proportion -------- 
    pop = population(subPop=[5,6,7], ploidy=1, loci=[1])
    arr = pop.arrGenotype()
    # mapArr is separate
    mapArr = list(arr)
    mapArr = range(pop.popSize())
    arr[:] = mapArr
    # step 1, split, the additional subpop will be put at the end
    pop.splitSubPopByProportion(1, [2/5.,3/5.])
    self.assertEqual(pop.subPopSizes(), (5,2,7,4) )
    # check if 7,8,9,10 is moved to subpopulation 3.
    # underlying genotype will *not* be sorted
    self.assertGenotype(pop, 0, [0,1,2,3,4])
    self.assertGenotype(pop, 1, [5,6])
    self.assertGenotype(pop, 2, [11,12,13,14,15,16,17])
    self.assertGenotype(pop, 3, [7,8,9,10])
    #
    # recover population
    pop.setSubPopStru([5,6,7], False)
    pop.arrGenotype()[:] = range(pop.popSize())
    SplitSubPop(pop, 1, proportions=[2/5.,3/5.], subPopID=[4,1], randomize=False)
    self.assertEqual(pop.subPopSizes(), (5,4,7,0,2))
    self.assertGenotype(pop, 0, [0,1,2,3,4])
    self.assertGenotype(pop, 1, [7,8,9,10])
    self.assertGenotype(pop, 2, [11,12,13,14,15,16,17])
    self.assertGenotype(pop, 3, [])
    self.assertGenotype(pop, 4, [5,6])
    # proportion does not add up to one?
    self.assertRaises(exceptions.ValueError,
      SplitSubPop, pop, 1, proportions=[2/3.,2/3.], subPopID=[4,1])
    # if given subPopID is already used?
    print "A warning should be issued"
    SplitSubPop(pop, 0, proportions=[2/5.,3/5.], subPopID=[2,3])
    #
    # split by proportion
   
  def testRemoveSubPops(self):
    'Testing function removeEmptySubPops, removeSubPops'
    pop = population(subPop=[0,1,0,2,3,0])
    self.assertEqual( pop.numSubPop(), 6)
    pop.removeEmptySubPops()
    self.assertEqual( pop.numSubPop(), 3)
    # remove subpop
    pop = population(subPop=[0,1,0,2,3,0], ploidy=1, loci=[1])
    arr = pop.arrGenotype()
    # mapArr is separate
    arr[:] = range(pop.popSize())
    pop.removeSubPops([1,2])
    self.assertEqual( pop.subPopSizes(), (0,2,3,0))
    # subpop will be shifted
    self.assertGenotype(pop, 0, [])
    self.assertGenotype(pop, 1, [1,2])
    self.assertGenotype(pop, 2, [3,4,5])
    self.assertGenotype(pop, 3, [])
    #
    pop.removeSubPops([2], shiftSubPopID=False)
    self.assertEqual( pop.subPopSizes(), (0,2,0,0))
    # subpop will be shifted
    self.assertGenotype(pop, 0, [])
    self.assertGenotype(pop, 1, [1,2])
    self.assertGenotype(pop, 2, [])
    self.assertGenotype(pop, 3, [])
    # should give warning
    print "A warning should be issued"
    pop.removeSubPops([8])
  
  def testRemoveIndividuals(self):
    'Testing function removeIndividuals'
    pop = population(subPop=[0,1,0,2,3,0], ploidy=1, loci=[1])
    arr = pop.arrGenotype()
    arr[:] = range(pop.popSize())
    #
    pop.removeIndividuals([2])
    self.assertEqual( pop.subPopSizes(), (0,1,0,1,3,0))
    # subpop will be shifted
    self.assertGenotype(pop, 0, [])
    self.assertGenotype(pop, 1, [0])
    self.assertGenotype(pop, 2, [])
    self.assertGenotype(pop, 3, [1])
    self.assertGenotype(pop, 4, [3,4,5])
    #
    pop.removeIndividuals([1], removeEmptySubPops=True)
    self.assertEqual( pop.subPopSizes(), (1,3))
    # subpop will be shifted
    self.assertGenotype(pop, 0, [0])
    self.assertGenotype(pop, 1, [3,4,5])
  
  def testMergeSubPops(self):
    'Testing function mergeSubPops'
    pop = population(subPop=[0,1,0,2,3,0], ploidy=1, loci=[1])
    arr = pop.arrGenotype()
    arr[:] = range(pop.popSize())
    #
    pop.mergeSubPops([1,2,4])
    self.assertEqual( pop.subPopSizes(), (0,4,0,2,0,0))
    # subpop will be shifted
    self.assertGenotype(pop, 0, [])
    self.assertGenotype(pop, 1, [0,3,4,5])
    self.assertGenotype(pop, 2, [])
    self.assertGenotype(pop, 3, [1,2])
    self.assertGenotype(pop, 4, [])
    self.assertGenotype(pop, 5, [])
    #
    pop.mergeSubPops([1,2,3], removeEmptySubPops=True)
    self.assertEqual( pop.subPopSizes(), (6,))
    # subpop will be shifted
    self.assertGenotype(pop, 0, [0,1,2,3,4,5])

  def testReorderSubPops(self):
    'Testing function reorderSubPops'
    pop = population(subPop=[1,2,3,4], ploidy=1, loci=[1])
    arr = pop.arrGenotype()
    arr[:] = range(pop.popSize())
    #
    pop.reorderSubPops(order=[1,3,0,2])
    self.assertEqual( pop.subPopSizes(), (2,4,1,3))
    # subpop will be shifted
    self.assertGenotype(pop, 0, [1,2])
    self.assertGenotype(pop, 1, [6,7,8,9])
    self.assertGenotype(pop, 2, [0])
    self.assertGenotype(pop, 3, [3,4,5])
    # by rank
    pop = population(subPop=[1,2,3,4], ploidy=1, loci=[1])
    arr = pop.arrGenotype()
    arr[:] = range(pop.popSize())
    #
    pop.reorderSubPops(rank=[1,3,0,2])
    self.assertEqual( pop.subPopSizes(), (3,1,4,2))
    # subpop will be shifted
    self.assertGenotype(pop, 0, [3,4,5])
    self.assertGenotype(pop, 1, [0])
    self.assertGenotype(pop, 2, [6,7,8,9])
    self.assertGenotype(pop, 3, [1,2])

  def testNewPopByIndInfo(self):
    'Testing function newPopByIndInfo'
    pop = population(subPop=[1,2,3,4], ploidy=1, loci=[1])
    arr = pop.arrGenotype()
    arr[:] = range(pop.popSize())
    oldPop = pop.clone()
    #
    pop1 = pop.newPopByIndID(info=[-1,0,1,1,2,1,1,2,-1,0])
    self.assertEqual( pop, oldPop)
    self.assertEqual( pop1.subPopSizes(), (2,4,2))
    # subpop will be shifted
    self.assertGenotype(pop1, 0, [1,9])
    self.assertGenotype(pop1, 1, [2,3,5,6])
    self.assertGenotype(pop1, 2, [4,7])
    # change new pop will be change old one
    pop1.individual(0).setAllele(1, 0)
    self.assertNotEqual(pop.individual(0).allele(0), 1)

  def testRemoveLoci(self):
    'Testing function removeLoci'
    pop = population(subPop=[1,2], ploidy=2, loci=[2,3,1])
    arr = pop.arrGenotype()
    arr[:] = range(pop.totNumLoci())*(pop.popSize()*pop.ploidy())
    pop.removeLoci(remove=[2])
    self.assertEqual( pop.numChrom(), 3)
    self.assertEqual( pop.numLoci(0), 2)
    self.assertEqual( pop.numLoci(1), 2)
    self.assertEqual( pop.numLoci(2), 1)
    self.assertEqual( pop.arrGenotype().count(2), 0)
    pop.removeLoci(remove=[4])
    self.assertEqual( pop.numChrom(), 2)
    self.assertEqual( pop.numLoci(0), 2)
    self.assertEqual( pop.numLoci(1), 2)
    self.assertEqual( pop.arrGenotype().count(5), 0)
    # keep
    pop.removeLoci(keep=[1,2])
    self.assertEqual( pop.numChrom(), 2)
    self.assertEqual( pop.numLoci(0), 1)
    self.assertEqual( pop.numLoci(1), 1)
    if alleleType() == 'binary':
      self.assertEqual( pop.arrGenotype().count(3), 0 )
      self.assertEqual( pop.arrGenotype().count(1), pop.popSize()*pop.ploidy()*2 )
    else:
      self.assertEqual( pop.arrGenotype().count(3), pop.popSize()*pop.ploidy() )
      self.assertEqual( pop.arrGenotype().count(1), pop.popSize()*pop.ploidy() )
    
  def testPopInfo(self):
    'Testing population info related functions'
    pop = population(subPop=[1,2])
    pop.setIndSubPopIDWithID()
    self.assertEqual(pop.exposeIndInfo(), [0,1,1])
    #
    pop.setIndSubPopID([3,4,5])
    self.assertEqual(pop.exposeIndInfo(), [3,4,5])
    self.assertRaises(exceptions.ValueError, 
      pop.setIndSubPopID, [2,3])      
    #
    pop.individual(1).setAffected(True)
    self.assertEqual(pop.exposeAffectedness(), [0,1,0])
    #
    
  def testArrGenotype(self):
    'Testing function arrGenotype'
    pop = population(loci=[1,2], subPop=[1,2])
    arr = pop.arrGenotype()
    self.assertEqual( len(arr), pop.genoSize()*pop.popSize())
    arr = pop.arrGenotype(1)
    self.assertEqual( len(arr), pop.genoSize()*pop.subPopSize(1))
    arr[0] = 1
    self.assertEqual( pop.individual(0,1).allele(0), 1)
    self.assertRaises(exceptions.IndexError, 
      pop.arrGenotype, 2)
    # arr assignment
    arr[:] = 1
    self.assertEqual( pop.individual(0,0).arrGenotype(), [0]*pop.genoSize())
    self.assertEqual( pop.individual(0,1).arrGenotype(), [1]*pop.genoSize())
    self.assertEqual( pop.individual(1,1).arrGenotype(), [1]*pop.genoSize())
        
  def testCompare(self):
    'Testing population comparison'
    pop = population(10, loci=[2])
    pop1 = population(10, loci=[2])
    self.assertEqual( pop == pop1, True)
    pop.individual(0).setAllele(1, 0)
    self.assertEqual( pop == pop1, False)
    pop1 = population(10, loci=[3])
    pop1.individual(0).setAllele(1, 0)
    # false becase of geno structure difference.
    self.assertEqual( pop == pop1, False)

  def testSaveLoadPopulation(self):
    'Testing save and load populations'
    if alleleType() != 'binary':
      pop = population(size=10, ploidy=2, loci=[5, 7], 
        subPop=[2, 8], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        maxAllele=4, alleleNames=['_','A','C','T','G']) 
      InitByFreq(pop, [.2, .3, .5])
    else:
      pop = population(size=10, ploidy=2, loci=[5, 7], 
        subPop=[2, 8], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
        alleleNames=['1','2']) 
      InitByFreq(pop, [.2, .8])
    # .xml format may not be available (under mac)
    for file in ['a.txt', 'a.bin', 'a.xml', 'a.txt.gz', 'a.bin.gz', 'a.xml.gz']:
      if (file == 'a.xml' or file == 'a.xml.gz') and not supportXML():
        continue
      pop.savePopulation(file, compress=False)
      assert os.path.isfile(file)
      pop1 = LoadPopulation(file)
      self.assertEqual(pop, pop1)
      pop.savePopulation(file, compress=True)
      assert os.path.isfile(file)
      pop1 = LoadPopulation(file)
      os.remove(file)
    # can load file with wrong extension
    # can load file with wrong extension
    pop.savePopulation('a.txt', format='bin')
    pop1 = LoadPopulation('a.txt')
    self.assertEqual(pop, pop1)
    os.remove('a.txt')
    #
    # save load several populations
    # make pop and pop1 different
    pop.individual(0).setAllele(0,1)
    pop1.individual(0).setAllele(1,1)
    self.assertNotEqual(pop, pop1)
    SavePopulations([pop, pop1], 'a.txt')
    (pop2, pop3) = LoadPopulations('a.txt')
    self.assertEqual(pop, pop2)
    self.assertEqual(pop1, pop3)
    self.assertNotEqual(pop, pop3)
    os.remove('a.txt')    
    
  def testPopVars(self):
    'Testing population variables'
    pop = population()
    self.assertEqual( pop.grp(), -1)
    self.assertEqual( pop.rep(), -1)
    # var will be copied?
    pop.dvars().x = 1
    pop1 = pop.clone()
    self.assertEqual( pop1.dvars().x, 1)
    pop1.dvars().y = 2
    self.assertEqual(pop.vars().has_key('y'), False)
    # test if a variable will be fully updated (a previous bug)
    # If a vector has 16 numbers, and then it will have a value
    # of 10 numbers, will the rest of the 6 numbers be removed?
    pop = population(1000, loci=[2,4])
    InitByFreq(pop, [.2, .3, .5])
    Stat(pop, alleleFreq=range(0,6))
    self.assertEqual(len(pop.dvars().alleleFreq), 6)
    pop.removeLoci(remove=[0,4])
    Stat(pop, alleleFreq=range(0,4))
    self.assertEqual(len(pop.dvars().alleleFreq), 4)

  def testAncestry(self):
    'Testing ancestral population related functions'
    pop = population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
    InitByFreq(pop, [.2,.8])
    gt = list(pop.arrGenotype())
    self.assertEqual(pop.ancestralDepth(), 0)
    pop1 = population(subPop=[2,3], loci=[2,3], ancestralDepth=2)
    InitByFreq(pop1, [.8,.2])
    gt1 = list(pop1.arrGenotype())
    # can not do, because of different genotype
    pop.pushAndDiscard(pop1)
    # pop1 should be empty now
    self.assertEqual(pop1.popSize(), 0)
    # pop should have one ancestry
    self.assertEqual(pop.ancestralDepth(), 1)
    # genotype of pop should be pop1
    self.assertEqual(pop.arrGenotype(), gt1)
    # use ancestry pop
    pop.useAncestralPop(1)
    self.assertEqual(pop.arrGenotype(), gt)
    # use back
    pop.useAncestralPop(0)
    self.assertEqual(pop.arrGenotype(), gt1)
    # can not push itself
    self.assertRaises(exceptions.ValueError,
      pop.pushAndDiscard, pop)
    # can not do, because of different genotype
    pop2 = population(subPop=[3,5], loci=[2])
    self.assertRaises(exceptions.ValueError,
      pop.pushAndDiscard, pop2)
    # [ gt1, gt ]
    # push more?
    pop2 = population(subPop=[3,5], loci=[2,3])
    InitByFreq(pop2, [.2,.8])
    gt2 = list(pop2.arrGenotype())
    pop3 = pop2.clone()
    pop.pushAndDiscard(pop2)
    # [ gt2, gt1, gt]
    # pop should have one ancestry
    self.assertEqual(pop2.popSize(), 0)
    self.assertEqual(pop.ancestralDepth(), 2)
    # 
    self.assertEqual(pop.arrGenotype(), gt2)
    pop.useAncestralPop(1)
    self.assertEqual(pop.arrGenotype(), gt1)
    pop.useAncestralPop(2)
    self.assertEqual(pop.arrGenotype(), gt)
    pop.useAncestralPop(0)
    #
    pop.pushAndDiscard(pop3)
    # [ gt2, gt2, gt1]
    self.assertEqual(pop3.popSize(), 0)
    self.assertEqual(pop.ancestralDepth(), 2)
    pop.useAncestralPop(1)
    self.assertEqual(pop.arrGenotype(), gt2)
    pop.useAncestralPop(2)
    self.assertEqual(pop.arrGenotype(), gt1)

if __name__ == '__main__':
  unittest.main()
