#!/usr/bin/env python
#
# Purpose:
#!/usr/bin/env python
#
# This is a unittest file for population object
#
# Bo Peng (bpeng@rice.edu)
# 
# Last Modified: Sep, 2005
# 
# 

from simuPOP import *
import unittest, os, sys

class TestPopulation(unittest.TestCase):
   
  def setUp(self):
    ' set up a population with given property'
    self.pop = population(size=100, ploidy=2, loci=[5, 7], \
      subPop=[20, 80], lociDist=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], \
      maxAllele=4, alleleNames=['_','A','C','T','G']) 

  def testNewPopulation1(self):
    ' no exception '
    population(size=100, ploidy=2, loci=[5, 7], subPop=[20, 80],
      lociDist=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
      maxAllele=4, alleleNames=['_','A','C','T','G'])
    
  def testNewPopulation2(self):
    self.assertRaises(ValueError, 
      population(size=1000, ploidy=4, loci=[5, 7]*20, \
      subPop=[20, 80]*20))

  def testNewPopulation3(self):
    self.assertRaises(ValueError, 
      population(subPop=[20,20], ploidy=0) )

  def testNewPopulation3(self):
    self.assertRaises(TypeError, 
      population(subPop=[2], lociDist=[1,2]) )
      
  def testNewPopulation4(self):
    population(subPop=[2], loci=[2], lociDist=[[1,2]])

  def testNewPopulation5(self):
    population(subPop=[2], loci=[2], lociDist=[[1,2]])

  def testNewPopulation6(self):
    self.assertRaises(TypeError,   
      population(subPop=[2], loci=2, lociDist=[[2,1]]) )
    
  def testNewPopulation6(self):
    self.assertRaises(ValueError, 
      population(subPop=[2], loci=[2], lociDist=[[2,1]]) )

  def testPopSize(self):
    self.assertEqual(self.pop.popSize(), 100)

  def testPloidy(self):
    self.assertEqual( self.pop.ploidy(), 2)

  def testPloidyName(self):
    self.assertEqual( self.pop.ploidyName(), 'diploid')
    
  def testNumChrom(self):
    self.assertEqual( self.pop.numChrom(), 2)
    
  def testNumLoci(self):
    self.assertEqual( self.pop.numLoci(0), 5)
    self.assertEqual( self.pop.numLoci(1), 7)
    self.assertRaises(IndexError,
      self.pop.numLoci(2) )
    
  def testLocusDist(self):
    self.assertEqual( self.pop.locusDist(10), 12)
    
  def testArrLociDist(self):
    self.assertEqual( len( self.pop.arrLociDist()), 12)
     
  def testChromBegin(self):
    self.assertEqual( self.pop.chromBegin(1), 5)
     
  def testChromEnd(self):
    self.assertEqual( self.pop.chromEnd(0), 5)
    
  def testLocusIndex(self):
    self.assertEqual( self.pop.absLocusIndex(1,5), 10)
    
  def testChromLocusPair(self):
    self.assertEqual( self.pop.chromLocusPair(10), (1,5) )
    
  def testLocusDist(self):
    self.assertEqual( self.pop.locusDist( a.absLocusIndex(1,2) ), 6)
    
  def testNumSubPop(self):
    self.assertEqual( self.pop.numSubPop(), 2)
    
  def testTotNumLoci(self):
    self.assertEqual( self.pop.totNumLoci(), 12)
    
  def testGenoSize(self):
    self.assertEqual( self.pop.genoSize(), 
      self.pop.totNumLoci()*self.pop.ploidy() )
    
  def testAbsIndIndex(self):
    self.assertEqual( self.pop.absIndIndex(1,1), 21)
    
  def testSubPopIndPair(self):
    self.assertEqual( self.pop.subPopIndPair(4), (0,4))
    
  def testSubPopSize(self):
    self.assertEqual( self.pop.subPopSize(1), 80)
    self.assertRaises( IndexError, self.pop.subPopSize, 2 )
    self.assertEqual( self.pop.subPopSizes(), (20,80) )
    
  def testSubPopBegin(self):
    self.assertEqual( self.pop.subPopBegin(1), 20)
    
  def testSubPopEnd(self):
    self.assertEqual( self.pop.subPopEnd(0), 20)
    
  def testSetIndInfoWithSubPopID(self):
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
      assert self.pop == pop1
      os.remove(file)
    
  def testArrGenotype(self):
    self.assertEqual( len( self.pop.arrGenotype() ),
      self.pop.genoSize() * self.pop.popSize() ) 
    self.assertEqual( len( self.pop.arrGenotype(0) ),
      self.pop.genoSize() * self.pop.popSize(0) ) 
    
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
