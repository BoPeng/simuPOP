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
        gt = list(pop.arrGenotype(subPop, True))
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

    def testIterator(self):
        'Testing the individual iterators'
        pop = population(loci=[1], subPop=[4,6])
        for ind in pop.individuals():
            ind.setAllele(1, 0)
        for ind in pop.individuals(1):
            ind.setAllele(2, 0)
        for ind in pop.individuals(0):
            self.assertEqual(ind.allele(0), 1)
        if alleleType() == 'binary':
            for ind in pop.individuals(1):
                self.assertEqual(ind.allele(0), 1)
        else:
            for ind in pop.individuals(1):
                self.assertEqual(ind.allele(0), 2)
     
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
        #
        pop = population(size=10, infoFields=['age'])
        pop.individual(1).setInfo(1, 'age')
        self.assertEqual(pop.individual(1).info('age'), 1)
        self.assertEqual(pop.individual(1).info('age'), 1)
        self.assertEqual(pop.individual(1).info('age'), 1)
        pop.setIndInfo(range(10), 'age')
        self.assertEqual(pop.individual(0).info('age'), 0)
        #print pop.indInfo('age', True)
        pop.setSubPopStru(newSubPopSizes=[2,8], allowPopSizeChange=False)
        for i in range(10):
            self.assertEqual(pop.individual(i).info('age'), i)
        

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
        # test if info is swapped
        pop = population(10, infoFields=['age'])
        pop.setIndInfo(range(10), 'age')
        pop1 = population(5, infoFields=['fitness'])
        pop1.setIndInfo(range(10,15), 'fitness')
        pop.swap(pop1)
        self.assertEqual(pop.infoField(0), 'fitness')
        self.assertEqual(pop1.infoField(0), 'age')
        for i in range(5):
            self.assertEqual(pop.individual(i).info('fitness'), i+10)
        for i in range(10):
            self.assertEqual(pop1.individual(i).info('age'), i)
        ###
        ###  ARRINDINFO is only valid for head node.
        ###
        if mpiRank() == 0:
            self.assertEqual(pop.arrIndInfo(True), range(10,15))
            self.assertEqual(pop1.arrIndInfo(True), range(10))
        else:
            self.assertRaises(exceptions.ValueError, pop.arrIndInfo, True)
            self.assertRaises(exceptions.ValueError, pop1.arrIndInfo, True)
        self.assertEqual(pop.popSize(), 5)
        self.assertEqual(pop1.popSize(), 10)

    def testSplitSubPop(self):
        'Testing function splitSubPop'
        pop = population(subPop=[5,6,7], ploidy=1, loci=[1])
        arr = pop.arrGenotype(True)
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
        pop.arrGenotype(True)[:] = range(pop.popSize())
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
        arr = pop.arrGenotype(True)
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
        pop.arrGenotype(True)[:] = range(pop.popSize())
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
        arr = pop.arrGenotype(True)
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
        # see if remove subPOp change info
        pop = population(subPop=[2,4,5], infoFields=['age'])
        self.assertEqual( pop.numSubPop(), 3)
        pop.setIndInfo(range(11), 0)
        if mpiRank() == 0:
            self.assertEqual(pop.arrIndInfo(True), range(11))
        pop.removeSubPops([1])
        if mpiRank() == 0:        
            self.assertEqual(pop.arrIndInfo(0, True), range(2))
            self.assertEqual(pop.arrIndInfo(1, True), range(6,11))
    
    def testRemoveIndividuals(self):
        'Testing function removeIndividuals'
        pop = population(subPop=[0,1,0,2,3,0], ploidy=1, loci=[1])
        arr = pop.arrGenotype(True)
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
        # see if remove individual change info
        pop = population(subPop=[2,4,5], infoFields=['age'])
        pop.setIndInfo(range(11), 0)
        pop.removeIndividuals([2,3,4,5])
        self.assertEqual(pop.subPopSizes(), (2,0,5))
        if mpiRank() == 0:
            self.assertEqual(pop.arrIndInfo(0, True), range(2))
            self.assertEqual(pop.arrIndInfo(2, True), range(6,11))
            # nothing in the middle left
            self.assertEqual(pop.arrIndInfo(True), range(2) + range(6,11))
    
    
    def testMergeSubPops(self):
        'Testing function mergeSubPops'
        pop = population(subPop=[0,1,0,2,3,0], ploidy=1, loci=[1])
        arr = pop.arrGenotype(True)
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
        # see if merging affect individual id.
        pop = population(subPop=[2,4,5], infoFields=['age'])
        pop.setIndInfo(range(11), 0)
        pop.mergeSubPops([0,2], removeEmptySubPops=True)
        self.assertEqual(pop.subPopSizes(), (7,4))
        # the order may be different
        if mpiRank() == 0:
            self.assertEqual(sum(pop.arrIndInfo(0, False)), sum(range(2)+range(6,11)))
            self.assertEqual(pop.arrIndInfo(1, True), range(2,6))
        

    def testReorderSubPops(self):
        'Testing function reorderSubPops'
        pop = population(subPop=[1,2,3,4], ploidy=1, loci=[1])
        arr = pop.arrGenotype(True)
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
        arr = pop.arrGenotype(True)
        arr[:] = range(pop.popSize())
        #
        pop.reorderSubPops(rank=[1,3,0,2])
        self.assertEqual( pop.subPopSizes(), (3,1,4,2))
        # subpop will be shifted
        self.assertGenotype(pop, 0, [3,4,5])
        self.assertGenotype(pop, 1, [0])
        self.assertGenotype(pop, 2, [6,7,8,9])
        self.assertGenotype(pop, 3, [1,2])
        # reorder does not change info
        pop = population(subPop=[1,2,3,4], ploidy=1, loci=[1], infoFields=['age'])
        pop.setIndInfo(range(10), 'age')
        #
        pop.reorderSubPops(rank=[1,3,0,2])
        self.assertEqual( pop.subPopSizes(), (3,1,4,2))
        newInfoSums = [sum([3,4,5]), sum([0]), sum([6,7,8,9]), sum([1,2])]
        if mpiRank() == 0:
            for i in range(4):
                self.assertEqual(sum(pop.arrIndInfo(i, True)), newInfoSums[i])
        

    def testNewPopByIndID(self):
        'Testing function newPopByIndInfo'
        pop = population(subPop=[1,2,3,4], ploidy=1, loci=[1])
        arr = pop.arrGenotype(True)
        arr[:] = range(pop.popSize())
        oldPop = pop.clone()
        #
        pop1 = pop.newPopByIndID(id=[-1,0,1,1,2,1,1,2,-1,0])
        self.assertEqual( pop, oldPop)
        self.assertEqual( pop1.subPopSizes(), (2,4,2))
        # subpop will be shifted
        self.assertGenotype(pop1, 0, [1,9])
        self.assertGenotype(pop1, 1, [2,3,5,6])
        self.assertGenotype(pop1, 2, [4,7])
        # change new pop will be change old one
        pop1.individual(0).setAllele(1, 0)
        self.assertNotEqual(pop.individual(0).allele(0), 1)
        # does new pop keeps info
        pop = population(subPop=[1,2,3,4], ploidy=1, loci=[1], infoFields=['age'])
        pop.setIndInfo(range(10),'age')
        pop1 = pop.newPopByIndID(id=[-1,8,7,6,5,4,3,2,1,-1], removeEmptySubPops=True)
        self.assertEqual(pop1.popSize(), 8)
        self.assertEqual(pop1.subPopSizes(), tuple([1]*8))
        for i in range(8):
            self.assertEqual(pop1.individual(i).info('age'), 8-i)
        

    def testRemoveLoci(self):
        'Testing function removeLoci'
        # FIXME: make sure to test a case when 
        # genotype has to be moved from one node to another 
        # in the MPI modules
        pop = population(subPop=[1,2], ploidy=2, loci=[2,3,1])
        arr = pop.arrGenotype(True)
        arr[:] = range(pop.totNumLoci())*(pop.popSize()*pop.ploidy())
        pop.removeLoci(remove=[2])
        self.assertEqual( pop.numChrom(), 3)
        self.assertEqual( pop.numLoci(0), 2)
        self.assertEqual( pop.numLoci(1), 2)
        self.assertEqual( pop.numLoci(2), 1)
        self.assertEqual( pop.arrGenotype(True).count(2), 0)
        pop.removeLoci(remove=[4])
        self.assertEqual( pop.numChrom(), 2)
        self.assertEqual( pop.numLoci(0), 2)
        self.assertEqual( pop.numLoci(1), 2)
        self.assertEqual( pop.arrGenotype(True).count(5), 0)
        # keep
        pop.removeLoci(keep=[1,2])
        self.assertEqual( pop.numChrom(), 2)
        self.assertEqual( pop.numLoci(0), 1)
        self.assertEqual( pop.numLoci(1), 1)
        if alleleType() == 'binary':
            self.assertEqual( pop.arrGenotype(True).count(3), 0 )
            self.assertEqual( pop.arrGenotype(True).count(1), pop.popSize()*pop.ploidy()*2 )
        else:
            self.assertEqual( pop.arrGenotype(True).count(3), pop.popSize()*pop.ploidy() )
            self.assertEqual( pop.arrGenotype(True).count(1), pop.popSize()*pop.ploidy() )
        
    def testArrGenotype(self):
        'Testing function arrGenotype'
        pop = population(loci=[1,2], subPop=[1,2])
        arr = pop.arrGenotype(True)
        self.assertEqual( len(arr), pop.genoSize()*pop.popSize())
        arr = pop.arrGenotype(1, True)
        self.assertEqual( len(arr), pop.genoSize()*pop.subPopSize(1))
        arr[0] = 1
        self.assertEqual( pop.individual(0,1).allele(0), 1)
        self.assertRaises(exceptions.IndexError, 
            pop.arrGenotype, 2, True)
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

        
    def testPopInfo(self):
        'Testing info-related functions'
        # create a population without info field
        pop = population(10)
        ind = pop.individual(0)
        self.assertRaises(exceptions.IndexError, ind.info, 0)
        # give name
        self.assertRaises(exceptions.ValueError, population, infoFields='age')
        pop = population(10, infoFields=['age'])
        self.assertEqual(pop.infoField(0), 'age')
        self.assertEqual(pop.infoFields(), ('age',))
        ind = pop.individual(0)
        self.assertRaises(exceptions.IndexError, ind.info, 1)
        # can set more names
        pop = population(10, infoFields=['age', 'fitness', 'trait1'])
        self.assertEqual(pop.infoField(0), 'age')
        self.assertEqual(pop.infoField(2), 'trait1')
        self.assertRaises(exceptions.IndexError, pop.infoField, 3)
        self.assertEqual(pop.infoFields(), ('age', 'fitness', 'trait1'))
        # can set and read each info
        pop.setIndInfo(range(10), 0)
        pop.setIndInfo(range(10,20), 1)
        pop.setIndInfo(range(20,30), 2)
        self.assertRaises(exceptions.IndexError, pop.setIndInfo, range(30,40), 3)
        for i in range(3):
            for j in range(10):
                self.assertEqual(pop.individual(j).info(i), i*10+j)
        ind = pop.individual(0)
        if mpiRank() == 0:
            self.assertEqual(ind.arrInfo(), (0,10,20))
            self.assertEqual(pop.arrIndInfo(True)[:8], (0,10,20,1,11,21,2,12))
            self.assertEqual(pop.arrIndInfo(0, True)[:8], (0,10,20,1,11,21,2,12))
        else:
            self.assertRaises(exceptions.ValueError, ind.arrInfo)
            self.assertRaises(exceptions.ValueError, pop.arrIndInfo, True)
            self.assertRaises(exceptions.ValueError, pop.arrIndInfo, 0, True)
        self.assertRaises(exceptions.IndexError, pop.arrIndInfo, 1, True)
        # access by name
        #pop.setIndInfo(range(30,40), 'sex')
        #pop.arrIndInfo('sex')
        self.assertRaises(exceptions.IndexError, ind.info, 'sex')
        ind.setInfo(18, 'age')
        self.assertEqual(ind.info('age'), 18)
        #
        pop = population(10, infoFields=['age', 'fitness'])
        self.assertEqual(pop.hasInfoField('age'), True)
        self.assertEqual(pop.hasInfoField('fitness'), True)
        self.assertEqual(pop.hasInfoField('misc'), False)
        self.assertEqual(pop.infoFields(), ('age', 'fitness'))
        self.assertEqual(pop.infoSize(), 2)
        ind = pop.individual(0)
        # set info
        ind.setInfo(2, 0)
        self.assertEqual(ind.info('age'), 2)
        # get info
        self.assertEqual(ind.info(0), 2)
        # create another field
        self.assertEqual(pop.infoIdx('age'), 0)
        # set values
        pop.setIndInfo([1, 2,3,4,5,6,7,8,9,10], 0)
        for i in range(10):
            self.assertEqual(pop.individual(i).info(0), i+1)
        self.assertEqual(pop.infoIdx('fitness'), 1)
        # adding fitness should not interfere with age.
        for i in range(10):
            self.assertEqual(pop.individual(i).info(0), i+1)
        # set info by name
        pop = population(10, infoFields=['age', 'fitness'])
        #print pop.indInfo('fitness', True)
        #print pop.indInfo('age', True)
        for i in range(10):
            pop.individual(i).setInfo(i+50, 'fitness')
            self.assertEqual(pop.individual(i).info('fitness'), i+50)
        pop.setIndInfo(range(50,60), 'fitness')
        for i in range(10):
            pop.individual(i).setInfo(i+50, 'fitness')
            self.assertEqual(pop.individual(i).info('fitness'), i+50)
        # 
        #  test indInfo
        pop = population(subPop=[4,6], infoFields=['age', 'fitness'])
        pop.setIndInfo(range(10), 'age')
        pop.setIndInfo(range(100, 110), 'fitness')
        self.assertEqual(pop.indInfo('age', True), tuple([float(x) for x in range(10)]))
        self.assertEqual(pop.indInfo('fitness', True), tuple([float(x) for x in range(100, 110)]))
        self.assertEqual(pop.indInfo('age', 1, True), tuple([float(x) for x in range(4, 10)]))
        self.assertEqual(pop.indInfo('fitness', 0, True), tuple([float(x) for x in range(100, 104)]))
        #
        # test reset info fields
        pop = population(size=10, infoFields=['age'])
        pop.setInfoFields(['age', 'fitness'])
        self.assertEqual(pop.infoSize(), 2)
        self.assertEqual(pop.infoFields(), ('age', 'fitness'))
        # set values
        pop.setIndInfo(range(10), 'age')
        pop.setIndInfo(range(100, 110), 'fitness')
        self.assertEqual(pop.indInfo('age', True), tuple([float(x) for x in range(10)]))
        self.assertEqual(pop.indInfo('fitness', True), tuple([float(x) for x in range(100, 110)]))
        # add an existing field
        pop.addInfoField('fitness')
        self.assertEqual(pop.infoSize(), 2)
        # add a new one.
        pop.addInfoField('misc')
        self.assertEqual(pop.infoSize(), 3)
        pop.setIndInfo(range(200, 210), 'fitness')
        self.assertEqual(pop.indInfo('age', True), tuple([float(x) for x in range(10)]))
        self.assertEqual(pop.indInfo('fitness', True), tuple([float(x) for x in range(200, 210)]))

        
    def testPopVars(self):
        'Testing population variables'
        # FIXME: currently variable is stored on all nodes.
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
        gt = list(pop.arrGenotype(True))
        self.assertEqual(pop.ancestralDepth(), 0)
        pop1 = population(subPop=[2,3], loci=[2,3], ancestralDepth=2)
        InitByFreq(pop1, [.8,.2])
        gt1 = list(pop1.arrGenotype(True))
        # can not do, because of different genotype
        pop.pushAndDiscard(pop1)
        # pop1 should be empty now
        self.assertEqual(pop1.popSize(), 0)
        # pop should have one ancestry
        self.assertEqual(pop.ancestralDepth(), 1)
        # genotype of pop should be pop1
        self.assertEqual(pop.arrGenotype(True), gt1)
        # use ancestry pop
        pop.useAncestralPop(1)
        self.assertEqual(pop.arrGenotype(True), gt)
        # use back
        pop.useAncestralPop(0)
        self.assertEqual(pop.arrGenotype(True), gt1)
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
        gt2 = list(pop2.arrGenotype(True))
        pop3 = pop2.clone()
        pop.pushAndDiscard(pop2)
        # [ gt2, gt1, gt]
        # pop should have one ancestry
        self.assertEqual(pop2.popSize(), 0)
        self.assertEqual(pop.ancestralDepth(), 2)
        # 
        self.assertEqual(pop.arrGenotype(True), gt2)
        pop.useAncestralPop(1)
        self.assertEqual(pop.arrGenotype(True), gt1)
        pop.useAncestralPop(2)
        self.assertEqual(pop.arrGenotype(True), gt)
        pop.useAncestralPop(0)
        #
        pop.pushAndDiscard(pop3)
        # [ gt2, gt2, gt1]
        self.assertEqual(pop3.popSize(), 0)
        self.assertEqual(pop.ancestralDepth(), 2)
        pop.useAncestralPop(1)
        self.assertEqual(pop.arrGenotype(True), gt2)
        pop.useAncestralPop(2)
        self.assertEqual(pop.arrGenotype(True), gt1)
        #
        # test if ind info is saved with ancestral populations
        pop = population(10, ancestralDepth=2, infoFields=['age', 'fitness'])
        pop.setIndInfo(range(10), 'age')
        pop.setIndInfo(range(10, 20), 'fitness')
        pop1 = population(20, infoFields=['age', 'fitness'])
        pop1.setIndInfo(range(100, 120), 'age')
        pop1.setIndInfo(range(110, 130), 'fitness')
        pop.pushAndDiscard(pop1)
        # test info
        self.assertEqual(pop.popSize(), 20)
        if mpiRank() == 0:
            self.assertEqual(pop.arrIndInfo(True)[:4], [100, 110, 101, 111])
        pop.useAncestralPop(1)
        self.assertEqual(pop.popSize(), 10)
        if mpiRank() == 0:
            self.assertEqual(pop.arrIndInfo(True)[:4], [0, 10, 1, 11])


    def testSaveLoadPopulation(self):
        'Testing save and load populations'
        if alleleType() != 'binary':
            pop = population(size=10, ploidy=2, loci=[5, 7], 
                subPop=[2, 8], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
                maxAllele=4, alleleNames=['_','A','C','T','G'],
                infoFields=['age', 'fitness']) 
            InitByFreq(pop, [.2, .3, .5])
            pop.setIndInfo(range(10), 'age')
            pop.setIndInfo(range(100, 110), 'fitness')
        else:
            pop = population(size=10, ploidy=2, loci=[5, 7], 
                subPop=[2, 8], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
                alleleNames=['1','2'], infoFields=['age', 'fitness']) 
            InitByFreq(pop, [.2, .8])
            pop.setIndInfo(range(10), 'age')
            pop.setIndInfo(range(100, 110), 'fitness')
        # .xml format may not be available (under mac)
        for file in ['a.txt', 'a.bin', 'a.xml', 'a.txt.gz', 'a.bin.gz', 'a.xml.gz']:
            if (file == 'a.xml' or file == 'a.xml.gz') and not supportXML():
                continue
            pop.savePopulation(file, compress=False)
            assert os.path.isfile(file), "File %s does not exist" % file
            pop1 = LoadPopulation(file)
            self.assertEqual(pop, pop1)
            pop.savePopulation(file, compress=True)
            assert os.path.isfile(file), "File %s does not exist" % file
            pop1 = LoadPopulation(file)
            self.assertEqual(pop, pop1)
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

    def testCrossSaveLoad(self):
        'Testing population saved by other modules'
        if alleleType() == 'binary':
            if os.path.isfile('test_ba.txt'):
                pop = LoadPopulation('test_ba.txt')
            else:
                pop = population(size=10, ploidy=2, loci=[5, 7], 
                    subPop=[2, 8], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
                    alleleNames=['1', '2'], infoFields=['age', 'fitness'], ancestralDepth=2) 
                InitByFreq(pop, [.2, .8])
                pop.setIndInfo(range(10), 'age')
                pop.setIndInfo(range(100, 110), 'fitness')
                simu = simulator(pop, randomMating(), rep=1)
                simu.evolve(ops=[], end=2)
                simu.population(0).savePopulation('test_ba.txt')
            # try to load file
            if os.path.isfile('test_std.txt'):
                pop1 = LoadPopulation('test_std.txt')
                self.assertEqual(pop, pop1)
            if os.path.isfile('test_la.txt'):
                pop1 = LoadPopulation('test_la.txt')
                self.assertEqual(pop, pop1)            
        elif alleleType() == 'short':
            if os.path.isfile('test_ba.txt'):
                pop = LoadPopulation('test_ba.txt')
                pop.savePopulation('test_std.txt')
            else:
                return
            if os.path.isfile('test_la.txt'):
                pop1 = LoadPopulation('test_la.txt')
                self.assertEqual(pop, pop1)
        else:
            if os.path.isfile('test_ba.txt'):
                pop = LoadPopulation('test_ba.txt')
                pop.savePopulation('test_la.txt')
            else:
                return
            if os.path.isfile('test_std.txt'):
                pop1 = LoadPopulation('test_std.txt')
                self.assertEqual(pop, pop1)


    def testMergePopulation(self):
        'Testing merge populations...'
        pop = population(subPop=[7, 3, 4], loci=[4, 5, 1])
        pop1 = population(subPop=[4, 5], loci=[4,5,1])
        pop2 = population(subPop=[4, 5], loci=[4,1])
        InitByFreq(pop, [.2, .3, .5])
        InitByFreq(pop1, [.5, .5])
        InitByFreq(pop2, [.2, .8])
        # merge error (number of subpop mismatch)
        self.assertRaises(exceptions.ValueError, pop.mergePopulation, pop2)
        # merge without subpop size change
        pop_ori = pop.clone()
        pop.mergePopulation(pop1)
        self.assertEqual(pop.subPopSizes(), (7,3,4,4,5))
        for sp in range(3):
            for i in range(pop.subPopSize(sp)):
                self.assertEqual(pop.individual(i, sp), pop_ori.individual(i, sp))
        for sp in range(2):
            for i in range(pop.subPopSize(3+sp)):
                self.assertEqual(pop.individual(i, sp+3), pop1.individual(i, sp))
        # merge with new subpop sizes
        pop = pop_ori.clone()
        # total size should not change
        self.assertRaises(exceptions.ValueError, pop.mergePopulation, pop1, newSubPopSizes=[5, 20])
        #
        pop.mergePopulation(pop1, newSubPopSizes=[9, 10,4])
        self.assertEqual(pop.subPopSizes(), (9, 10, 4))
        for i in range(pop_ori.popSize()):
            self.assertEqual(pop.individual(i), pop_ori.individual(i))
        for i in range(pop1.popSize()):
            self.assertEqual(pop.individual(i+pop_ori.popSize()), pop1.individual(i))
        #
        # test the Merge function
        pop = population(subPop=[7, 3, 4], loci=[4, 5, 1])
        pop1 = population(subPop=[4, 5], loci=[4,5,1])
        pop2 = population(subPop=[4, 5], loci=[4,1])
        InitByFreq(pop, [.2, .3, .5])
        pop_ori = pop.clone()
        InitByFreq(pop1, [.5, .5])
        pop1_ori = pop1.clone()
        InitByFreq(pop2, [.2, .8])
        pop2_ori = pop2.clone()
        #self.assertRaises(exceptions.ValueError, MergePopulations, [pop, pop2])
        #
        mp = MergePopulations(pops=[pop, pop1, pop1])
        # populaition not changed
        self.assertEqual(pop, pop_ori)
        self.assertEqual(pop1, pop1_ori)
        pop.mergePopulation(pop1)
        pop.mergePopulation(pop1)
        self.assertEqual(pop, mp)
        #
        # test for merge of ancestral gen
        pop = population(subPop=[7, 3, 4], loci=[4, 5, 1])
        pop.setAncestralDepth(-1)
        pop1 = population(subPop=[4, 5], loci=[4, 5, 1])
        pop.pushAndDiscard(pop1)
        #
        pop1 = pop.clone()
        #
        pop2 = MergePopulations([pop, pop1])
        self.assertEqual(pop2.ancestralDepth(), 1)
        self.assertEqual(pop2.subPopSizes(), (4,5,4,5))
        pop2.useAncestralPop(1)
        self.assertEqual(pop2.subPopSizes(), (7,3,4,7,3,4))
        # test for keepAncestralPops (FIXME)
        
    def testMergePopulationByLoci(self):
        'Testing merge populations...'
        pop = population(subPop=[7, 3, 4], loci=[4, 5, 1])
        pop1 = population(subPop=[4, 5], loci=[4,5,1])
        pop2 = population(subPop=[4, 5], loci=[4,1])
        InitByFreq(pop, [.2, .3, .5])
        InitByFreq(pop1, [.5, .5])
        InitByFreq(pop2, [.2, .8])
        # merge error (number of subpop mismatch)
        self.assertRaises(exceptions.ValueError, pop.mergePopulationByLoci, pop2)
        # merge without subpop size change
        pop_ori = pop.clone()
        pop.mergePopulationByLoci(pop1)
        self.assertEqual(pop.subPopSizes(), (7,3,4,4,5))
        for sp in range(3):
            for i in range(pop.subPopSize(sp)):
                self.assertEqual(pop.individual(i, sp), pop_ori.individual(i, sp))
        for sp in range(2):
            for i in range(pop.subPopSize(3+sp)):
                self.assertEqual(pop.individual(i, sp+3), pop1.individual(i, sp))
        # merge with new subpop sizes
        pop = pop_ori.clone()
        # total size should not change
        self.assertRaises(exceptions.ValueError, pop.mergePopulationByLoci, pop1, newSubPopSizes=[5, 20])
        #
        pop.mergePopulationByLoci(pop1, newSubPopSizes=[9, 10,4])
        self.assertEqual(pop.subPopSizes(), (9, 10, 4))
        for i in range(pop_ori.popSize()):
            self.assertEqual(pop.individual(i), pop_ori.individual(i))
        for i in range(pop1.popSize()):
            self.assertEqual(pop.individual(i+pop_ori.popSize()), pop1.individual(i))
        #
        # test the Merge function
        pop = population(subPop=[7, 3, 4], loci=[4, 5, 1])
        pop1 = population(subPop=[4, 5], loci=[4,5,1])
        pop2 = population(subPop=[4, 5], loci=[4,1])
        InitByFreq(pop, [.2, .3, .5])
        pop_ori = pop.clone()
        InitByFreq(pop1, [.5, .5])
        pop1_ori = pop1.clone()
        InitByFreq(pop2, [.2, .8])
        pop2_ori = pop2.clone()
        #self.assertRaises(exceptions.ValueError, MergePopulationByLocis, [pop, pop2])
        #
        mp = MergePopulationByLocis(pops=[pop, pop1, pop1])
        # populaition not changed
        self.assertEqual(pop, pop_ori)
        self.assertEqual(pop1, pop1_ori)
        pop.mergePopulationByLoci(pop1)
        pop.mergePopulationByLoci(pop1)
        self.assertEqual(pop, mp)
        #
        # test for merge of ancestral gen
        pop = population(subPop=[7, 3, 4], loci=[4, 5, 1])
        pop.setAncestralDepth(-1)
        pop1 = population(subPop=[4, 5], loci=[4, 5, 1])
        pop.pushAndDiscard(pop1)
        #
        pop1 = pop.clone()
        #
        pop2 = MergePopulationByLocis([pop, pop1])
        self.assertEqual(pop2.ancestralDepth(), 1)
        self.assertEqual(pop2.subPopSizes(), (4,5,4,5))
        pop2.useAncestralPop(1)
        self.assertEqual(pop2.subPopSizes(), (7,3,4,7,3,4))


    def testResizePopulation(self):
        'Testing population resize...'
        pop = population(subPop=[7,3,4], loci=[4,5,1])
        InitByFreq(pop, [.2, .3, .5])
        pop1 = pop.clone()
        pop2 = pop.clone()
        # resize error (number of subpop mismatch)
        self.assertRaises(exceptions.ValueError, pop1.resize, [5,5])
        # resize without propagation
        pop1.resize([5, 5, 8], propagate=False)
        for sp in range(pop1.numSubPop()):
            for i in range(min(pop1.subPopSize(sp), pop.subPopSize(sp))):
                self.assertEqual(pop1.individual(i, sp), pop.individual(i, sp))
            for i in range(min(pop1.subPopSize(sp), pop.subPopSize(sp)), pop1.subPopSize(sp)):
                self.assertEqual(pop1.individual(i, sp).arrGenotype(), [0]*20)
        # resize with propagation
        pop2.resize([5, 5, 8], propagate=True)
        for sp in range(pop1.numSubPop()):
            for i in range(pop2.subPopSize(sp)):
                self.assertEqual(pop2.individual(i, sp), pop.individual(i%pop.subPopSize(sp), sp))


if __name__ == '__main__':
    unittest.main()
