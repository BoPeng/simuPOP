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
import unittest, os, sys, exceptions, random, copy

class TestPopulation(unittest.TestCase):
    # define a few functions to create basic populations
    def getPop(self, scramble=False, VSP=False, size=[20, 80], loci = [1, 2],
            ancGen=0, *arg, **kwargs):
        pop = population(size=size, ploidy=2, loci=loci, infoFields=['x'],
            ancGen=ancGen, *arg, **kwargs)
        pop.setGenotype([random.randint(1, 5) for x in range(pop.popSize()*pop.ploidy())])
        pop.setIndInfo([random.random() for x in range(pop.popSize())], 'x')
        for i in range(ancGen):
            pop.push(self.getPop(size=size, loci=loci, ancGen=0, *arg, **kwargs))
        InitSex(pop)
        if VSP:
            pop.setVirtualSplitter(sexSplitter())
        if scramble:
            pop.scramble()
        return pop

    def testAbsIndIndex(self):
        'Testing population::absIndIndex(idx, subPop), popSize()'
        pop = self.getPop()
        # ind, subPop
        self.assertEqual(pop.absIndIndex(1, 1), 21)
        self.assertEqual(pop.absIndIndex(10, 0), 10)
        self.assertRaises(exceptions.IndexError, pop.absIndIndex, 0, 2 )
        self.assertEqual(pop.popSize(), 100)

    def testSubPop(self):
        'Testing population::subPopBegin(subPop), subPopEnd(subPop), numSubPop()'
        'subPopSize(subPop), subPopSizes(), subPopIndPair(idx)'
        pop = self.getPop()
        self.assertEqual(pop.subPopBegin(1), 20)
        self.assertRaises(exceptions.IndexError, pop.subPopBegin, 2 )
        self.assertEqual(pop.subPopEnd(0), 20)
        self.assertRaises(exceptions.IndexError, pop.subPopEnd, 2 )
        self.assertEqual(pop.numSubPop(), 2)
        self.assertEqual(pop.subPopSize(0), 20)
        self.assertEqual(pop.subPopSize(1), 80)
        self.assertRaises(exceptions.IndexError, pop.subPopSize, 2 )
        self.assertEqual(pop.subPopSizes(), (20, 80) )
        self.assertEqual(pop.subPopIndPair(21), (1, 1) )
        self.assertRaises(exceptions.IndexError, pop.subPopIndPair, 200 )
        Stat(pop, numOfMale=True)
        pop.setVirtualSplitter(sexSplitter())
        self.assertEqual(pop.subPopSize([1, 0]), pop.dvars(1).numOfMale)
        self.assertEqual(pop.subPopSize([1, 1]), pop.dvars(1).numOfFemale)

    def testVirtualSubPop(self):
        'Testing population::numVirtualSubPop(), setVirtualSplitter(splitter), subPopName(subPop)'
        pop = population(1000, infoFields=['x'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(infoSplitter('x', values=range(10, 15)))
        self.assertEqual(pop.numVirtualSubPop(), 5)
        self.assertEqual(pop.subPopName(0), "unnamed")
        self.assertEqual(pop.subPopName([0, 0]), "unnamed - x = 10")
        self.assertEqual(pop.subPopName([0, 1]), "unnamed - x = 11")
        self.assertEqual(pop.subPopName([0, 4]), "unnamed - x = 14")
        self.assertRaises(exceptions.IndexError, pop.subPopName, 1)
        self.assertRaises(exceptions.IndexError, pop.subPopName, [0, 5])
        # with given names
        pop = population(size=[200, 500], infoFields=['x'], subPopNames=['A', 'B'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(infoSplitter('x', values=range(10, 15)))
        self.assertEqual(pop.numVirtualSubPop(), 5)
        self.assertEqual(pop.subPopName(0), "A")
        self.assertEqual(pop.subPopName(1), "B")
        self.assertRaises(exceptions.IndexError, pop.subPopName, 2)
        self.assertEqual(pop.subPopName([0, 0]), "A - x = 10")
        self.assertEqual(pop.subPopName([0, 1]), "A - x = 11")
        self.assertEqual(pop.subPopName([1, 4]), "B - x = 14")
        self.assertRaises(exceptions.IndexError, pop.subPopName, [0, 5])

    def testIndividuals(self):
        'Testing function population::individuals(), individuals(subPop), individual(idx, subPop=0)'
        def testAllInd(pop):
            self.assertEqual(len(list(pop.individuals())), pop.popSize())
            self.assertEqual(len(list(pop.individuals(0))), pop.subPopSize(0))
            self.assertEqual(len(list(pop.individuals(1))), pop.subPopSize(1))
        testAllInd(self.getPop())
        testAllInd(self.getPop(True))
        testAllInd(self.getPop(False, True))
        testAllInd(self.getPop(True, True))
        pop = population([20, 80], loci = [5, 7], infoFields=['x'])
        pop.individual(0).setAllele(1, 0)
        self.assertEqual(pop.individual(0).allele(0), 1)

    def testGenotype(self):
        'Testing population::genotype(), genotype(subPop)'
        pop = population(loci=[1, 2], size=[1, 2])
        arr = pop.genotype()
        self.assertEqual(len(arr), pop.genoSize()*pop.popSize())
        arr = pop.genotype(1)
        self.assertEqual(len(arr), pop.genoSize()*pop.subPopSize(1))
        self.assertRaises(exceptions.IndexError, pop.genotype, 2)

    def testSetGenotype(self):
        'Testing population::setGenotype(geno), setGenotype(geno, subPop)'
        pop = population(loci=[1, 2], size=[1, 2])
        self.assertRaises(exceptions.IndexError, pop.setGenotype, [1], 2)
        pop.setGenotype([1, 2, 3])
        self.assertEqual(pop.individual(0).genotype(), [1, 2, 3, 1, 2, 3])
        self.assertEqual(pop.individual(1).genotype(), [1, 2, 3, 1, 2, 3])
        self.assertEqual(pop.individual(2).genotype(0), [1, 2, 3])
        self.assertEqual(pop.individual(2).genotype(1), [1, 2, 3])
        pop.setGenotype([2, 4], 1)
        self.assertEqual(pop.individual(0).genotype(), [1, 2, 3, 1, 2, 3])
        self.assertEqual(pop.individual(1).genotype(0), [2, 4, 2])
        self.assertEqual(pop.individual(1).genotype(1), [4, 2, 4])
        self.assertEqual(pop.individual(2).genotype(), [2, 4, 2, 4, 2, 4])

    def testAncestor(self):
        'Testing population::ancestor(idx, gen), ancestor(idx, subPop, gen)'
        pop = population([100, 200], loci=[10, 20], infoFields=['x', 'y'],
            ancGen=5)
        InitByFreq(pop, [0.2, 0.8])
        for ind in pop.individuals():
            ind.setInfo(random.randint(4, 10), 'x')
            ind.setInfo(random.randint(10, 100), 'y')
        pop1 = population([200, 100], loci=[10, 20], infoFields=['x', 'y'])
        InitByFreq(pop1, [0.5, 0.5])
        for ind in pop1.individuals():
            ind.setInfo(random.randint(4, 10), 'x')
            ind.setInfo(random.randint(10, 100), 'y')
        pop_c = pop.clone()
        pop.push(pop1)
        for idx, ind in enumerate(pop_c.individuals()):
            self.assertEqual(ind, pop.ancestor(idx, 1))
            self.assertEqual(ind.info('x'), pop.ancestor(idx, 1).info('x'))
            self.assertEqual(ind.info('y'), pop.ancestor(idx, 1).info('y'))
        self.assertRaises(exceptions.IndexError, pop.ancestor, 10000, 2)
        self.assertRaises(exceptions.IndexError, pop.ancestor, 10000, 3)
        for idx, ind in enumerate(pop_c.individuals(0)):
            self.assertEqual(ind, pop.ancestor(idx, 0, 1))

    def testAncestralGens(self):
        'Testing population::ancestralGens(), setAncestralDepth(depth), useAncestralGen(idx)'
        pop = population(size=[3, 5], loci=[2, 3],  infoFields=['x'])
        InitByFreq(pop, [.2, .8])
        pop.setIndInfo([random.random() for x in range(8)], 'x')
        pop.setAncestralDepth(-1)
        gt = list(pop.genotype())
        inf = pop.indInfo('x')
        self.assertEqual(pop.ancestralGens(), 0)
        pop1 = population(size=[2, 3], loci=[2, 3], ancGen=2,  infoFields=['x'])
        InitByFreq(pop1, [.8, .2])
        pop1.setIndInfo([random.random() for x in range(8)], 'x')
        gt1 = list(pop1.genotype())
        inf1 = pop1.indInfo('x')
        pop.push(pop1)
        self.assertEqual(pop.ancestralGens(), 1)
        self.assertEqual(pop.genotype(), gt1)
        self.assertEqual(pop.indInfo('x'), inf1)
        # subPopSize, indInfo
        self.assertEqual(pop.subPopSize(0), 2)
        self.assertEqual(pop.subPopSize(1), 3)
        pop.useAncestralGen(1)
        self.assertEqual(pop.genotype(), gt)
        self.assertEqual(pop.indInfo('x'), inf)
        pop.useAncestralGen(0)
        self.assertEqual(pop.genotype(), gt1)
        self.assertEqual(pop.indInfo('x'), inf1)
        pop2 = population(size=[3, 5], loci=[2, 3], infoFields=['x'])
        pop2.setIndInfo([random.random() for x in range(8)], 'x')
        inf2 = pop2.indInfo('x')
        InitByFreq(pop2, [.2, .8])
        gt2 = list(pop2.genotype())
        pop.push(pop2)
        self.assertEqual(pop.ancestralGens(), 2)
        self.assertEqual(pop.genotype(), gt2)
        self.assertEqual(pop.indInfo('x'), inf2)
        pop.useAncestralGen(1)
        self.assertEqual(pop.genotype(), gt1)
        self.assertEqual(pop.indInfo('x'), inf1)
        pop.useAncestralGen(2)
        self.assertEqual(pop.genotype(), gt)
        self.assertEqual(pop.indInfo('x'), inf)
        # out of bound ancestral generation number
        self.assertRaises(exceptions.ValueError, pop.useAncestralGen, 3 )

    def testAddChrom(self):
        'Testing population::addChrom'
        pop = self.getPop(chromNames=['c1', 'c2'], lociPos=[[1], [3, 5]], lociNames = ['l1', 'l2', 'l3'], ancGen=5)
        pop1 = pop.clone()
        pop.addChrom([7, 8, 9], ['l4', 'l5', 'l6'], 'c3')
        self.assertEqual(pop.numChrom(), 3)
        self.assertEqual(pop.chromNames(), ('c1', 'c2', 'c3'))
        self.assertEqual(pop.numLoci(), (1, 2, 3))
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            for idx in range(pop.popSize()):
                ind = pop.individual(idx)
                ind1 = pop1.individual(idx)
                for ch in range(2):
                    self.assertEqual(ind.genotype(0, ch), ind1.genotype(0, ch))
                    self.assertEqual(ind.genotype(1, ch), ind1.genotype(1, ch))
                    # new chromosome has zero values
                    self.assertEqual(ind.genotype(0, 2), [0]*3)  # new alleles are zero
                    self.assertEqual(ind.genotype(1, 2), [0]*3)  # new alleles are zero
        # lociPos is not ordered
        self.assertRaises(exceptions.ValueError,  pop.addChrom, [13, 12, 11], ['l4', 'l5', 'l6'], 'c3')
        # given loci names are not unique.
        self.assertRaises(exceptions.ValueError,  pop.addChrom, [11, 12, 13], ['l4', 'l5', 'l6'], 'c3')
        # # given chromsome name is not unique.
        self.assertRaises(exceptions.ValueError,  pop.addChrom, [11, 12, 13], ['l4', 'l5', 'l6'], 'c4')

    def testAddChromFromPop(self):
        'Testing population::addChromFromPop(pop)'
        pop = population(size=100, ploidy=2, loci=[1, 2], chromNames=["c1", "c2"], lociNames = ['l1', 'l2', 'l3'])
        pop2 = pop.clone()
        pop1 = population(size=100, ploidy=2, loci=[2, 3], chromNames=["c3", "c4"],
            lociNames = ['l4', 'l5', 'l6', 'l7', 'l8'])
        pop.addChromFromPop(pop1)
        self.assertEqual(pop.numChrom(), 4)
        self.assertEqual(pop.chromNames(), ('c1', 'c2', 'c3', 'c4'))
        self.assertEqual(pop.numLoci(), (1, 2, 2, 3))
        for i in range(100):
            ind = pop.individual(i)
            ind1 = pop2.individual(i)
            ind2 = pop1.individual(i)
            for loc in range(3):
                self.assertEqual(ind.allele(loc), ind1.allele(loc))
            for loc in range(5):
                self.assertEqual(ind.allele(loc+3), ind2.allele(loc))
        pop = population(size=100, ploidy=2, loci=[1, 2])
        pop1 = population(size=100, ploidy =2, loci=[1, 2])
        self.assertRaises(exceptions.ValueError, pop.addChromFromPop, pop1)
        pop2 = population(size=200, ploidy=2, loci=[2, 3], chromNames=["c3", "c4"],
            lociNames = ['l4', 'l5', 'l6', 'l7', 'l8'])
        self.assertRaises(exceptions.ValueError, pop.addChromFromPop, pop2)

    def testAddIndFromPop(self):
        'Testing population::addIndFromPop(pop)'
        pop = self.getPop(ancGen=3)
        pop1 = self.getPop(ancGen=3)
        pop.setIndInfo([random.randint(4, 10) for x in range(pop.popSize())], 'x')
        pop.addIndFromPop(pop1)
        self.assertEqual(pop.numSubPop(), 4)
        self.assertEqual(pop.subPopSizes(), (20, 80, 20, 80))
        for i in range(100):
            self.assertEqual(pop.individual(100+i), pop1.individual(i))
        pop1 = self.getPop(ancGen=2)
        # different numbers of ancestral generations
        self.assertRaises(exceptions.ValueError, pop.addIndFromPop, pop1)
        pop1 = population(size=100, ploidy=2, loci=[1, 2, 3])
        # different genotype structure
        self.assertRaises(exceptions.ValueError, pop.addIndFromPop, pop1)
        # Test scrambled populations
        pop = self.getPop(scramble=True, ancGen=3)
        pop1 = self.getPop(scramble=True, ancGen=3)
        pop.addIndFromPop(pop1)
        for i in range(100):
            self.assertEqual(pop.individual(100+i), pop1.individual(i))

    def testAddLociFromPop(self):
        'Testing population::addLociFromPop(pop)'
        pop = self.getPop(chromNames=["c1", "c2"], ancGen=5, lociPos=[[1], [2, 5]], lociNames = ['l1', 'l2', 'l3'])
        pop1 = pop.clone()
        pop2 = self.getPop(chromNames=["c3", "c4"], ancGen=5, lociPos=[[4], [3, 6]], lociNames = ['l4', 'l5', 'l6'])
        pop.addLociFromPop(pop2);
        self.assertEqual(pop.numLoci(), (2, 4))
        self.assertEqual(pop.lociPos(), (1, 4, 2, 3, 5, 6))
        self.assertEqual(pop.chromNames(), ('c1', 'c2'))
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            pop2.useAncestralGen(gen)
            for idx in range(pop.popSize()):
                ind = pop.individual(idx)
                inds = [pop1.individual(idx), pop2.individual(idx)]
                # i: index in population
                # src: the source population
                # j: index in source population
                for i, src, j in [(0, 0, 0), (1, 1, 0), (2, 0, 1), (3, 1, 1), (4, 0, 2), (5, 1, 2)]:
                    for p in range(pop.ploidy()):
                        self.assertEqual(ind.allele(i, p), inds[src].allele(j, p))


    def testAddLoci(self):
        'Testing population::addLoci(chrom, pos, names=[])'
        pop = self.getPop(chromNames=["c1", "c2"], ancGen=5, lociPos=[[1], [3, 5]], lociNames = ['l1', 'l2', 'l3'])
        pop1 = pop.clone()
        pop.addLoci([1], [6], ['l4'])
        self.assertEqual(pop.numLoci(), (1, 3))
        self.assertEqual(pop.lociPos(), (1, 3, 5, 6))
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            for idx in range(pop.popSize()):
                ind = pop.individual(idx)
                ind1 = pop1.individual(idx)
                for loc in range(3):
                    self.assertEqual(ind.allele(loc), ind1.allele(loc))
                self.assertEqual(ind.allele(3), 0)
        self.assertRaises(exceptions.ValueError, pop.addLoci, [2], [6], ['l4'])

    def testDeepcopy(self):
        'Testing deepcopy of population'
        pop = self.getPop(True, False, ancGen=3)
        InitByFreq(pop, [0.2, 0.8])
        # shallow copy
        pop1 = pop
        InitByFreq(pop1, [0.8, 0.2])
        self.assertEqual(pop, pop1)
        # deep copy
        pop1 = pop.clone()
        self.assertEqual(pop, pop1)
        InitByFreq(pop1, [0.5, 0.5])
        self.assertNotEqual(pop, pop1)
        # using Python copy.copy
        pop1 = copy.copy(pop)
        self.assertEqual(pop, pop1)
        InitByFreq(pop1, [0.5, 0.5])
        self.assertEqual(pop, pop1)
        # using Python copy.deepcopy
        pop1 = copy.deepcopy(pop)
        self.assertEqual(pop, pop1)
        InitByFreq(pop1, [0.5, 0.5])
        self.assertNotEqual(pop, pop1)

    def testExtract(self):
        'Testing population::Extract(field=None, loci=None, infoFields=None, ancGen =-1)'
        pop = population(size=[3, 5], loci=[2, 3],  infoFields=['x', 'y'])
        for ind in pop.individuals():
            n = random.randint(-1, 5)
            ind.setInfo(n, 'x')
            ind.setInfo(n + 10, 'y')
            ind.setGenotype([n+1])
        pop1 = pop.extract(field='x')
        for sp in range(5):
            for ind in pop1.individuals(sp):
                self.assertEqual(ind.info('x'), sp)
                self.assertEqual(ind.info('y'), sp + 10)
                self.assertEqual(ind.genotype(), [sp+1]*(pop.totNumLoci()*pop.ploidy()))

    def testMergeSubPops(self):
        'Testing population::mergeSubPops(subpops=[])'
        pop = self.getPop(size=[100, 20, 30, 80, 50, 60])
        pop1 = pop.clone()
        pop.mergeSubPops([1, 2, 4])
        self.assertEqual(pop.subPopSize(1), pop1.subPopSize(1)+pop1.subPopSize(2)+pop1.subPopSize(4))
        for (oldsp, newsp) in [(0, 0), (3, 2), (5, 3)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop.subPopSize(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop.individual(idx, newsp))


    def testRemoveSubPops(self):
        'Testing population::removeEmptySubPops(), removeSubPops()'
        pop = self.getPop(size=[0, 100, 0, 20, 30, 0, 50])
        pop1 = pop.clone()
        self.assertEqual( pop.numSubPop(), 7)
        pop.removeSubPops([x for x in range(7) if pop.subPopSize(x) == 0])
        self.assertEqual( pop.numSubPop(), 4)
        self.assertEqual( pop.subPopSizes(), (100, 20, 30, 50))
        for (oldsp, newsp) in [(1, 0), (3, 1), (4, 2), (6, 3)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop.subPopSize(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop.individual(idx, newsp))
        # remove subpop
        pop2 = pop.clone()
        pop.removeSubPops([1, 2])
        self.assertEqual( pop.subPopSizes(), (100, 50))
        for (oldsp, newsp) in [(0, 0), (3, 1)]:  # map of old and new id.
            self.assertEqual(pop2.subPopSize(oldsp), pop.subPopSize(newsp))
            for idx in range(pop2.subPopSize(oldsp)):
                self.assertEqual(pop2.individual(idx, oldsp), pop.individual(idx, newsp))
        #pop.removeSubPops([8])


    def testRemoveIndividuals(self):
        'Testing population::removeIndividuals(inds)'
        pop = self.getPop(size =[20, 100, 30])
        pop1 = pop.clone()
        pop.removeIndividuals([15])
        self.assertEqual(pop.subPopSizes(), (19, 100, 30))
        for idx in range(15):
            self.assertEqual(pop1.individual(idx), pop.individual(idx))
        for idx in range(15, pop.popSize()):
            self.assertEqual(pop1.individual(idx+1), pop.individual(idx))


    def testRemoveLoci(self):
        'Testing population::removeLoci(loci=[], keep=[])'
        pop = self.getPop(size=[1, 2], loci=[2, 3, 1], ancGen=5)
        pop.removeLoci([2])
        pop1 = pop.clone()
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            for idx in range(pop.popSize()):
                ind = pop.individual(idx)
                ind1 = pop1.individual(idx)
                for loc in range(2):
                    self.assertEqual(ind.allele(loc), ind1.allele(loc))
                for loc in range(2, 5):
                    self.assertEqual(ind.allele(loc), ind1.allele(loc+1))



if __name__ == '__main__':
    unittest.main()


