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
    def getPop(self, scramble=False, VSP=False):
        pop = population(size=[20, 80], ploidy=2, loci=[5, 7],
            lociPos=[ [2, 3, 4, 5, 6], [2, 4, 6, 8, 10, 12, 14]],
            alleleNames=['_', 'A', 'C', 'T', 'G'],
            infoFields=['a', 'b'])
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
        'Testing population::numVirtualSubPop(), setVirtualSplitter(splitter), virtualSubPopName(subPop)'
        pop = population(1000, infoFields=['x'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(infoSplitter('x', values=range(10, 15)))
        self.assertEqual(pop.numVirtualSubPop(), 5)
        infos = list(pop.indInfo('x'))
        self.assertEqual(pop.virtualSubPopName(0), "x = 10")
        self.assertEqual(pop.virtualSubPopName(1), "x = 11")
        self.assertEqual(pop.virtualSubPopName(4), "x = 14")

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
        pop = population(1000, infoFields=['x'])
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
            ancestralDepth=5)
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
        pop.pushAndDiscard(pop1)
        for idx, ind in enumerate(pop_c.individuals()):
            self.assertEqual(ind, pop.ancestor(idx, 1))
            self.assertEqual(ind.info('x'), pop.ancestor(idx, 1).info('x'))
            self.assertEqual(ind.info('y'), pop.ancestor(idx, 1).info('y'))
        self.assertRaises(exceptions.IndexError, pop.ancestor, 10000, 2)
        self.assertRaises(exceptions.IndexError, pop.ancestor, 10000, 3)
        for idx, ind in enumerate(pop_c.individuals(0)):
            self.assertEqual(ind, pop.ancestor(idx, 0, 1))

    def testDeepcopy(self):
        'Testing deepcopy of population'
        pop = population(10, loci=[2])
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

if __name__ == '__main__':
    unittest.main()
