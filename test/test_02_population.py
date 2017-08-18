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
  
import unittest, os, sys, random, copy
from simuOpt import setOptions
setOptions(quiet=True)
new_argv = []
for arg in sys.argv:
    if arg in ['short', 'long', 'binary', 'mutant', 'lineage']:
        setOptions(alleleType = arg)
    elif arg.startswith('-j'):
        setOptions(numThreads = int(arg[2:]))
    else:
        new_argv.append(arg) 

sys.argv=new_argv
from simuPOP import *

class TestPopulation(unittest.TestCase):
    # define a few functions to create basic populations
    def getPop(self, VSP=False, size=[20, 80], loci = [1, 2], infoFields=['x'],
            ancGen=0, *arg, **kwargs):
        pop = Population(size=size, ploidy=2, loci=loci, infoFields=infoFields,
            ancGen=ancGen, *arg, **kwargs)
        pop.setGenotype([random.randint(1, 5) for x in range(pop.popSize()*pop.ploidy())])
        for info in infoFields:
            pop.setIndInfo([random.random() for x in range(pop.popSize())], info)
        for i in range(ancGen):
            pop.push(self.getPop(size=size, loci=loci, infoFields=infoFields, ancGen=0, *arg, **kwargs))
        initSex(pop)
        if VSP:
            pop.setVirtualSplitter(SexSplitter())
        return pop

    def testAbsIndIndex(self):
        'Testing Population::absIndIndex(idx, subPop), popSize()'
        pop = self.getPop()
        # ind, subPop
        self.assertEqual(pop.absIndIndex(1, 1), 21)
        self.assertEqual(pop.absIndIndex(10, 0), 10)
        self.assertRaises(IndexError, pop.absIndIndex, 0, 2 )
        self.assertEqual(pop.popSize(), 100)

    def testSubPop(self):
        'Testing Population::subPopBegin(subPop), subPopEnd(subPop), numSubPop()'
        'subPopSize(subPop), subPopSizes(), subPopIndPair(idx)'
        pop = self.getPop()
        self.assertEqual(pop.subPopBegin(1), 20)
        self.assertRaises(IndexError, pop.subPopBegin, 2 )
        self.assertEqual(pop.subPopEnd(0), 20)
        self.assertRaises(IndexError, pop.subPopEnd, 2 )
        self.assertEqual(pop.numSubPop(), 2)
        self.assertEqual(pop.subPopSize(0), 20)
        self.assertEqual(pop.subPopSize(1), 80)
        self.assertRaises(IndexError, pop.subPopSize, 2 )
        self.assertEqual(pop.subPopSizes(), (20, 80) )
        self.assertEqual(pop.subPopIndPair(21), (1, 1) )
        self.assertRaises(IndexError, pop.subPopIndPair, 200 )
        stat(pop, numOfMales=True, vars=['numOfMales_sp', 'numOfFemales_sp'])
        pop.setVirtualSplitter(SexSplitter())
        self.assertEqual(pop.subPopSize([1, 0]), pop.dvars(1).numOfMales)
        self.assertEqual(pop.subPopSize([1, 1]), pop.dvars(1).numOfFemales)

    def testVirtualSubPop(self):
        'Testing Population::numVirtualSubPop(), setVirtualSplitter(splitter), subPopName(subPop)'
        pop = Population(1000, infoFields=['x'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(InfoSplitter('x', values=list(range(10, 15))))
        self.assertEqual(pop.numVirtualSubPop(), 5)
        self.assertEqual(pop.subPopName(0), "")
        self.assertEqual(pop.subPopName([0, 0]), "x = 10")
        self.assertEqual(pop.subPopName([0, 1]), "x = 11")
        self.assertEqual(pop.subPopName([0, 4]), "x = 14")
        self.assertRaises(IndexError, pop.subPopName, 1)
        self.assertRaises(IndexError, pop.subPopName, [0, 5])
        # this function accepts vsp name
        pop.individuals([0, 'x = 10'])
        # with given names
        pop = Population(size=[200, 500], infoFields=['x'], subPopNames=['A', 'B'])
        for ind in pop.individuals(): ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(InfoSplitter('x', values=list(range(10, 15))))
        self.assertEqual(pop.numVirtualSubPop(), 5)
        self.assertEqual(pop.subPopName(0), "A")
        self.assertEqual(pop.subPopName(1), "B")
        self.assertRaises(IndexError, pop.subPopName, 2)
        self.assertEqual(pop.subPopName([0, 0]), "A - x = 10")
        self.assertEqual(pop.subPopName([0, 1]), "A - x = 11")
        self.assertEqual(pop.subPopName([1, 4]), "B - x = 14")
        self.assertRaises(IndexError, pop.subPopName, [0, 5])

    def testPopSize(self):
        'Testing Population.popSize by male, female, pair'
        pop = self.getPop(size=[80, 20, 30, 50], ancGen=5)
        pop.mergeSubPops([1, 2])
        pop.removeIndividuals([1, 59, 130])
        initSex(pop)
        pop.setVirtualSplitter(SexSplitter())
        self.assertEqual(pop.popSize(), 177)
        self.assertEqual(pop.popSize(1), 180)
        self.assertEqual(pop.subPopSizes(), (78, 50, 49))
        self.assertEqual(pop.subPopSize(1), 50)
        self.assertEqual(pop.subPopSizes(1), (80, 20, 30, 50))
        self.assertEqual(pop.subPopSize(1, 1), 20)
        pop.useAncestralGen(2)
        self.assertEqual(pop.popSize(), 180)
        self.assertEqual(pop.subPopSizes(), (80, 20, 30, 50))
        self.assertEqual(pop.popSize(2), 180)
        self.assertEqual(pop.subPopSizes(2), (80, 20, 30, 50))
        self.assertEqual(pop.popSize(3), 180)
        self.assertEqual(pop.subPopSizes(3), (80, 20, 30, 50))
        self.assertEqual(pop.subPopSize(1, 3), 20)
        # virtual
        stat(pop, numOfMales=True, subPops=[1])
        self.assertEqual(pop.subPopSize((1,0)), pop.dvars().numOfMales)
        pop.useAncestralGen(0)
        self.assertEqual(pop.subPopSize((1,0), ancGen=2), pop.dvars().numOfMales)

    def testPopSizeBySex(self):
        'Testing Population.popSize, Population.subPopSizes and Population.subPopSize'
        pop = self.getPop(size=[90, 10, 30, 50], ancGen=5)
        initSex(pop, sex=[MALE, FEMALE, MALE])
        pop.useAncestralGen(2)
        initSex(pop, sex=[MALE, FEMALE, FEMALE])
        pop.useAncestralGen(0)
        self.assertEqual(pop.popSize(sex=MALE_ONLY), 120)
        self.assertEqual(pop.popSize(sex=FEMALE_ONLY), 60)
        self.assertEqual(pop.popSize(sex=PAIR_ONLY), 60)
        self.assertEqual(pop.popSize(2, sex=MALE_ONLY), 60)
        self.assertEqual(pop.popSize(2, sex=FEMALE_ONLY), 120)
        self.assertEqual(pop.popSize(2, sex=PAIR_ONLY), 60)
        # VSP
        pop.setVirtualSplitter(SexSplitter())
        self.assertEqual(pop.subPopSize(2, sex=MALE_ONLY), 20)
        self.assertEqual(pop.subPopSize(2, sex=FEMALE_ONLY), 10)
        self.assertEqual(pop.subPopSize(2, sex=PAIR_ONLY), 10)
        self.assertEqual(pop.subPopSize(2, 2, sex=MALE_ONLY), 10)
        self.assertEqual(pop.subPopSize(2, 2, sex=FEMALE_ONLY), 20)
        self.assertEqual(pop.subPopSize(2, 2, sex=PAIR_ONLY), 10)
        #
        self.assertEqual(pop.subPopSize((2, 'Male'), sex=MALE_ONLY), 20)
        self.assertEqual(pop.subPopSize((2, 'Male'), sex=FEMALE_ONLY), 0)
        self.assertEqual(pop.subPopSize((2, 'Male'), sex=PAIR_ONLY), 0)
        self.assertEqual(pop.subPopSize((0, 'Female'), 2, sex=MALE_ONLY), 0)
        self.assertEqual(pop.subPopSize((0, 'Female'), 2, sex=FEMALE_ONLY), 60)
        self.assertEqual(pop.subPopSize((0, 'Female'), 2, sex=PAIR_ONLY), 0)


    def testLociPos(self):
        'Testing lociPos parameter of Population::Population'
        # test for a bug that condier the following two numbers are the same.
        Population(loci=2, lociPos=[29.114998502, 29.114998525])

    def testSubPopName(self):
        'Testing Population::setSubPopName(name, subPop), subPopByName(subPop)'
        pop = self.getPop(size=[80, 20, 30, 50], ancGen=5)
        pop.setSubPopName('A', 0)
        pop.setSubPopName('B', 1)
        pop.setSubPopName('B', 2)
        pop.setSubPopName('C', 3)
        self.assertEqual(pop.subPopName(0), 'A')
        self.assertEqual(pop.subPopName(1), 'B')
        self.assertEqual(pop.subPopName(2), 'B')
        self.assertEqual(pop.subPopName(3), 'C')
        self.assertEqual(pop.subPopByName('A'), 0)
        self.assertEqual(pop.subPopByName('B'), 1)
        self.assertEqual(pop.subPopByName('C'), 3)
        self.assertRaises(ValueError, pop.subPopByName, 'D')

    def testIndividuals(self):
        'Testing function Population::individuals(), individuals(subPop), individual(idx, subPop=0)'
        def testAllInd(pop):
            self.assertEqual(len(list(pop.individuals())), pop.popSize())
            self.assertEqual(len(list(pop.individuals(0))), pop.subPopSize(0))
            self.assertEqual(len(list(pop.individuals(1))), pop.subPopSize(1))
        testAllInd(self.getPop())
        testAllInd(self.getPop(True))
        pop = Population([20, 80], loci = [5, 7], infoFields=['x'])
        pop.individual(0).setAllele(1, 0)
        self.assertEqual(pop.individual(0).allele(0), 1)

    def testGenotype(self):
        'Testing Population::genotype(), genotype(subPop)'
        pop = Population(loci=[1, 2], size=[1, 2])
        arr = pop.genotype()
        self.assertEqual(len(arr), pop.genoSize()*pop.popSize())
        arr = pop.genotype(1)
        self.assertEqual(len(arr), pop.genoSize()*pop.subPopSize(1))
        self.assertRaises(IndexError, pop.genotype, 2)



    def testSetGenotype(self):
        'Testing Population::setGenotype(geno), setGenotype(geno, subPop)'
        pop = Population(loci=[1, 2], size=[1, 2])
        self.assertRaises(IndexError, pop.setGenotype, [1], 2)
        if moduleInfo()['alleleType'] == 'binary':
            pop.setGenotype([0, 1, 0])
            self.assertEqual(pop.individual(0).genotype(), [0, 1, 0, 0, 1, 0])
            self.assertEqual(pop.individual(1).genotype(), [0, 1, 0, 0, 1, 0])
            self.assertEqual(pop.individual(2).genotype(0), [0, 1, 0])
            self.assertEqual(pop.individual(2).genotype(1), [0, 1, 0])
            pop.setGenotype([1, 0], 1)
            self.assertEqual(pop.individual(0).genotype(), [0, 1, 0, 0, 1, 0])
            self.assertEqual(pop.individual(1).genotype(0), [1, 0, 1])
            self.assertEqual(pop.individual(1).genotype(1), [0, 1, 0])
            self.assertEqual(pop.individual(2).genotype(), [1, 0, 1, 0, 1, 0])
            # virtual subpopulation
            pop = self.getPop(size = 100, VSP=True)
            self.assertEqual(pop.numSubPop(), 1)
            self.assertEqual(pop.numVirtualSubPop(), 2)
            pop.setGenotype([5], [0, 0])
            pop.setGenotype([6], [0, 1])
            for idx, ind in enumerate(pop.individuals([0, 0])):
                self.assertEqual(ind.allele(idx%6), 1)
            for idx, ind in enumerate(pop.individuals([0, 1])):
                self.assertEqual(ind.allele(idx%6), 1)
        else:
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
            # virtual subpopulation
            pop = self.getPop(size = 100, VSP=True)
            self.assertEqual(pop.numSubPop(), 1)
            self.assertEqual(pop.numVirtualSubPop(), 2)
            pop.setGenotype([5], [0, 0])
            pop.setGenotype([6], [0, 1])
            for idx, ind in enumerate(pop.individuals([0, 0])):
                self.assertEqual(ind.allele(idx%6), 5)
            for idx, ind in enumerate(pop.individuals([0, 1])):
                self.assertEqual(ind.allele(idx%6), 6)

    def testAncestor(self):
        'Testing Population::ancestor(idx, gen), ancestor(idx, gen, subPop), push(pop)'
        pop = Population([100, 200], loci=[10, 20], infoFields=['x', 'y'],
            ancGen=5)
        initSex(pop)
        initGenotype(pop, freq=[0.2, 0.8])
        for ind in pop.individuals():
            ind.setInfo(random.randint(4, 10), 'x')
            ind.setInfo(random.randint(10, 100), 'y')
        pop1 = Population([200, 100], loci=[10, 20], infoFields=['x', 'y'])
        initSex(pop1)
        initGenotype(pop1, freq= [0.5, 0.5])
        for ind in pop1.individuals():
            ind.setInfo(random.randint(4, 10), 'x')
            ind.setInfo(random.randint(10, 100), 'y')
        pop_c = pop.clone()
        pop.push(pop1)
        for idx, ind in enumerate(pop_c.individuals()):
            self.assertEqual(ind, pop.ancestor(idx, 1))
            self.assertEqual(ind.info('x'), pop.ancestor(idx, 1).info('x'))
            self.assertEqual(ind.info('y'), pop.ancestor(idx, 1).info('y'))
        self.assertRaises(IndexError, pop.ancestor, 2, 10000)
        self.assertRaises(IndexError, pop.ancestor, 3, 10000)
        for idx, ind in enumerate(pop_c.individuals(0)):
            self.assertEqual(ind, pop.ancestor(idx, 1, 0))

    def testAncestralGens(self):
        'Testing Population::ancestralGens(), setAncestralDepth(depth), useAncestralGen(idx)'
        pop = Population(size=[3, 5], loci=[2, 3],  infoFields=['x'])
        initSex(pop)
        initGenotype(pop, freq=[.2, .8])
        pop.setIndInfo([random.random() for x in range(8)], 'x')
        pop.setAncestralDepth(-1)
        gt = list(pop.genotype())
        inf = pop.indInfo('x')
        self.assertEqual(pop.ancestralGens(), 0)
        pop1 = Population(size=[2, 3], loci=[2, 3], ancGen=2,  infoFields=['x'])
        initSex(pop1)
        initGenotype(pop1, freq= [.8, .2])
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
        pop2 = Population(size=[3, 5], loci=[2, 3], infoFields=['x'])
        pop2.setIndInfo([random.random() for x in range(8)], 'x')
        inf2 = pop2.indInfo('x')
        initSex(pop2)
        initGenotype(pop2, freq= [.2, .8])
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
        self.assertRaises(ValueError, pop.useAncestralGen, 3 )
        # setAncestralDepth
        pop = self.getPop(ancGen = 5)
        pop.setAncestralDepth(3)
        self.assertEqual(pop.ancestralGens(), 3)

    def testAddChrom(self):
        'Testing Population::addChrom'
        pop = self.getPop(chromNames=['c1', 'c2'], lociPos=[1, 3, 5], lociNames = ['l1', 'l2', 'l3'], ancGen=5)
        pop1 = pop.clone()
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            pop.setLineage(1)
            pop1.setLineage(2)
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
                    #
                    if moduleInfo()['alleleType'] == 'lineage':
                        self.assertEqual(ind.lineage(0, ch), [1]*ind.numLoci(ch))
                        self.assertEqual(ind.lineage(0, ch), [1]*ind.numLoci(ch))
                        self.assertEqual(ind.lineage(0, 2), [0]*3)
                        self.assertEqual(ind.lineage(1, 2), [0]*3)
        # lociPos is not ordered
        self.assertRaises(ValueError,  pop.addChrom, [13, 12, 11], ['l4', 'l5', 'l6'], 'c3')
        # given loci names are not unique.
        self.assertRaises(ValueError,  pop.addChrom, [11, 12, 13], ['l4', 'l5', 'l6'], 'c3')
        # # given chromsome name is not unique.
        self.assertRaises(ValueError,  pop.addChrom, [11, 12, 13], ['l4', 'l5', 'l6'], 'c4')

    def testAddChromFrom(self):
        'Testing Population::addChromFrom(pop)'
        pop = Population(size=100, ploidy=2, loci=[1, 2], chromNames=["c1", "c2"], lociNames = ['l1', 'l2', 'l3'])
        pop2 = pop.clone()
        pop1 = Population(size=100, ploidy=2, loci=[2, 3], chromNames=["c3", "c4"],
            lociNames = ['l4', 'l5', 'l6', 'l7', 'l8'])
        pop.addChromFrom(pop1)
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
        pop = Population(size=100, ploidy=2, loci=[1, 2])
        pop1 = Population(size=200, ploidy=2, loci=[2, 3], chromNames=["c3", "c4"],
            lociNames = ['l4', 'l5', 'l6', 'l7', 'l8'])
        # population size is different
        self.assertRaises(ValueError, pop.addChromFrom, pop1)
        # see what happens to alleleNames
        pop1 = Population(size=100, ploidy=2, loci=[1, 2], chromNames=["c1", "c2"],
            lociNames = ['l1', 'l2', 'l3'], alleleNames=['A', 'B'])
        pop2 = Population(size=100, ploidy=2, loci=[2, 3], chromNames=["c3", "c4"],
            lociNames = ['l4', 'l5', 'l6', 'l7', 'l8'])
        pop1.addChromFrom(pop2)
        self.assertEqual(pop1.alleleNames(0), ('A', 'B'))
        self.assertEqual(pop1.alleleNames(2), ('A', 'B'))
        self.assertEqual(pop1.alleleNames(3), ())
        self.assertEqual(pop1.alleleNames(4), ())
        self.assertEqual(pop1.alleleNames(7), ())
        #
        pop1 = Population(size=100, ploidy=2, loci=[1, 2], chromNames=["c1", "c2"],
            lociNames = ['l1', 'l2', 'l3'], alleleNames=['A', 'B'])
        pop2 = Population(size=100, ploidy=2, loci=[2], chromNames=["c3"],
            lociNames = ['l4', 'l5'],
            alleleNames=[['E', 'F'], ['C', 'D']])
        pop1.addChromFrom(pop2)
        self.assertEqual(pop1.alleleNames(0), ('A', 'B'))
        self.assertEqual(pop1.alleleNames(2), ('A', 'B'))
        self.assertEqual(pop1.alleleNames(3), ('E', 'F'))
        self.assertEqual(pop1.alleleNames(4), ('C', 'D'))

    def testAddIndFrom(self):
        'Testing Population::addIndFrom(pop)'
        pop = self.getPop(ancGen=3)
        pop1 = self.getPop(ancGen=3)
        pop.setIndInfo([random.randint(4, 10) for x in range(pop.popSize())], 'x')
        pop.addIndFrom(pop1)
        self.assertEqual(pop.numSubPop(), 4)
        self.assertEqual(pop.subPopSizes(), (20, 80, 20, 80))
        for i in range(100):
            self.assertEqual(pop.individual(100+i), pop1.individual(i))
        pop1 = self.getPop(ancGen=2)
        # different numbers of ancestral generations
        self.assertRaises(ValueError, pop.addIndFrom, pop1)
        pop1 = Population(size=100, ploidy=2, loci=[1, 2, 3])
        # different genotype structure
        self.assertRaises(ValueError, pop.addIndFrom, pop1)

    def testAddLociFrom(self):
        'Testing Population::addLociFrom(pop)'
        pop = self.getPop(chromNames=["c1", "c2"], ancGen=5, lociPos=[1, 2, 5], lociNames = ['l1', 'l2', 'l3'])
        pop1 = pop.clone()
        pop2 = self.getPop(chromNames=["c3", "c4"], ancGen=5, lociPos=[4, 3, 6], lociNames = ['l4', 'l5', 'l6'])
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            pop2.useAncestralGen(gen)
            pop.setLineage(1)
            pop1.setLineage(2)
            pop2.setLineage(3)
        pop.addLociFrom(pop2);
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
                # src: the source Population
                # j: index in source Population
                for i, src, j in [(0, 0, 0), (1, 1, 0), (2, 0, 1), (3, 1, 1), (4, 0, 2), (5, 1, 2)]:
                    for p in range(pop.ploidy()):
                        self.assertEqual(ind.allele(i, p), inds[src].allele(j, p))
                for i in range(0, 2, 4):
                    for p in range(pop.ploidy()):
                        if moduleInfo()['alleleType'] == 'lineage':
                            self.assertEqual(ind.alleleLineage(i, p), 1)
                        else:
                            self.assertEqual(ind.alleleLineage(i, p), 0)
                for i in range(1, 3, 5):
                    for p in range(pop.ploidy()):
                        if moduleInfo()['alleleType'] == 'lineage':
                            self.assertEqual(ind.alleleLineage(i, p), 3)
                        else:
                            self.assertEqual(ind.alleleLineage(i, p), 0)
        # allele names
        pop = self.getPop(chromNames=["c1", "c2"], ancGen=5, lociPos=[1, 2, 5],
            lociNames = ['l1', 'l2', 'l3'], alleleNames=['A'])
        pop2 = self.getPop(chromNames=["c3", "c4"], ancGen=5, lociPos=[4, 3, 6],
            lociNames = ['l4', 'l5', 'l6'], alleleNames=[['B'], ['C', 'D'], ['E']])
        pop.addLociFrom(pop2);
        self.assertEqual(pop.alleleNames(0), ('A',))
        self.assertEqual(pop.alleleNames(1), ('B',))
        self.assertEqual(pop.alleleNames(2), ('A',))
        self.assertEqual(pop.alleleNames(3), ('C', 'D'))
        self.assertEqual(pop.alleleNames(4), ('A',))
        self.assertEqual(pop.alleleNames(5), ('E',))

    def testAddLociFromByName(self):
        'Testing Population::addLociFrom(pop, byName=True)'
        pop = self.getPop(chromNames=["c1", "c2"], ancGen=5, lociPos=[1, 2, 5],
            lociNames = ['l1', 'l2', 'l3'])
        pop1 = pop.clone()
        pop2 = self.getPop(chromNames=["c2", "c4"], ancGen=5, lociPos=[4, 3, 6],
            lociNames = ['l4', 'l5', 'l6'])
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            pop2.useAncestralGen(gen)
            pop.setLineage(1)
            pop1.setLineage(2)
            pop2.setLineage(3)
        pop.addLociFrom(pop2, byName=True);
        self.assertEqual(pop.numLoci(), (1, 3, 2))
        self.assertEqual(pop.lociPos(), (1, 2, 4, 5, 3, 6))
        self.assertEqual(pop.chromNames(), ('c1', 'c2', 'c4'))
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            pop2.useAncestralGen(gen)
            for idx in range(pop.popSize()):
                ind = pop.individual(idx)
                inds = [pop1.individual(idx), pop2.individual(idx)]
                # i: index in population
                # src: the source Population
                # j: index in source Population
                for i, src, j in [(0, 0, 0), (1, 0, 1), (2, 1, 0), (3, 0, 2), (4, 1, 1), (5, 1, 2)]:
                    for p in range(pop.ploidy()):
                        self.assertEqual(ind.allele(i, p), inds[src].allele(j, p))
                for i in range(0, 1, 3):
                    for p in range(pop.ploidy()):
                        if moduleInfo()['alleleType'] == 'lineage':
                            self.assertEqual(ind.alleleLineage(i, p), 1)
                        else:
                            self.assertEqual(ind.alleleLineage(i, p), 0)
                for i in range(2, 4, 5):
                    for p in range(pop.ploidy()):
                        if moduleInfo()['alleleType'] == 'lineage':
                            #???
                            self.assertEqual(ind.alleleLineage(i, p), 3)
                        else:
                            self.assertEqual(ind.alleleLineage(i, p), 0)
        # allele names
        pop = self.getPop(chromNames=["c1", "c2"], ancGen=5, lociPos=[1, 2, 5],
            lociNames = ['l1', 'l2', 'l3'], alleleNames=['A'])
        pop2 = self.getPop(chromNames=["c3", "c2"], ancGen=5, lociPos=[4, 3, 6],
            lociNames = ['l4', 'l5', 'l6'], alleleNames=[['B'], ['C', 'D'], ['E']])
        pop.addLociFrom(pop2, True);
        self.assertEqual(pop.alleleNames(0), ('A',))
        self.assertEqual(pop.alleleNames(1), ('A',))
        self.assertEqual(pop.alleleNames(2), ('C', 'D'))
        self.assertEqual(pop.alleleNames(3), ('A',))
        self.assertEqual(pop.alleleNames(4), ('E',))
        self.assertEqual(pop.alleleNames(5), ('B',))



    def testAddLoci(self):
        'Testing Population::addLoci(chrom, pos, names=[])'
        # special cases where the destination chromosome is empty
        # empty chromosome is the last
        pop = Population(size=5, ploidy = 2, loci = [1,0], lociPos = [0.4],
            chromNames = ['One', 'Two'], lociNames = ['selSite'])
        initGenotype(pop, freq=[0.5, 0.5])
        g1 = [x.allele(0, 0) for x in pop.individuals()] + [x.allele(0, 1) for x in pop.individuals()]
        pop.addLoci(chrom = 1, pos = 0.1, lociNames = '')
        self.assertEqual(pop.numChrom(), 2)
        self.assertEqual(pop.numLoci(0), 1)
        self.assertEqual(pop.numLoci(1), 1)
        self.assertEqual(pop.locusName(0), 'selSite')
        self.assertEqual(pop.locusName(1), '')
        g1_after = [x.allele(0, 0) for x in pop.individuals()] + [x.allele(0, 1) for x in pop.individuals()]
        self.assertEqual(g1, g1_after)
        g2 = [x.allele(1, 0) for x in pop.individuals()] + [x.allele(1, 1) for x in pop.individuals()]
        self.assertEqual(g2, [0]*10)
        #
        # empty chromosome is the first
        pop = Population(size=5, ploidy = 2, loci = [0,1], lociPos = [0.4],
            chromNames = ['One', 'Two'], lociNames = ['selSite'])
        initGenotype(pop, freq=[0.5, 0.5])
        g1 = [x.allele(0, 0) for x in pop.individuals()] + [x.allele(0, 1) for x in pop.individuals()]
        pop.addLoci(chrom = 0, pos = 0.1, lociNames = '')
        self.assertEqual(pop.numChrom(), 2)
        self.assertEqual(pop.numLoci(0), 1)
        self.assertEqual(pop.numLoci(1), 1)
        self.assertEqual(pop.locusName(0), '')
        self.assertEqual(pop.locusName(1), 'selSite')
        g1_after = [x.allele(1, 0) for x in pop.individuals()] + [x.allele(1, 1) for x in pop.individuals()]
        self.assertEqual(g1, g1_after)
        g2 = [x.allele(0, 0) for x in pop.individuals()] + [x.allele(0, 1) for x in pop.individuals()]
        self.assertEqual(g2, [0]*10)
        #
        # empty chromosome is in the middle
        pop = Population(size=5, ploidy = 2, loci = [1,0,0,1], lociPos = [0.4, 0.4],
            chromNames = ['One', 'Two', 'Three', 'Four'], lociNames = ['on_ch1', 'on_ch4'])
        initGenotype(pop, freq=[0.5, 0.5])
        g1 = [x.allele(0, 0) for x in pop.individuals()] + [x.allele(0, 1) for x in pop.individuals()]
        g2 = [x.allele(1, 0) for x in pop.individuals()] + [x.allele(1, 1) for x in pop.individuals()]
        pop.addLoci(chrom = 2, pos = 0.1, lociNames = 'new')
        self.assertEqual(pop.numChrom(), 4)
        self.assertEqual(pop.numLoci(0), 1)
        self.assertEqual(pop.numLoci(1), 0)
        self.assertEqual(pop.numLoci(2), 1)
        self.assertEqual(pop.numLoci(3), 1)
        self.assertEqual(pop.locusName(0), 'on_ch1')
        self.assertEqual(pop.locusName(1), 'new')
        self.assertEqual(pop.locusName(2), 'on_ch4')
        g1_after = [x.allele(0, 0) for x in pop.individuals()] + [x.allele(0, 1) for x in pop.individuals()]
        g2_after = [x.allele(2, 0) for x in pop.individuals()] + [x.allele(2, 1) for x in pop.individuals()]
        self.assertEqual(g1, g1_after)
        self.assertEqual(g2, g2_after)
        g3_after = [x.allele(1, 0) for x in pop.individuals()] + [x.allele(1, 1) for x in pop.individuals()]
        self.assertEqual(g3_after, [0]*10)
        #
        pop = self.getPop(size = 100, chromNames=["c1", "c2"], ancGen=5, lociPos=[1, 3, 5], lociNames = ['l1', 'l2', 'l3'])
        pop1 = pop.clone()
        newpos = pop.addLoci([0, 1, 1], [2, 6, 7], ['l4', 'l5', 'l6'])
        self.assertEqual(pop.numLoci(), (2, 4))
        self.assertEqual(pop.lociPos(), (1, 2, 3, 5, 6, 7))
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            for idx in range(pop.popSize()):
                ind = pop.individual(idx)
                ind1 = pop1.individual(idx)
                # i: index in population
                # j: index in source Population
                for i, j in [(0, 0), (2, 1), (3, 2)]:
                    for p in range(pop.ploidy()):
                        self.assertEqual(ind.allele(i, p), ind1.allele(j, p))
                for k in newpos:
                    self.assertEqual(ind.allele(k), 0)
        self.assertRaises(ValueError, pop.addLoci, [2], [8], ['l7'])
        #

    def testDeepcopy(self):
        'Testing deepcopy of population'
        pop = self.getPop(False, ancGen=3)
        initSex(pop)
        initGenotype(pop, freq=[0.2, 0.8])
        # shallow copy
        pop1 = pop
        initSex(pop1)
        initGenotype(pop1, freq= [0.8, 0.2])
        self.assertEqual(pop, pop1)
        # deep copy
        pop1 = pop.clone()
        self.assertEqual(pop, pop1)
        initSex(pop1)
        initGenotype(pop1, freq= [0.5, 0.5])
        self.assertNotEqual(pop, pop1)
        # using Python copy.copy
        pop1 = copy.copy(pop)
        self.assertEqual(pop, pop1)
        initSex(pop1)
        initGenotype(pop1, freq= [0.5, 0.5])
        self.assertEqual(pop, pop1)
        # using Python copy.deepcopy
        pop1 = copy.deepcopy(pop)
        self.assertEqual(pop, pop1)
        initSex(pop1)
        initGenotype(pop1, freq= [0.5, 0.5])
        self.assertNotEqual(pop, pop1)


    def testMergeSubPops(self):
        'Testing Population::MergeSubPops(subpops=[])'
        pop = self.getPop(size=[100, 20, 30, 80, 50, 60], subPopNames=['A', 'B', 'C', 'D', 'E', 'F'])
        pop1 = pop.clone()
        pop.mergeSubPops([1, 2, 4])
        self.assertEqual(pop.subPopSize(1), pop1.subPopSize(1)+pop1.subPopSize(2)+pop1.subPopSize(4))
        for (oldsp, newsp) in [(0, 0), (3, 2), (5, 3)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop.subPopSize(newsp))
            self.assertEqual(pop1.subPopName(oldsp), pop.subPopName(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop.individual(idx, newsp))
        # set new name to merged subpopulation
        pop = self.getPop(size=[100, 20, 30, 80, 50, 60], subPopNames=['A', 'B', 'C', 'D', 'E', 'F'])
        sp = pop.mergeSubPops([2, 1, 4], name='new')
        self.assertEqual(sp, 1)
        self.assertEqual(pop.subPopName(sp), 'new')
        self.assertEqual(pop.subPopNames(), ('A', 'new', 'D', 'F'))
        self.assertEqual(pop.subPopSize(1), pop1.subPopSize(1)+pop1.subPopSize(2)+pop1.subPopSize(4))
        #
        # merge to a specified subpopulation
        pop = self.getPop(size=[100, 20, 30, 80, 50, 60], subPopNames=['A', 'B', 'C', 'D', 'E', 'F'])
        pop1 = pop.clone()
        self.assertRaises(ValueError, pop.mergeSubPops, [1, 3, 4], toSubPop=2)
        pop.mergeSubPops([1, 3, 4], toSubPop=3)
        self.assertEqual(pop.subPopSize(2), pop1.subPopSize(1)+pop1.subPopSize(3)+pop1.subPopSize(4))
        for (oldsp, newsp) in [(0, 0), (2, 1), (5, 3)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop.subPopSize(newsp))
            self.assertEqual(pop1.subPopName(oldsp), pop.subPopName(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop.individual(idx, newsp))
        # set new name to merged subpopulation
        pop = self.getPop(size=[100, 20, 30, 80, 50, 60], subPopNames=['A', 'B', 'C', 'D', 'E', 'F'])
        pop1 = pop.clone()
        sp = pop.mergeSubPops([2, 1, 4], toSubPop=4, name='new')
        for (oldsp, newsp) in [(0, 0), (3, 1), (5, 3)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop.subPopSize(newsp))
            self.assertEqual(pop1.subPopName(oldsp), pop.subPopName(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop.individual(idx, newsp))
        self.assertEqual(sp, 2)
        self.assertEqual(pop.subPopName(sp), 'new')
        self.assertEqual(pop.subPopNames(), ('A', 'D', 'new', 'F'))
        self.assertEqual(pop.subPopSize(2), pop1.subPopSize(1)+pop1.subPopSize(2)+pop1.subPopSize(4))

    def testRemoveSubPops(self):
        'Testing Population::removeSubPops()'
        pop = self.getPop(size=[0, 100, 0, 20, 30, 0, 50], subPopNames=['A', 'B', 'C', 'D', 'E', 'F', 'G'])
        initSex(pop)
        initGenotype(pop, freq=[0.5, 0.5])
        pop1 = pop.clone()
        self.assertEqual(pop.numSubPop(), 7)
        pop.removeSubPops([x for x in range(7) if pop.subPopSize(x) == 0])
        self.assertEqual(pop.numSubPop(), 4)
        self.assertEqual(pop.subPopSizes(), (100, 20, 30, 50))
        for (oldsp, newsp) in [(1, 0), (3, 1), (4, 2), (6, 3)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop.subPopSize(newsp))
            self.assertEqual(pop1.subPopName(oldsp), pop.subPopName(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop.individual(idx, newsp))
        # remove subpop
        pop2 = pop.clone()
        pop.removeSubPops([1, 2])
        self.assertEqual(pop.subPopSizes(), (100, 50))
        for (oldsp, newsp) in [(0, 0), (3, 1)]:  # map of old and new id.
            self.assertEqual(pop2.subPopSize(oldsp), pop.subPopSize(newsp))
            self.assertEqual(pop2.subPopName(oldsp), pop.subPopName(newsp))
            for idx in range(pop2.subPopSize(oldsp)):
                self.assertEqual(pop2.individual(idx, oldsp), pop.individual(idx, newsp))
        self.assertRaises(IndexError, pop.removeSubPops, [8])
        # accept single input
        pop.removeSubPops(0)
        # 
        # now  for virtual subpopulation
        pop = self.getPop(size=[0, 100, 0, 20], subPopNames=['A', 'B', 'C', 'D'])
        initGenotype(pop, freq=[0.5, 0.5])
        initSex(pop)
        pop.setVirtualSplitter(SexSplitter())
        numFemale = pop.subPopSize([1,1])
        pop.removeSubPops([(1,0), 2])
        self.assertEqual(pop.numSubPop(), 3)
        self.assertEqual(pop.subPopSizes(), (0, numFemale, 20))
        for ind in pop.individuals(1):
            self.assertEqual(ind.sex(), FEMALE)
        # continue...
        pop.removeSubPops([(1,1), 2])
        self.assertEqual(pop.numSubPop(), 2)
        self.assertEqual(pop.subPopSizes(), (0, 0))
        #
        # test if allele lineage is correctly removed
        if moduleInfo()['alleleType'] != 'lineage':
            return
        pop = Population([10]*10, loci=10)
        pop.lineage()[:] = range(2000)
        pop.removeSubPops(range(1, 10, 2))
        # 0, 1, .... 20 | ... | 199
        # removed 200 ... 399
        # 400 ... 599
        for sp in range(5):
            for idx, ind in enumerate(pop.individuals(sp)):
                self.assertEqual(ind.lineage(), range(200*2*sp + idx*20, 200*2*sp + idx*20 + 20))
        
    def testRemoveIndividuals(self):
        'Testing Population::removeIndividuals(inds)'
        pop = self.getPop(size =[20, 100, 30], subPopNames=['sp1', 'sp2', 'sp3'])
        pop1 = pop.clone()
        pop.removeIndividuals([15])
        self.assertEqual(pop.subPopSizes(), (19, 100, 30))
        for idx in range(15):
            self.assertEqual(pop1.individual(idx), pop.individual(idx))
        for idx in range(15, pop.popSize()):
            self.assertEqual(pop1.individual(idx+1), pop.individual(idx))
        # accept single input
        pop.removeIndividuals(2)
        # 1) pop.removeIndividuals([500]) should yield an exception.
        pop = pop1.clone()
        self.assertRaises(IndexError, pop.removeIndividuals, 500)
        # 2) pop.removeIndividuals([]) should not change anything (self.assertEqual(pop, pop1))
        pop = pop1.clone()
        pop.removeIndividuals([])
        self.assertEqual(pop, pop1)
        # 3) pop.removeIndividuals(range(15, 25)) ...
        pop = pop1.clone()
        inds = list(range(15, 25))
        random.shuffle(inds)
        pop.removeIndividuals(inds)
        self.assertEqual(pop.subPopSizes(), (15, 95, 30))
        for idx in range(15):
            self.assertEqual(pop1.individual(idx), pop.individual(idx))
        for idx in range(24, pop.popSize()):
            self.assertEqual(pop1.individual(idx+10), pop.individual(idx))
        # 4) pop.removeIndividuals(range(15, 125)) removes the middle subpopulation
        #    and some individuals in subpopulation 0? Check if subpopulation name is handled correctly.
        pop = pop1.clone()
        inds = list(range(15, 125))
        random.shuffle(inds)
        pop.removeIndividuals(inds)
        self.assertEqual(pop.subPopSizes(), (15, 0, 25))
        for idx in range(15):
            self.assertEqual(pop1.individual(idx), pop.individual(idx))
        for idx in range(15, pop.popSize()):
            self.assertEqual(pop1.individual(idx+110), pop.individual(idx))
        self.assertEqual(pop.subPopNames(), pop1.subPopNames())
        # 5) pop.removeIndividuals(range(pop.subPopBegin(1), pop.subPopEnd(1))) removes the middle subpopulation.
        #    Check if subpopulation name is handled correctly.
        pop = pop1.clone()
        inds = list(range(pop.subPopBegin(1), pop.subPopEnd(1)))
        random.shuffle(inds)
        pop.removeIndividuals(inds)
        self.assertEqual(pop.subPopSizes(), (20, 0, 30))
        for idx in range(20):
            self.assertEqual(pop1.individual(idx), pop.individual(idx))
        for idx in range(21, pop.popSize()):
            self.assertEqual(pop1.individual(idx+100), pop.individual(idx))
        self.assertEqual(pop.subPopNames(), pop1.subPopNames())
        # 6) pop.removeIndividuals(range(pop.popSize())) removes all individuals in this population.
        pop = pop1.clone()
        inds = list(range(0, 150))
        random.shuffle(inds)
        pop.removeIndividuals(inds)
        self.assertEqual(pop.subPopSizes(), (0, 0, 0))
        self.assertEqual(pop.subPopNames(), pop1.subPopNames())
        #
        #
        # by ID?
        pop = self.getPop(size=[100, 200], loci=[2, 3, 1], ancGen=5,
            infoFields=['ind_id'])
        for gen in range(6):
            pop.useAncestralGen(gen)
            initGenotype(pop, freq=[0.5, 0.5])
        pop.useAncestralGen(0)
        IdTagger().reset(1)
        tagID(pop)
        exclude = set([random.randint(1, 1800) for x in range(600)])
        pop1 = pop.clone()
        pop1.removeIndividuals(IDs=list(exclude))
        sz = []
        sz1 = []
        for gen in range(6):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            sz.append(pop.popSize())
            sz1.append(pop1.popSize())
            id = set(pop.indInfo('ind_id'))
            id1 = set(pop1.indInfo('ind_id'))
            for e in exclude:
                self.assertEqual(e in id1, False)
            for ind in pop1.individuals():
                self.assertEqual(ind, pop.indByID(ind.ind_id), gen)
        self.assertEqual(sum(sz), sum(sz1) + len(exclude))
        # remove multiple individual
        pop = Population(10, infoFields='x')
        pop.setIndInfo([1, 2, 2, 3, 4, 5, 2, 3, 4, 3], 'x')
        pop.removeIndividuals(IDs=2, idField='x')
        self.assertEqual(pop.popSize(), 7)
        self.assertEqual(pop.indInfo('x'), (1, 3, 4, 5, 3, 4, 3))
        pop.removeIndividuals(IDs=[2,3,4], idField='x')
        self.assertEqual(pop.popSize(), 2)
        self.assertEqual(pop.indInfo('x'), (1, 5))
        # by filter function
        pop = Population(10, infoFields='x')
        pop.setIndInfo([1, 2, 2, 3, 4, 5, 2, 3, 4, 3], 'x')
        pop.removeIndividuals(filter=lambda ind: ind.x in [3, 4])
        self.assertEqual(pop.popSize(), 5)
        self.assertEqual(pop.indInfo('x'), (1, 2, 2, 5, 2))
        # test if allele lineage is correctly removed
        if moduleInfo()['alleleType'] != 'lineage':
            return
        pop = Population([10]*10, loci=10)
        pop.lineage()[:] = range(2000)
        pop.removeIndividuals(range(0, 100,2))
        lin = list(pop.lineage())
        self.assertEqual(len(lin), 1000)
        # the lineage of the remaining individuals are correctly handled
        for idx,ind in enumerate(pop.individuals()):
            self.assertEqual(ind.lineage(), range((idx*2+1)*20, (idx*2+2)*20))

    def testExtractSubPops(self):
        'Testing Population::extractSubPops()'
        pop = self.getPop(size=[0, 100, 0, 20, 30, 0, 50], subPopNames=['A', 'B', 'C', 'D', 'E', 'F', 'G'])
        initSex(pop)
        initGenotype(pop, freq=[0.5, 0.5])
        initLineage(pop, range(10))
        self.assertEqual(pop.numSubPop(), 7)
        pop1 = pop.extractSubPops([x for x in range(7) if pop.subPopSize(x) != 0])
        self.assertEqual(pop1.numSubPop(), 4)
        self.assertEqual(pop1.subPopSizes(), (100, 20, 30, 50))
        for (oldsp, newsp) in [(1, 0), (3, 1), (4, 2), (6, 3)]:  # map of old and new id.
            self.assertEqual(pop.subPopSize(oldsp), pop1.subPopSize(newsp))
            self.assertEqual(pop.subPopName(oldsp), pop1.subPopName(newsp))
            for idx in range(pop.subPopSize(oldsp)):
                self.assertEqual(pop.individual(idx, oldsp), pop1.individual(idx, newsp))
        # extract subpop
        pop2 = pop1.extractSubPops([1, 2])
        self.assertEqual(pop2.subPopSizes(), (20, 30))
        for (oldsp, newsp) in [(1, 0), (2, 1)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop2.subPopSize(newsp))
            self.assertEqual(pop1.subPopName(oldsp), pop2.subPopName(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop2.individual(idx, newsp))
        self.assertRaises(IndexError, pop.extractSubPops, [8])
        # accept single input
        pop.extractSubPops(0)
        # 
        # now  for virtual subpopulation
        pop = self.getPop(size=[0, 100, 0, 20], subPopNames=['A', 'B', 'C', 'D'])
        initGenotype(pop, freq=[0.5, 0.5])
        initSex(pop)
        pop.setVirtualSplitter(SexSplitter())
        numMale = pop.subPopSize([1,0])
        pop = pop.extractSubPops([(1,0), 3])
        self.assertEqual(pop.numSubPop(), 2)
        self.assertEqual(pop.subPopSizes(), (numMale, 20))
        for ind in pop.individuals(0):
            self.assertEqual(ind.sex(), MALE)
        # continue...
        pop1 = pop.extractSubPops([(0,1), 1])
        self.assertEqual(pop1.numSubPop(), 2)
        self.assertEqual(pop1.subPopSizes(), (0, 20))
        # remove multiple individual
        pop = Population(10, infoFields='x')
        pop.setIndInfo([1, 2, 2, 3, 4, 5, 2, 3, 4, 3], 'x')
        pop1 = pop.extractIndividuals(IDs=2, idField='x')
        self.assertEqual(pop1.popSize(), 3)
        self.assertEqual(pop1.indInfo('x'), (2, 2, 2))
        pop1 = pop.extractIndividuals(IDs=[2,3,4], idField='x')
        self.assertEqual(pop1.popSize(), 8)
        self.assertEqual(pop1.indInfo('x'), (2, 2, 3, 4, 2, 3, 4, 3))
        # by filter function
        pop = Population(10, infoFields='x')
        pop.setIndInfo([1, 2, 2, 3, 4, 5, 2, 3, 4, 3], 'x')
        pop1 = pop.extractIndividuals(filter=lambda ind: ind.x in [3, 4])
        self.assertEqual(pop1.popSize(), 5)
        self.assertEqual(pop1.indInfo('x'), (3, 4, 3, 4, 4))

    def testExtractSubPops(self):
        'Testing Population::extractSubPops()'
        pop = self.getPop(size=[0, 100, 0, 20, 30, 0, 50], subPopNames=['A', 'B', 'C', 'D', 'E', 'F', 'G'])
        initSex(pop)
        initGenotype(pop, freq=[0.5, 0.5])
        self.assertEqual(pop.numSubPop(), 7)
        pop1 = pop.extractSubPops([x for x in range(7) if pop.subPopSize(x) != 0])
        self.assertEqual(pop1.numSubPop(), 4)
    
    def testRearrangedExtractSubPops(self):
        'Testing Population::extractSubPops(subPops, true)'
        pop = self.getPop(size=[0, 100, 0, 20, 30, 0, 50], subPopNames=['A', 'B', 'C', 'D', 'E', 'F', 'G'])
        initSex(pop)
        initGenotype(pop, freq=[0.5, 0.5])
        self.assertEqual(pop.numSubPop(), 7)
        pop1 = pop.extractSubPops([x for x in range(6, 0, -1) if pop.subPopSize(x) != 0], True)
        self.assertEqual(pop1.numSubPop(), 4)
        self.assertEqual(pop1.subPopSizes(), (50, 30, 20, 100))
        for (oldsp, newsp) in [(6, 0), (4, 1), (3, 2), (1, 3)]:  # map of old and new id.
            self.assertEqual(pop.subPopSize(oldsp), pop1.subPopSize(newsp))
            self.assertEqual(pop.subPopName(oldsp), pop1.subPopName(newsp))
            for idx in range(pop.subPopSize(oldsp)):
                self.assertEqual(pop.individual(idx, oldsp), pop1.individual(idx, newsp))
        # extract subpop
        pop2 = pop1.extractSubPops([2, 1], True)
        self.assertEqual(pop2.subPopSizes(), (20, 30))
        for (oldsp, newsp) in [(2, 0), (1, 1)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop2.subPopSize(newsp))
            self.assertEqual(pop1.subPopName(oldsp), pop2.subPopName(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop2.individual(idx, newsp))
        self.assertRaises(IndexError, pop.extractSubPops, [8])
        # accept single input
        pop.extractSubPops(0, True)
        # 
        # now for virtual subpopulation
        pop = self.getPop(size=[0, 100, 0, 20], subPopNames=['A', 'B', 'C', 'D'])
        initGenotype(pop, freq=[0.5, 0.5])
        initSex(pop)
        pop.setVirtualSplitter(SexSplitter())
        numMale = pop.subPopSize([1,0])
        pop = pop.extractSubPops([3, (1,0), (1,1)], True)
        self.assertEqual(pop.numSubPop(), 3)
        self.assertEqual(pop.subPopSizes(), (20, numMale, 100-numMale))
        for ind in pop.individuals(1):
            self.assertEqual(ind.sex(), MALE)
        for ind in pop.individuals(2):
            self.assertEqual(ind.sex(), FEMALE)
        # continue...
        pop1 = pop.extractSubPops([(0,1), 1], True)
        self.assertEqual(pop1.numSubPop(), 2)
        self.assertEqual(pop1.subPopSize(1), numMale)


    def testExtractIndividuals(self):
        'Testing Population::removeIndividuals(inds)'
        pop = self.getPop(size =[20, 100, 30], subPopNames=['sp1', 'sp2', 'sp3'])
        initSex(pop)
        initGenotype(pop, freq=[0.4, 0.6])
        pop1 = pop.extractIndividuals()
        self.assertEqual(pop1.subPopSizes(), (0, 0, 0))
        self.assertEqual(pop1.subPopNames(), ('sp1', 'sp2', 'sp3'))
        pop1 = pop.extractIndividuals([15, 110, 120, 121])
        self.assertEqual(pop1.subPopSizes(), (1, 1, 2))
        for idx,oldidx in enumerate([15, 110, 120, 121]):
            self.assertEqual(pop1.individual(idx), pop.individual(oldidx))
        # accept single input
        pop.extractIndividuals(2)
        # 1) pop.extractIndividuals([500]) should yield an exception.
        self.assertRaises(IndexError, pop.extractIndividuals, 500)
        #
        # FIXME: Needs more tests
        #
        # by ID?
        pop = self.getPop(size=[100, 200], loci=[2, 3, 1], ancGen=5,
            infoFields=['ind_id'])
        for gen in range(6):
            pop.useAncestralGen(gen)
            initGenotype(pop, freq=[0.5, 0.5])
        pop.useAncestralGen(0)
        tagID(pop, reset=True)
        include = set([random.randint(1, 1800) for x in range(600)])
        pop1 = pop.extractIndividuals(IDs=list(include))
        sz1 = []
        for gen in range(6):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            sz1.append(pop1.popSize())
            id1 = set(pop1.indInfo('ind_id'))
            for e in id1:
                self.assertEqual(e in include, True)
            for ind in pop1.individuals():
                self.assertEqual(ind, pop.indByID(ind.ind_id), gen)
        self.assertEqual(sum(sz1), len(include))

    def testRemoveLoci(self):
        'Testing Population::removeLoci(loci=[], keep=[])'
        # Fixme: test loci, and keep, and test unordered parameters
        pop = self.getPop(size=[1, 2], loci=[2, 3, 1], ancGen=5)
        pop1 = pop.clone()
        # FIXME: test remove multiple loci from multiple chromosomes,
        # which may not be in order
        pop.removeLoci(2)
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
        # test the possibility of using loci names to remove loci
        pop = self.getPop(size=[1, 2], loci=[2, 3, 2], lociNames=['a%d' % x for x in range(7)], ancGen=5)
        pop1 = pop.clone()
        # FIXME: test remove multiple loci from multiple chromosomes,
        # which may not be in order
        pop.removeLoci('a2')
        pop.removeLoci(['a4', 'a5'])
        self.assertEqual(pop.numLoci(), (2,1,1))
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            for idx in range(pop.popSize()):
                ind = pop.individual(idx)
                ind1 = pop1.individual(idx)
                for loc in range(2):
                    self.assertEqual(ind.allele(loc), ind1.allele(loc))
                self.assertEqual(ind.allele(2), ind1.allele(3))
                self.assertEqual(ind.allele(3), ind1.allele(6))
        #
        #  testing remove the last locus
        if moduleInfo()['alleleType'] == 'binary':
            pop = Population(100, loci=10)
            initGenotype(pop, haplotypes=range(10))
            pop.removeLoci(9)
            self.assertEqual(pop.individual(0).genotype(), [0]+[1]*8 + [0] + [1]*8)
            #
            pop.removeLoci(2)
            self.assertEqual(pop.individual(0).genotype(), [0]+[1]*7 + [0] + [1]*7)
            #
            pop.addLoci(0, 2.8, 'test')
            self.assertEqual(pop.individual(0).genotype(), [0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1])
            pop.removeLoci('test')
            self.assertEqual(pop.individual(0).genotype(), [0] + [1]*7 + [0] + [1]*7)
        else:
            pop = Population(100, loci=10)
            initGenotype(pop, haplotypes=range(10))
            pop.removeLoci(9)
            self.assertEqual(pop.individual(0).genotype(), list(range(9)) + list(range(9)))
            #
            pop.removeLoci(2)
            self.assertEqual(pop.individual(0).genotype(), [0, 1, 3, 4, 5, 6, 7, 8, 0, 1, 3, 4, 5, 6, 7, 8])
            #
            pop.addLoci(0, 2.8, 'test')
            self.assertEqual(pop.individual(0).genotype(), [0, 1, 0, 3, 4, 5, 6, 7, 8, 0, 1, 0, 3, 4, 5, 6, 7, 8])
            pop.removeLoci('test')
            self.assertEqual(pop.individual(0).genotype(), [0, 1, 3, 4, 5, 6, 7, 8, 0, 1, 3, 4, 5, 6, 7, 8])


    def testRecodeAlleles(self):
        'Testing Population::recodeAlleles(alleles, loci)'
        pop = self.getPop(size=[10, 20], loci=[4, 5], ancGen=0, lociNames=['a%d' % x for x in range(9)])
        initSex(pop)
        initGenotype(pop, freq=[.2, .8])
        old = list(pop.genotype())
        # switch 0 and 1
        pop.recodeAlleles([1, 0])
        new = list(pop.genotype())
        for x,y in zip(old, new):
            self.assertEqual(x + y, 1)
        # clear to 0
        pop.recodeAlleles([0, 0], alleleNames=[['A']])
        self.assertEqual(pop.genotype(), [0]*(pop.totNumLoci()*pop.popSize()*pop.ploidy()))
        # use a function?
        def func(allele, locus):
            return allele + locus
        #
        pop.recodeAlleles(func, loci=1)
        self.assertEqual(pop.genotype(), [0, 1, 0, 0, 0, 0, 0, 0, 0] * (pop.popSize()*pop.ploidy()))
        # recode specified loci.
        pop = self.getPop(size=[10, 20], loci=[4, 5], ancGen=0, lociNames=['a%d' % x for x in range(9)])
        initSex(pop)
        initGenotype(pop, freq=[0, .2, .8])
        pop.recodeAlleles([0]*5, loci=(2, 4))
        for ind in pop.individuals():
            for loc in range(9):
                for p in range(2):
                    if loc in [2, 4]:
                        self.assertEqual(ind.allele(loc, p), 0)
                    else:
                        self.assertNotEqual(ind.allele(loc, p), 0)
        # using loci names
        pop.recodeAlleles([1]*10, loci=['a2', 'a4'])
        for ind in pop.individuals():
            for loc in ['a2', 'a4']:
                for p in range(2):
                    self.assertEqual(ind.allele(pop.locusByName(loc), p), 1)
        # FIXME: recode ancestral generations.

    def testResize(self):
        'Testing Population::resize(newSubPopSizes, propagate=false)'
        pop = self.getPop(size=[100, 20, 30], loci=[4, 5, 1])
        initSex(pop)
        initGenotype(pop, freq=[.2, .3, .5])
        pop1 = pop.clone()
        pop2 = pop.clone()
        # number of subpop mismatch
        self.assertRaises(ValueError, pop1.resize, [50, 50])
        # resize without propagation
        pop1.resize([50, 50, 80], propagate=False)
        for sp in range(pop1.numSubPop()):
            for i in range(min(pop1.subPopSize(sp), pop.subPopSize(sp))):
                self.assertEqual(pop1.individual(i, sp), pop.individual(i, sp))
            for i in range(min(pop1.subPopSize(sp), pop.subPopSize(sp)), pop1.subPopSize(sp)):
                self.assertEqual(pop1.individual(i, sp).genotype(), [0]*20)
        # resize with Population
        pop2.resize([50, 50, 80], propagate=True)
        for sp in range(pop1.numSubPop()):
            for i in range(pop2.subPopSize(sp)):
                self.assertEqual(pop2.individual(i, sp), pop.individual(i%pop.subPopSize(sp), sp))
        # resize from empty subpopulation?
        pop = self.getPop(size=[100, 0, 30, 0], loci=[4, 5, 1])
        self.assertEqual(pop.subPopSizes(), (100, 0, 30, 0))
        pop.resize([100, 20, 50, 20])
        self.assertEqual(pop.subPopSizes(), (100, 20, 50, 20))

    def testSplitSubPop(self):
        'Testing Population::splitSubPop(subPop, sizes)'
        pop = Population(size=[100, 80, 50], subPopNames=['A', 'B', 'C'])
        pop1 = pop.clone()
        self.assertRaises(ValueError, pop.splitSubPop, 1, [20, 70])
        ids = pop.splitSubPop(1, [20, 60])
        self.assertEqual(ids, (1, 2))
        self.assertEqual(pop1.subPopSize(1), pop.subPopSize(1)+pop.subPopSize(2))
        self.assertEqual(pop1.subPopName(1), pop.subPopName(1))
        self.assertEqual(pop1.subPopName(1), pop.subPopName(2))
        for idx in range(20):
            self.assertEqual(pop1.individual(idx, 1), pop. individual(idx, 1))
        for idx in range(20, 80):
            self.assertEqual(pop1.individual(idx, 1), pop. individual(idx-20, 2))
        for (oldsp, newsp) in [(0, 0), (2, 3)]:  # map of old and new id.
            self.assertEqual(pop1.subPopSize(oldsp), pop.subPopSize(newsp))
            self.assertEqual(pop1.subPopName(oldsp), pop.subPopName(newsp))
            for idx in range(pop1.subPopSize(oldsp)):
                self.assertEqual(pop1.individual(idx, oldsp), pop.individual(idx, newsp))
        # assign new names to split subpopulation
        pop = Population(size=[100, 80, 50])
        self.assertRaises(ValueError, pop.splitSubPop, 1, [20, 70], names=['A1'])
        ids = pop.splitSubPop(1, [20, 60], names=['A1', 'A2'])
        self.assertEqual(ids, (1, 2))
        self.assertEqual(pop.subPopName(1), 'A1')
        self.assertEqual(pop.subPopName(2), 'A2')
        self.assertEqual(pop.subPopNames(), ('', 'A1', 'A2', ''))
        # split by proportion
        pop = Population(size=[100, 80, 50])
        self.assertRaises(ValueError, pop.splitSubPop, 1, [0.2, 0.7])
        pop.splitSubPop(1, [0.3, 0.7])
        self.assertEqual(pop.subPopSizes(), (100, 24, 56, 50))
        pop.splitSubPop(0, [0.333, 0.667])
        self.assertEqual(pop.subPopSizes(), (33, 67, 24, 56, 50))

    def testSetSubPopByIndInfo(self):
        'Testing Population::setSubPopByIndInfo(field)'
        pop = self.getPop(subPopNames=['A', 'B'])
        for ind in pop.individuals():
            n = random.randint(-1, 5)
            ind.setInfo(n, 'x')
        pop1 = pop.setSubPopByIndInfo('x')
        self.assertEqual(pop.numSubPop(), 6)
        self.assertEqual(pop.subPopName(0), 'A')
        self.assertEqual(pop.subPopName(1), 'B')
        for i in range(2, 6):
            self.assertEqual(pop.subPopName(i), '')
        # apply this function to an empty information would crash simuPOP (issue #19)
        pop = Population(size=0, infoFields='a')
        pop.setSubPopByIndInfo('a')

    def testSortIndividuals(self):
        'Testing Population::sortIndividuals(infoFields)'
        pop = self.getPop(size=[1000, 2000], infoFields=['a', 'b'])
        initInfo(pop, lambda: random.randint(1, 5), infoFields=['a', 'b'])
        pop.sortIndividuals('a')
        for sp in range(2):
            for i in range(1, pop.subPopSize(sp)):
                self.assertTrue(pop.individual(i-1, sp).a <= pop.individual(i, sp).a)
        self.assertTrue(pop.individual(999).a > pop.individual(0, 1).a)
        # sorting in reverse order
        initInfo(pop, lambda: random.randint(1, 5), infoFields=['a', 'b'])
        pop.sortIndividuals('a', reverse=True)
        for sp in range(2):
            for i in range(1, pop.subPopSize(sp)):
                self.assertTrue(pop.individual(i-1, sp).a >= pop.individual(i, sp).a)
        self.assertTrue(pop.individual(999).a < pop.individual(0, 1).a)


            
    def testAddInfoFields(self):
        'Testing Population::addInfoFields(fields, init=0)'
        pop = self.getPop()
        pop1 = pop.addInfoFields('fitness', 5.0)
        self.assertEqual(pop.infoSize(), 2)
        pop1 = pop.addInfoFields('misc', 6.0)
        self.assertEqual(pop.infoSize(), 3)
        self.assertEqual(pop.indInfo('fitness'), tuple([5.0]*pop.popSize()))
        self.assertEqual(pop.indInfo('misc'), tuple([6.0]*pop.popSize()))
        pop = self.getPop()
        pop.addInfoFields(['x', 'fitness', 'misc'],  2.0)
        self.assertEqual(pop.infoSize(), 3)
        self.assertEqual(pop.indInfo('fitness'), tuple([2.0]*pop.popSize()))
        self.assertEqual(pop.indInfo('misc'), tuple([2.0]*pop.popSize()))
        # info field x is re-initialized
        self.assertEqual(pop.indInfo('x'), tuple([2.0]*pop.popSize()))
        # add again, but reinitialize ...
        pop.addInfoFields(['x', 'x', 'misc'], 5)
        self.assertEqual(pop.infoSize(), 3)
        self.assertEqual(pop.indInfo('misc'), tuple([5.0]*pop.popSize()))


    def testIndInfo(self):
        'Testing Population::indInfo(idx), indInfo(name), indInfo(idx, subPop)'
        'indInfo(name, subPop), setIndInfo(values, idx), setIndInfo(values, name)'
        'setIndInfo(values, idx, subPop), setIndInfo(values, name, subPop)'
        # no VSP, set and read info
        def testSetAndRead(pop):
            pop.setIndInfo([3], 'x')
            for idx, ind in enumerate(pop.individuals()):
                self.assertEqual(ind.info('x'), 3)
            self.assertEqual(pop.indInfo('x'), tuple([3]*pop.popSize()))
            self.assertEqual(pop.indInfo(0), tuple([3]*pop.popSize()))
            self.assertEqual(pop.indInfo('x', 0), tuple([3]*pop.subPopSize(0)))
            #
            pop.setIndInfo([1, 2], 'x', 0)
            pop.setIndInfo([3, 4], 0, 1)
            for idx, ind in enumerate(pop.individuals(0)):
                if idx % 2 == 0:
                    self.assertEqual(ind.info('x'), 1)
                else:
                    self.assertEqual(ind.info('x'), 2)
            self.assertEqual(pop.indInfo('x', 0), tuple([1, 2]*(int(pop.subPopSize(0)/2))))
            self.assertEqual(pop.indInfo(0, 1), tuple([3, 4]*(int(pop.subPopSize(1)/2))))
        #
        testSetAndRead(self.getPop())
        testSetAndRead(self.getPop(True))
        # test for virtual subpopulation
        def testVSPSetAndRead(pop):
            pop.setIndInfo([1, 2], 'x', [1, 0])
            pop.setIndInfo([3], 0, [1, 1])
            for idx, ind in enumerate(pop.individuals([1, 0])):
                self.assertEqual(ind.sex(), MALE)
                if idx % 2 == 0:
                    self.assertEqual(ind.info('x'), 1)
                else:
                    self.assertEqual(ind.info('x'), 2)
            for idx, ind in enumerate(pop.individuals([1, 1])):
                self.assertEqual(ind.sex(), FEMALE)
                self.assertEqual(ind.info('x'), 3)

            self.assertEqual(pop.indInfo('x', [1, 0]), tuple(([1, 2]*pop.subPopSize(1))[:pop.subPopSize([1, 0])]))
            self.assertEqual(pop.indInfo(0, [1, 1]), tuple([3]*pop.subPopSize([1, 1])))
        #
        self.assertRaises(ValueError, testVSPSetAndRead, self.getPop())
        testVSPSetAndRead(self.getPop(VSP=True))

    def testSetInfoFields(self):
        'Testing Population::setInfoFields(fields, init=0)'
        pop = self.getPop()
        pop1 = pop.setInfoFields(['fitness', 'misc'],  3)
        self.assertEqual(pop.infoSize(), 2)
        # info field x is removed
        self.assertEqual(pop.indInfo('fitness'), tuple([3]*pop.popSize()))
        self.assertEqual(pop.indInfo('misc'), tuple([3]*pop.popSize()))
        # allow set duplicated info fields
        pop1 = pop.setInfoFields(['fitness', 'fitness'],  1)
        self.assertEqual(pop.infoSize(), 1)

    def testClone(self):
        'Testing Population::clone()'
        pop = self.getPop(ancGen = 5)
        pop1 = pop.clone()
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            self.assertEqual(pop, pop1)

    def testSave(self):
        'Testing Population::save(filename)'
        pop = self.getPop(ancGen=5, infoFields=['a', 'b'])
        for gen in range(pop.ancestralGens(), -1, -1):
            initGenotype(pop, freq=[0.3, 0.7])
            initSex(pop)
            initInfo(pop, lambda:random.randint(0, 40), infoFields=['a', 'b'])
        pop.save("popout")
        pop1 = loadPopulation("popout")
        self.assertEqual(pop, pop1)
        self.assertEqual(pop.indInfo('a'), pop1.indInfo('a'))
        self.assertEqual(pop.indInfo('b'), pop1.indInfo('b'))
        #
        stat(pop, alleleFreq=list(range(pop.totNumLoci())))
        a = pop.dvars().alleleFreq[0][1]
        pop.save("popout")
        pop1 = loadPopulation("popout")
        self.assertEqual(a, pop1.dvars().alleleFreq[0][1])
        self.assertEqual(pop, pop1)
        #
        # testing the save of a population with non-pickleable objects
        pop.dvars().module_os = os
        self.assertTrue('module_os' in pop.vars())
        # module_os is not saved because it is a module
        pop.save('popout')
        pop1 = loadPopulation('popout')
        self.assertFalse('module_os' in pop1.vars())
        os.remove('popout')

    def testCrossPlatformLoad(self):
        'Testing loading populations created from other platform and allele types'
        localFile = 'sample_%d_%s_v3.pop' % ( \
            moduleInfo()['wordsize'],
            {'short': 'std', 
             'binary': 'ba',
             'long': 'la',
             'mutant': 'mu',
             'lineage': 'lin'
            }[moduleInfo()['alleleType']])
        if not os.path.isfile(localFile):
            print('Creating local pop file')
            pop = Population(10000, loci=100, infoFields=['a', 'ind_id'])
            initGenotype(pop, genotype=[0, 1, 1, 1, 0, 1, 1])
            initInfo(pop, values=[1, 2, 3, 4, 5], infoFields='a')
            stat(pop, alleleFreq=ALL_AVAIL)
            pop.save(localFile)
                
        for version in [0, 1, 2, 3]:
            for plat in [32, 64]:
                for mod in ['std', 'la', 'ba', 'mu', 'lin']:
                    if version == 0 and mod in ['lin', 'mu']:
                        continue
                    if version == 0:
                        popname = 'sample_%d_%s.pop' % (plat, mod)
                    else:
                        popname = 'sample_%d_%s_v1.pop' % (plat, mod)
                        if not os.path.isfile(popname):
                            print('Missing testing population name: %s' % popname)
                            continue
                    pop = Population()
                    #print('%s %s %s' % (version, plat, mod))
                    try:
                        pop = loadPopulation(popname)
                    except:
                        pass
                    self.assertEqual(pop.popSize(), 10000)
                    self.assertEqual(list(pop.indInfo('a')),
                        [1, 2, 3, 4, 5] * int(10000 / 5))
                    self.assertEqual(pop.genotype(),
                        ([0, 1, 1, 1, 0, 1, 1] * int(10000*100*2/7+1))[:10000*100*2])
                    if '_v3' in popname:
                        self.assertTrue(isinstance(pop.dvars().alleleFreq, defdict))

    def testVars(self):
        'Testing Population::vars(), vars(subPop), dvars(), dvars(subPop)'
        pop = self.getPop(size=1000, loci=[2, 4])
        initSex(pop)
        initGenotype(pop, freq=[.2, .3, .5])
        stat(pop, alleleFreq=list(range(0, 6)))
        pop1 = pop.clone()
        self.assertEqual(len(pop.vars()["alleleFreq"]), 6)
        self.assertEqual(len(pop.dvars().alleleFreq), 6)
        self.assertEqual(len(pop1.vars()["alleleFreq"]), 6)
        self.assertEqual(len(pop1.dvars().alleleFreq), 6)
        # with subPop
        pop = self.getPop(size=[20, 80], loci=[2, 4])
        initSex(pop)
        initGenotype(pop, freq=[.2, .3, .5])
        stat(pop, alleleFreq=list(range(0, 6)), vars='alleleFreq_sp')
        pop1 = pop.clone()
        self.assertEqual(len(pop.vars(0)["alleleFreq"]), 6)
        self.assertEqual(len(pop.dvars(1).alleleFreq), 6)
        self.assertEqual(len(pop1.vars(0)["alleleFreq"]), 6)
        self.assertEqual(len(pop1.dvars(1).alleleFreq), 6)

    def testSexSplitter(self):
        'Testing SexSplitter::SexSplitter()'
        pop = Population(size=[20, 80])
        initSex(pop)
        initGenotype(pop, freq=[0.4, 0.6])
        stat(pop, numOfMales=True, vars=['numOfMales_sp', 'numOfFemales_sp'])
        pop.setVirtualSplitter(SexSplitter())
        self.assertEqual(pop.subPopSize([1, 0]), pop.dvars(1).numOfMales)
        self.assertEqual(pop.subPopSize([1, 1]), pop.dvars(1).numOfFemales)
        self.assertEqual(pop.subPopName([1, 0]), 'Male')
        self.assertEqual(pop.subPopName([1, 1]), 'Female')
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.sex(), MALE)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.sex(), FEMALE)
        # test nested virtual subpopulation
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.sex(), MALE)
            for ind1 in pop.individuals([0, 1]):
                self.assertEqual(ind1.sex(), FEMALE)
        numMale = 0
        numFemale = 0
        for ind in pop.individuals(1):
            if ind.sex() == MALE:
                numMale += 1
            else:
                numFemale += 1
        #print numMale, numFemale
        self.assertEqual(numMale == 0, False)
        self.assertEqual(numFemale == 0, False)

    def testAffectionSplitter(self):
        'Testing AffectionSplitter::AffectionSplitter()'
        pop = Population(size=[20, 80], loci=[1, 2])
        initSex(pop)
        initGenotype(pop, freq=[0.4, 0.6])
        maPenetrance(pop, loci=0, wildtype=0, penetrance=[0.2, 0.4, 0.8])
        stat(pop, numOfAffected=True, vars=['numOfAffected_sp', 'numOfUnaffected_sp'])
        pop.setVirtualSplitter(AffectionSplitter())
        self.assertEqual(pop.subPopSize([1, 1]), pop.dvars(1).numOfAffected)
        self.assertEqual(pop.subPopSize([1, 0]), pop.dvars(1).numOfUnaffected)
        self.assertEqual(pop.subPopName([1, 0]), 'Unaffected')
        self.assertEqual(pop.subPopName([1, 1]), 'Affected')
        for ind in pop.individuals([1, 1]):
            self.assertEqual(ind.affected(), True)
        for ind in pop.individuals([1, 0]):
            self.assertEqual(ind.affected(), False)
        numAffected = 0
        numUnaffected = 0
        for ind in pop.individuals(1):
            if ind.affected():
                numAffected += 1
            else:
                numUnaffected += 1
        self.assertEqual(numAffected == 0, False)
        self.assertEqual(numUnaffected == 0, False)

    def testInfoSplitter(self):
        'Testing InfoSplitter::InfoSplitter(field, values=[], cutoff=[])'
        pop = Population(1000, infoFields=['x'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(InfoSplitter('x', values=list(range(10, 15))))
        self.assertEqual(pop.numVirtualSubPop(), 5)
        infos = list(pop.indInfo('x'))
        self.assertEqual(pop.subPopName([0, 0]), "x = 10")
        self.assertEqual(pop.subPopName([0, 1]), "x = 11")
        self.assertEqual(pop.subPopName([0, 4]), "x = 14")
        self.assertEqual(pop.subPopSize([0, 0]), infos.count(10))
        self.assertEqual(pop.subPopSize([0, 1]), infos.count(11))
        self.assertEqual(pop.subPopSize([0, 2]), infos.count(12))
        self.assertEqual(pop.subPopSize([0, 3]), infos.count(13))
        self.assertEqual(pop.subPopSize([0, 4]), infos.count(14))
        for i in range(5):
            for ind in pop.individuals([0, i]):
                self.assertEqual(ind.info('x'), 10+i)
        # test cutoff
        pop.setVirtualSplitter(InfoSplitter('x', cutoff=[11.5, 13.5]))
        self.assertEqual(pop.subPopName([0, 0]), "x < 11.5")
        self.assertEqual(pop.subPopName([0, 1]), "11.5 <= x < 13.5")
        self.assertEqual(pop.subPopName([0, 2]), "x >= 13.5")
        self.assertEqual(pop.subPopSize([0, 0]), infos.count(10) + infos.count(11))
        self.assertEqual(pop.subPopSize([0, 1]), infos.count(12) + infos.count(13))
        self.assertEqual(pop.subPopSize([0, 2]),
            sum([infos.count(x) for x in range(14, 21)]))
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.info('x') < 11.5, True)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(11.5 <= ind.info('x') < 13.5, True)
        for ind in pop.individuals([0, 2]):
            self.assertEqual(ind.info('x') >=13.5, True)
        # test range
        pop.setVirtualSplitter(InfoSplitter('x', ranges=[[11.5, 13.5], [9.5, 12.5]]))
        self.assertEqual(pop.numVirtualSubPop(), 2)
        self.assertEqual(pop.subPopName([0, 0]), "11.5 <= x < 13.5")
        self.assertEqual(pop.subPopName([0, 1]), "9.5 <= x < 12.5")
        self.assertEqual(pop.subPopSize([0, 0]), infos.count(12) + infos.count(13))
        self.assertEqual(pop.subPopSize([0, 1]), infos.count(10) + infos.count(11) + infos.count(12))
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.info('x') >= 11.5 and ind.info('x') < 13.5, True)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(9.5 <= ind.info('x') < 12.5, True)

    def testProportionSplitter(self):
        'Testing ProportionSplitter::ProportionSplitter(proportions=[])'
        pop = Population(10)
        pop.setVirtualSplitter(ProportionSplitter([0.01]*100))
        for i in range(100):
            self.assertEqual(pop.subPopName([0, i]), "Prop 0.01")
            if i != 99:
                self.assertEqual(pop.subPopSize([0, i]), 0)
            else:
                # the last vsp is specially treated to avoid such problem.
                self.assertEqual(pop.subPopSize([0, i]), 10)
        #
        pop = Population(1000)
        pop.setVirtualSplitter(ProportionSplitter([0.4, 0.6]))
        self.assertEqual(pop.subPopSize([0, 0]), 400)
        self.assertEqual(pop.subPopSize([0, 1]), 600)

    def testRangeSplitter(self):
        'Testing RangeSplitter::RangeSplitter(ranges)'
        pop = Population(100)
        pop.setVirtualSplitter(RangeSplitter(ranges=[[10, 20], [80, 200]]))
        self.assertEqual(pop.subPopName([0, 0]), "Range [10, 20)")
        self.assertEqual(pop.subPopName([0, 1]), "Range [80, 200)")
        self.assertEqual(pop.subPopSize([0, 0]), 10)
        self.assertEqual(pop.subPopSize([0, 1]), 20)

    def testGenotypeSplitter(self):
        'Testing GenotypeSplitter::GenotypeSplitter(loci(or locus), alleles, phase=False)'
        pop = Population(1000, loci=[2, 3])
        initSex(pop)
        initGenotype(pop, freq=[0.3, 0.7])
        pop.setVirtualSplitter(GenotypeSplitter(loci=1, alleles=[[0, 0], [1, 0]], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 0")
        self.assertEqual(pop.subPopName([0, 1]), "Genotype 1: 1 0")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0), 0)
            self.assertEqual(ind.allele(1, 1), 0)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0), 1)
            self.assertEqual(ind.allele(1, 1), 0)
        stat(pop, genoFreq=[1])
        self.assertEqual(pop.subPopSize([0, 0]), pop.dvars().genoNum[1][(0,0)])
        self.assertEqual(pop.subPopSize([0, 1]), pop.dvars().genoNum[1][(1,0)])
        # non-phased case
        pop.setVirtualSplitter(GenotypeSplitter(loci=1, alleles=[[0, 0], [1, 0]], phase=False))
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0), 0)
            self.assertEqual(ind.allele(1, 1), 0)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0)==0 or ind.allele(1, 0)==1, True)
            self.assertEqual(ind.allele(1, 1)==0 or ind.allele(1, 1)==1, True)
        # multiple loci
        pop.setVirtualSplitter(GenotypeSplitter(loci=[0, 1], alleles=[0, 1, 1, 1], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 0, 1: 0 1 1 1")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(0, 0), 0)
            self.assertEqual(ind.allele(0, 1), 1)
            self.assertEqual(ind.allele(1, 0), 1)
            self.assertEqual(ind.allele(1, 1), 1)
        # multiple genotype at the same locus
        pop.setVirtualSplitter(GenotypeSplitter(loci=1, alleles=[0, 1, 1, 1], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 1 1 1")
        for ind in pop.individuals([0, 0]):
             self.assertEqual(ind.allele(1, 0)==1  or ind.allele(1, 0)==0, True)
             self.assertEqual(ind.allele(1, 1), 1)
        # haploid case
        pop = Population(1000, ploidy = 1, loci=[2, 3])
        initSex(pop)
        initGenotype(pop, freq=[0.3, 0.7])
        pop.setVirtualSplitter(GenotypeSplitter(loci=1, alleles=[[0, 1], [2]], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 1")
        self.assertEqual(pop.subPopName([0, 1]), "Genotype 1: 2")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0)==1 or ind.allele(1, 0)==0, True)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0), 2)
        # use of names
        pop = Population(1000, ploidy = 1, loci=[2, 3], lociNames=['a%d' % x for x in range(5)])
        initSex(pop)
        initGenotype(pop, freq=[0.3, 0.7])
        pop.setVirtualSplitter(GenotypeSplitter(loci='a1', alleles=[[0, 1], [2]], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype a1: 0 1")
        self.assertEqual(pop.subPopName([0, 1]), "Genotype a1: 2")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0)==1 or ind.allele(1, 0)==0, True)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0), 2)


    def testCombinedSplitter(self):
        'Testing CombinedSplitter:: CombinedSplitter(splitters=[])'
        pop = Population(1000, loci=[2, 3])
        initSex(pop)
        initGenotype(pop, freq=[0.3, 0.7])
        pop.setVirtualSplitter(CombinedSplitter([
            GenotypeSplitter(loci=1, alleles=[[0, 0], [1, 0]], phase=True),
            SexSplitter()]))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 0")
        self.assertEqual(pop.subPopName([0, 1]), "Genotype 1: 1 0")
        self.assertEqual(pop.subPopName([0, 2]), "Male")
        self.assertEqual(pop.subPopName([0, 3]), "Female")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0), 0)
            self.assertEqual(ind.allele(1, 1), 0)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0), 1)
            self.assertEqual(ind.allele(1, 1), 0)
        for ind in pop.individuals([0, 2]):
            self.assertEqual(ind.sex(), MALE)
        for ind in pop.individuals([0, 3]):
            self.assertEqual(ind.sex(), FEMALE)
        stat(pop, numOfMales=True, vars='numOfFemales_sp')
        self.assertEqual(pop.subPopSize([0, 3]), pop.dvars(0).numOfFemales)
        #
        # combined splitter with vspMap
        #
        pop = Population(1000, loci=[2, 3])
        initSex(pop)
        initGenotype(pop, freq=[0.3, 0.7])
        pop.setVirtualSplitter(CombinedSplitter([
            GenotypeSplitter(loci=1, alleles=[[0, 0], [1, 0]], phase=True),
            SexSplitter()], vspMap=[[0,2], [1], [3]]))
        self.assertEqual(pop.numVirtualSubPop(), 3)
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 0 or Male")
        self.assertEqual(pop.subPopName([0, 1]), "Genotype 1: 1 0")
        self.assertEqual(pop.subPopName([0, 2]), "Female")
        for ind in pop.individuals([0, 0]):
            self.assertTrue((ind.allele(1, 0) == 0 and ind.allele(1, 1) == 0) or ind.sex() == MALE)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0), 1)
            self.assertEqual(ind.allele(1, 1), 0)
        for ind in pop.individuals([0, 2]):
            self.assertEqual(ind.sex(), FEMALE)
        #
        pop = Population(1000, loci=[2], infoFields='a')
        initInfo(pop, random.randint(0, 3), infoFields='a')
        pop.setVirtualSplitter(CombinedSplitter([InfoSplitter(field='a', values=list(range(4)))], vspMap=[[0,2], [1,3]]))
        self.assertEqual(pop.numVirtualSubPop(), 2)
        self.assertEqual(pop.subPopName([0, 0]), "a = 0 or a = 2")
        for ind in pop.individuals([0,0]):
            self.assertTrue(ind.info('a') in [0, 2])
        for ind in pop.individuals([0,1]):
            self.assertTrue(ind.info('a') in [1, 3])
        self.assertEqual(pop.subPopSize([0,0]) + pop.subPopSize([0,1]), pop.popSize())


    def testProductSplitter(self):
        'Testing CombinedSplitter::ProductSplitter(splitters=[])'
        pop = Population(1000, loci=[2, 3])
        initSex(pop)
        initGenotype(pop, freq=[0.3, 0.7])
        pop.setVirtualSplitter(ProductSplitter([
            GenotypeSplitter(loci=1, alleles=[[0, 0], [1, 0], [0, 1], [1, 1]], phase=True),
            SexSplitter()]))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 0, Male")
        self.assertEqual(pop.subPopName([0, 1]), "Genotype 1: 0 0, Female")
        self.assertEqual(pop.subPopName([0, 2]), "Genotype 1: 1 0, Male")
        self.assertEqual(pop.subPopName([0, 3]), "Genotype 1: 1 0, Female")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0), 0)
            self.assertEqual(ind.allele(1, 1), 0)
            self.assertEqual(ind.sex(), MALE)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0), 0)
            self.assertEqual(ind.allele(1, 1), 0)
            self.assertEqual(ind.sex(), FEMALE)
        for ind in pop.individuals([0, 2]):
            self.assertEqual(ind.allele(1, 0), 1)
            self.assertEqual(ind.allele(1, 1), 0)
            self.assertEqual(ind.sex(), MALE)
        for ind in pop.individuals([0, 3]):
            self.assertEqual(ind.allele(1, 0), 1)
            self.assertEqual(ind.allele(1, 1), 0)
            self.assertEqual(ind.sex(), FEMALE)
        stat(pop, numOfMales=True)
        for x in range(8):
            self.assertTrue(pop.subPopSize([0,x]) > 0)
        self.assertEqual(sum([pop.subPopSize([0,x]) for x in range(0, 8, 2)]), pop.dvars().numOfMales)
        self.assertEqual(sum([pop.subPopSize([0,x]) for x in range(1, 8, 2)]), pop.dvars().numOfFemales)

    def testIndByID(self):
        'Testing Population::indByID()'
        pop = self.getPop(size=[200]*4, ancGen=3, infoFields=['ind_id'])
        IdTagger().reset(1)
        tagID(pop)
        for i in range(400):
            id = random.randint(1, 800*4)
            ind = pop.indByID(id)
            self.assertEqual(ind.ind_id, id)
        self.assertRaises(IndexError, pop.indByID, 8000)
 
    def testIdentifyFamilies(self):
        'Testing Pedigree::identifyFamily'
        pop = Population(100, infoFields=['ind_id', 'ped_id'], ancGen=1)
        tagID(pop, reset=True)
        ped = Pedigree(pop, fatherField='', motherField='', infoFields=ALL_AVAIL)
        pedSize = ped.identifyFamilies(pedField='ped_id')
        self.assertEqual(pedSize, tuple([1]*100))
        self.assertEqual(ped.indInfo('ped_id'), tuple(range(100)))
        pop.evolve(
            matingScheme=RandomSelection(ops=[
                CloneGenoTransmitter(), IdTagger()]),
            gen = 1
        )
        ped = Pedigree(pop, fatherField='', motherField='', infoFields=ALL_AVAIL)
        pedSize = ped.identifyFamilies(pedField='ped_id')
        self.assertEqual(pedSize, tuple([1]*200))
        #
        pop = Population(100, infoFields=['ind_id', 'father_id', 'ped_id'], ancGen=1)
        tagID(pop, reset=True)
        pop.evolve(
            matingScheme=RandomSelection(ops=[
                CloneGenoTransmitter(), IdTagger(),
                    PedigreeTagger(infoFields='father_id')]),
            gen = 1
        )
        ped = Pedigree(pop, motherField='', infoFields=ALL_AVAIL)
        pedSize = ped.identifyFamilies(pedField='ped_id')
        self.assertEqual(sum(pedSize), 200)
        for idx, sz in enumerate(pedSize):
            if sz > 1:
                p = ped.extractIndividuals(IDs=idx, idField='ped_id')
                self.assertEqual(len(list(p.allIndividuals())), sz)
        #

    def testIdentifyAncestors(self):
        'Testing pedigree::identifyAncestors'
        pop = Population(100, infoFields=['ind_id', 'father_id'], ancGen=1)
        tagID(pop, reset=True)
        pop.evolve(
            matingScheme=RandomSelection(ops=[
                CloneGenoTransmitter(), IdTagger(),
                    PedigreeTagger(infoFields='father_id')]),
            gen = 1
        )
        ped = Pedigree(pop, motherField='', infoFields=ALL_AVAIL)
        IDs = ped.identifyAncestors()
        self.assertTrue(len(IDs) > 100)
        # two parents
        pop = Population(100, infoFields=['ind_id', 'father_id', 'mother_id'], ancGen=3)
        tagID(pop, reset=True)
        pop.evolve(
            initOps = InitSex(),
            matingScheme=RandomMating(ops=[
                MendelianGenoTransmitter(),
                IdTagger(),
                PedigreeTagger()]),
            gen = 5
        )
        pop.asPedigree()
        IDs = pop.identifyAncestors()
        self.assertTrue(len(IDs) > 200)
        # ancestors of selected parents
        IDs = pop.identifyAncestors(501)
        self.assertTrue(len(IDs) > 1 + 2 + 4)

    def testIdentifyOffspring(self):
        'Testing pedigree::offspring'
        pop = Population(100, infoFields=['ind_id', 'father_id'], ancGen=1)
        tagID(pop, reset=True)
        pop.evolve(
            matingScheme=RandomSelection(ops=[
                CloneGenoTransmitter(), IdTagger(),
                    PedigreeTagger(infoFields='father_id')]),
            gen = 1
        )
        ped = Pedigree(pop, motherField='', infoFields=ALL_AVAIL)
        IDs = ped.identifyOffspring(1)
        # two parents
        pop = Population(100, infoFields=['ind_id', 'father_id', 'mother_id'], ancGen=3)
        tagID(pop, reset=True)
        pop.evolve(
            initOps = InitSex(),
            matingScheme=RandomMating(ops=[
                MendelianGenoTransmitter(),
                IdTagger(),
                PedigreeTagger()]),
            gen = 5
        )
        pop.asPedigree()
        pop.useAncestralGen(3)
        anc = pop.indInfo('ind_id')[:10]
        IDs = pop.identifyOffspring(anc)
        len(IDs) > 20

    def testDescribeEvolProcess(self):
        'Testing population::evolve(dryrun=True'
        pop = Population(100, loci=3)
        tmp = sys.stdout 
        with open(os.devnull, 'w') as sys.stdout:
            pop.evolve(preOps=InitSex(),
                matingScheme=RandomMating(), dryrun=True)
        sys.stdout = tmp

    def testLineage(self):
        if moduleInfo()['alleleType'] != 'lineage':
            return
        IdTagger().reset(1)
        pop = Population(100, infoFields='ind_id', loci=10)
        tagID(pop)
        self.assertEqual(pop.lineage(), [0] * 2000)
        for ind in pop.allIndividuals():
            self.assertEqual(ind.lineage(), [0] * 20)
        # set lineage per loci
        initLineage(pop, range(20))
        self.assertEqual(pop.lineage(), list(range(20)) * 100)
        for ind in pop.allIndividuals():
            self.assertEqual(ind.lineage(), list(range(20)))
        # set lineage per ploidy
        initLineage(pop, list(range(2)), mode=PER_PLOIDY)
        self.assertEqual(pop.lineage(), (([0] * 10) + ([1] * 10)) * 100)
        initLineage(pop, list(range(200)), mode=PER_PLOIDY)
        self.assertEqual(pop.lineage(), sum([10 * [i] for i in range(200)], []))
        # set lineage per individual
        initLineage(pop, list(range(100)), mode=PER_INDIVIDUAL)
        self.assertEqual(pop.lineage(), sum([20 * [i] for i in range(100)], []))

       
    def testAllIndividuals(self):
        'Testing population::allIndividuals'
        pop = Population([100]*10, loci=3, ancGen=-1)
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        pop1 = pop.clone()
        pop1.mergeSubPops(list(range(5)))
        self.assertEqual(pop1.numSubPop(), 6)
        pop.push(pop1)
        #
        sex = [x.sex() for x in pop.allIndividuals()]
        self.assertEqual(len(sex), 2000)
        sex = [x.sex() for x in pop.allIndividuals(ancGens=1)]
        self.assertEqual(len(sex), 1000)
        # virtual subpopulation
        sex = [x.sex() for x in pop.allIndividuals(subPops=[(0,0)], ancGens=1)]
        self.assertEqual(len(sex), 50)
        def goThrough():
            [x.sex() for x in pop.allIndividuals(subPops=[(x,0) for x in range(10)], ancGens=0)]
        self.assertRaises(IndexError, goThrough)
        # this is OK
        sex = [x.sex() for x in pop.allIndividuals(subPops=[(x,0) for x in range(10)], ancGens=1)]
        for ind in pop.allIndividuals(subPops=[(x,0) for x in range(10)], ancGens=1):
            self.assertEqual(ind.sex(), MALE)
        self.assertEqual(len(sex), 500)
        # this is also OK
        sex = [x.sex() for x in pop.allIndividuals(subPops=[(ALL_AVAIL,0)])]
        self.assertEqual(len(sex), 1000)
        sex = [x.sex() for x in pop.allIndividuals(subPops=[(ALL_AVAIL,1)], ancGens=0)]
        for ind in pop.allIndividuals(subPops=[(ALL_AVAIL,1)]):
            self.assertEqual(ind.sex(), FEMALE)
        self.assertEqual(len(sex), 500)
        # this is also OK
        sex = [x.sex() for x in pop.allIndividuals(subPops=[(0, ALL_AVAIL)])]
        self.assertEqual(len(sex), 600)
        self.assertEqual(sex[:250], [MALE]*250)
        self.assertEqual(sex[250:500], [FEMALE]*250)
        self.assertEqual(sex[500:550], [MALE]*50)
        self.assertEqual(sex[550:600], [FEMALE]*50)
        # this is also OK
        sex = [x.sex() for x in pop.allIndividuals(subPops=[(ALL_AVAIL, ALL_AVAIL)])]
        self.assertEqual(len(sex), 2000)
        self.assertEqual(sex[:250], [MALE]*250)
        self.assertEqual(sex[250:500], [FEMALE]*250)
        self.assertEqual(sex[500:550], [MALE]*50)
        self.assertEqual(sex[550:600], [FEMALE]*50)

    def testIDOnlyPedigree(self):
        'Testing pedigrees with only individual ID'
        pop = Population(500, ancGen=1, infoFields='ind_id')
        initSex(pop)
        tagID(pop)
        pop.asPedigree(fatherField='', motherField='')
        pop.save('ind.ped')
        ped = loadPedigree('ind.ped')
        os.remove('ind.ped')


    def testSaveLoadPedigree(self):
        'Testing function loadPedigree'
        pop = Population(500, infoFields=['ind_id', 'father_id'], ancGen=-1, loci=[1])
        tagID(pop, reset=True)
        pop.evolve(
            initOps = InitGenotype(freq=[0.5, 0.5]),
            matingScheme=RandomSelection(ops=[
                CloneGenoTransmitter(), IdTagger(),
                    PedigreeTagger(infoFields='father_id', output='>>test.ped')]),
            gen = 10
        )
        #
        pop.asPedigree(motherField='')
        pop.save('test1.ped', loci=0)
        ped1 = loadPedigree('test1.ped')
        ped1.save('test2.ped', loci=0)
        ped2 = loadPedigree('test2.ped')
        self.assertEqual(ped1, ped2)
        #
        pop.save('test1.ped', loci=0, infoFields='ind_id')
        ped1 = loadPedigree('test1.ped', infoFields='ind_id1')
        ped1.save('test2.ped', loci=0, infoFields='ind_id')
        ped2 = loadPedigree('test2.ped', infoFields='ind_id1')
        self.assertEqual(ped1, ped2)
        #
        ped = loadPedigree('test.ped', motherField='')
        self.assertEqual(ped.ancestralGens(), 10)
        for gen in range(11):
            ped.useAncestralGen(gen)
            if gen == 10:
                self.assertTrue(ped.popSize() < 500)
            else:
                self.assertEqual(ped.popSize(), 500)
        self.assertEqual(ped.individual(0).father_id, 0)
        ped.useAncestralGen(0)
        self.assertNotEqual(ped.individual(0).father_id, 0)
        # two parents
        pop = Population(500, loci=[2], ancGen=-1, infoFields=['ind_id', 'father_id', 'mother_id'])
        tagID(pop, reset=True)
        pop.evolve(
            initOps = [
                InitSex(),
                InitGenotype(freq=[0.5, 0.5]),
            ],
            matingScheme=RandomMating(ops=[
                MendelianGenoTransmitter(),
                IdTagger(),
                PedigreeTagger(output='>>test.ped')]),
            gen = 20
        )
        #
        pop.asPedigree()
        pop.save('test1.ped', loci=0)
        ped1 = loadPedigree('test1.ped')
        ped1.save('test2.ped', loci=0)
        ped2 = loadPedigree('test2.ped')
        self.assertEqual(ped1, ped2)
        #
        pop.save('test1.ped', loci=0, infoFields='ind_id')
        ped1 = loadPedigree('test1.ped', infoFields='ind_id1')
        ped1.save('test2.ped', loci=0, infoFields='ind_id')
        ped2 = loadPedigree('test2.ped', infoFields='ind_id1')
        self.assertEqual(ped1, ped2)
        #
        ped = loadPedigree('test.ped')
        self.assertEqual(ped.ancestralGens(), 20)
        for gen in range(21):
            ped.useAncestralGen(gen)
            if gen == 20:
                self.assertTrue(ped.popSize() < 500)
            else:
                self.assertEqual(ped.popSize(), 500)
        self.assertEqual(ped.individual(0).father_id, 0)
        self.assertEqual(ped.individual(0).mother_id, 0)
        ped.useAncestralGen(0)
        self.assertNotEqual(ped.individual(0).father_id, 0)
        self.assertNotEqual(ped.individual(0).mother_id, 0)
        # cleanup
        for file in ['test.ped', 'test1.ped', 'test2.ped']:
            os.remove(file)

    def testDiscardIf(self):
        'Testing operator DiscardIf'
        pop = Population(1000, loci=2, infoFields=['a', 'b'])
        initInfo(pop, [1, 2, 3], infoFields='a')
        initInfo(pop, [1, 2, 3, 4], infoFields='b')
        discardIf(pop, 'a==b')
        for ind in pop.individuals():
            self.assertNotEqual(ind.a, ind.b)
        def func(a):
            return a <= 1
        discardIf(pop, func)
        for ind in pop.individuals():
            self.assertTrue(ind.a > 1)
        self.assertTrue(pop.popSize() < 1000)
        # testing virtual subpopulation
        pop.setVirtualSplitter(InfoSplitter(field='b', values=[3]))
        discardIf(pop, True, subPops=[(0,0)])
        for ind in pop.individuals():
            self.assertTrue(ind.b != 3)
        #
        # test discard by probability
        pop = Population(1000)
        discardIf(pop, cond='0.5')
        self.assertTrue(pop.popSize() > 450)
        self.assertTrue(pop.popSize() < 550)
        # probability from expression
        pop = Population(1000, loci=2, infoFields='a')
        initInfo(pop, [0.2], infoFields='a')
        discardIf(pop, cond='a+0.05')
        self.assertTrue(pop.popSize() > 700)
        self.assertTrue(pop.popSize() < 800)
        #
        pop = Population([1000, 1000], loci=2, infoFields='a')
        initInfo(pop, [0.2], infoFields='a', subPops=0)
        initInfo(pop, [0.8], infoFields='a', subPops=1)
        discardIf(pop, cond='a', subPops=0)
        self.assertTrue(pop.subPopSize(0) > 750)
        self.assertTrue(pop.subPopSize(0) < 850)
        self.assertEqual(pop.subPopSize(1), 1000)
        #
        discardIf(pop, cond='a', subPops=1)
        self.assertTrue(pop.subPopSize(1) > 150)
        self.assertTrue(pop.subPopSize(1) < 250)
        #
        #
        pop = Population(1000, loci=2, infoFields='a')
        initInfo(pop, [0.2], infoFields='a')
        discardIf(pop, cond=lambda a: a + 0.1)
        self.assertTrue(pop.popSize() > 650)
        self.assertTrue(pop.popSize() < 750)
        #
        pop = Population([1000, 1000], loci=2, infoFields='a')
        initInfo(pop, [0.2], infoFields='a', subPops=0)
        initInfo(pop, [0.8], infoFields='a', subPops=1)
        discardIf(pop, cond=lambda a: a, subPops=0)
        self.assertTrue(pop.subPopSize(0) > 750)
        self.assertTrue(pop.subPopSize(0) < 850)
        self.assertEqual(pop.subPopSize(1), 1000)
        #
        discardIf(pop, cond=lambda a: a-0.3, subPops=1)
        self.assertTrue(pop.subPopSize(1) > 450)
        self.assertTrue(pop.subPopSize(1) < 550)
        # error handling
        self.assertRaises(Exception, discardIf, pop, cond='a2')


    def testMutants(self):
        'Testing function Population.mutants'
        pop = Population([4,6], loci=20)
        pop.setGenotype([0,0,1])
        # mutants are at
        # 0 0 1=2 0 0 1=5 0 0 1 0 0 1 0 0 1 0 0 1=17 0 0 <- 6
        # 1=0 0 0 1=3 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1=18 0 <- 7
        # 0 1=1 ....................................1=19 <- 7
        # 0 0 1=2  ......................................<- 6
        mutants = list(pop.mutants())
        self.assertEqual(len(mutants), 3 * 40 + 13)
        self.assertEqual([x[0] for x in mutants][:13], [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38])
        self.assertEqual(len(mutants), pop.genotype().count(1))
        #
        mutants = list(pop.mutants(1))
        self.assertEqual(len(mutants), 80)
        self.assertEqual([x[0] for x in mutants][-13:], [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38])
        self.assertEqual(len(mutants), pop.genotype(1).count(1))
        #
        pop = Population(loci=[1, 2], size=[1, 2])
        arr = list(pop.mutants())
        self.assertEqual(len(arr), 0) 
        arr = list(pop.mutants(1))
        self.assertEqual(len(arr), 0)
        #
        pop.setGenotype([2, 3, 4])
        arr = list(pop.mutants())
        self.assertEqual(len(arr), pop.genoSize()*pop.popSize())
        arr = list(pop.mutants(1))
        self.assertEqual(len(arr), pop.genoSize()*pop.subPopSize(1))
        #
        pop.setGenotype([2, 0, 4, 0, 5])
        arr = list(pop.mutants())
        self.assertEqual(len(arr), 11)
        arr = list(pop.mutants(1))
        self.assertEqual(len(arr), 7)
        # set back
        pop.setGenotype([2, 0, 0, 0, 5])
        arr = list(pop.mutants())
        self.assertEqual(len(arr), 7)
        arr = list(pop.mutants(1))
        self.assertEqual(len(arr), 4)
        #
        
if __name__ == '__main__':
    unittest.main()


