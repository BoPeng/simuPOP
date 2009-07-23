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
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions, random, copy

class TestPopulation(unittest.TestCase):
    # define a few functions to create basic populations
    def getPop(self, VSP=False, size=[20, 80], loci = [1, 2], infoFields=['x'],
            ancGen=0, *arg, **kwargs):
        pop = population(size=size, ploidy=2, loci=loci, infoFields=infoFields,
            ancGen=ancGen, *arg, **kwargs)
        pop.setGenotype([random.randint(1, 5) for x in range(pop.popSize()*pop.ploidy())])
        for info in infoFields:
            pop.setIndInfo([random.random() for x in range(pop.popSize())], info)
        for i in range(ancGen):
            pop.push(self.getPop(size=size, loci=loci, infoFields=infoFields, ancGen=0, *arg, **kwargs))
        InitSex(pop)
        if VSP:
            pop.setVirtualSplitter(sexSplitter())
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
        Stat(pop, numOfMale=True, vars=['numOfMale_sp', 'numOfFemale_sp'])
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
        self.assertEqual(pop.subPopName(0), "")
        self.assertEqual(pop.subPopName([0, 0]), "x = 10")
        self.assertEqual(pop.subPopName([0, 1]), "x = 11")
        self.assertEqual(pop.subPopName([0, 4]), "x = 14")
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

    def testSubPopName(self):
        'Testing population::setSubPopName(name, subPop), subPopByName(subPop)'
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
        self.assertRaises(exceptions.ValueError, pop.subPopByName, 'D')

    def testIndividuals(self):
        'Testing function population::individuals(), individuals(subPop), individual(idx, subPop=0)'
        def testAllInd(pop):
            self.assertEqual(len(list(pop.individuals())), pop.popSize())
            self.assertEqual(len(list(pop.individuals(0))), pop.subPopSize(0))
            self.assertEqual(len(list(pop.individuals(1))), pop.subPopSize(1))
        testAllInd(self.getPop())
        testAllInd(self.getPop(True))
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
        if AlleleType() == 'binary':
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
        'Testing population::ancestor(idx, gen), ancestor(idx, gen, subPop), push(pop)'
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
        self.assertRaises(exceptions.IndexError, pop.ancestor, 2, 10000)
        self.assertRaises(exceptions.IndexError, pop.ancestor, 3, 10000)
        for idx, ind in enumerate(pop_c.individuals(0)):
            self.assertEqual(ind, pop.ancestor(idx, 1, 0))

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
        # setAncestralDepth
        pop = self.getPop(ancGen = 5)
        pop.setAncestralDepth(3)
        self.assertEqual(pop.ancestralGens(), 3)

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

    def testAddChromFrom(self):
        'Testing population::addChromFrom(pop)'
        pop = population(size=100, ploidy=2, loci=[1, 2], chromNames=["c1", "c2"], lociNames = ['l1', 'l2', 'l3'])
        pop2 = pop.clone()
        pop1 = population(size=100, ploidy=2, loci=[2, 3], chromNames=["c3", "c4"],
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
        pop = population(size=100, ploidy=2, loci=[1, 2])
        pop1 = population(size=100, ploidy =2, loci=[1, 2])
        self.assertRaises(exceptions.ValueError, pop.addChromFrom, pop1)
        pop2 = population(size=200, ploidy=2, loci=[2, 3], chromNames=["c3", "c4"],
            lociNames = ['l4', 'l5', 'l6', 'l7', 'l8'])
        self.assertRaises(exceptions.ValueError, pop.addChromFrom, pop2)

    def testAddIndFrom(self):
        'Testing population::addIndFrom(pop)'
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
        self.assertRaises(exceptions.ValueError, pop.addIndFrom, pop1)
        pop1 = population(size=100, ploidy=2, loci=[1, 2, 3])
        # different genotype structure
        self.assertRaises(exceptions.ValueError, pop.addIndFrom, pop1)

    def testAddLociFrom(self):
        'Testing population::addLociFrom(pop)'
        pop = self.getPop(chromNames=["c1", "c2"], ancGen=5, lociPos=[[1], [2, 5]], lociNames = ['l1', 'l2', 'l3'])
        pop1 = pop.clone()
        pop2 = self.getPop(chromNames=["c3", "c4"], ancGen=5, lociPos=[[4], [3, 6]], lociNames = ['l4', 'l5', 'l6'])
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
                # src: the source population
                # j: index in source population
                for i, src, j in [(0, 0, 0), (1, 1, 0), (2, 0, 1), (3, 1, 1), (4, 0, 2), (5, 1, 2)]:
                    for p in range(pop.ploidy()):
                        self.assertEqual(ind.allele(i, p), inds[src].allele(j, p))

    def testAddLoci(self):
        'Testing population::addLoci(chrom, pos, names=[])'
        pop = self.getPop(size = 100, chromNames=["c1", "c2"], ancGen=5, lociPos=[[1], [3, 5]], lociNames = ['l1', 'l2', 'l3'])
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
                # j: index in source population
                for i, j in [(0, 0), (2, 1), (3, 2)]:
                    for p in range(pop.ploidy()):
                        self.assertEqual(ind.allele(i, p), ind1.allele(j, p))
                for k in newpos:
                    self.assertEqual(ind.allele(k), 0)
        self.assertRaises(exceptions.ValueError, pop.addLoci, [2], [8], ['l7'])

    def testDeepcopy(self):
        'Testing deepcopy of population'
        pop = self.getPop(False, ancGen=3)
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


    def comparePop(self, pop1, pop2, inds, loci=None, infoFields=None, gen=0):
        #inds should be original pop1 index of individuals in pop2
        pop1.useAncestralGen(gen)
        pop2.useAncestralGen(gen)
        if loci is None:
            for idx, ind in enumerate(inds):
                if infoFields is None:
                    self.assertEqual(pop1.individual(ind), pop2.individual(idx))
                else:
                    self.assertEqual(pop1.individual(ind).genotype(), pop2.individual(idx).genotype())
                    for info in infoFields:
                        self.assertEqual(pop1.individual(ind).info(info), pop2.individual(idx).info(info))
        else:
            for idx, ind in enumerate(inds):
                if infoField is None:
                    self.assertEqual(pop1.individual(ind), pop2.individual(idx))
                else:
                    self.assertEqual(pop1.individual(ind).genotype(), pop2.individual(idx).genotype())
                    for info in infoFields:
                        self.assertEqual(pop1.individual(ind).info(info), pop2.individual(idx).info(info))
                for idx1, loc in enumerate(loci):
                    for p in range(pop1.ploidy()):
                        self.assertEqual(pop1.individual(ind).allele(loc, p), pop2.individual(idx).allele(idx1, p))
            
        # for idx, ind in enumerate(inds):
        #    self.assertEqual(pop1.individual(ind), pop2.individual(idx))

        
               
    def newtestExtract(self):
        pop1 = self.getPop(size=[200, 400, 200], loci=[20, 40, 20], ancGen=3, infoFields=['x'])
        for i in pop1.individuals():
            n = random.choice([-1, -3, -5, 0, 1, 2, 3, 4, 5])
            i.setInfo(n, 'x')
        pop2 = pop1.extract(field='x')
        for gen in range(4):
            pop1.useAncestralGen(gen)
            pop2.useAncestralGen(gen)
            inds = []
            for j in range(0,6):
                for n in range(pop1.popSize()):
                    ind = pop1.individual(n)
                    if ind.info('x') == j:
                        inds.append(n)
            self.comparePop(pop1, pop2, inds, infoFields=['x'], gen=gen)
        #if ancGen == -1:
        #    self.assertEqual(pop1.ancestralGens(), pop2.ancestralGens())
        #    topGen = pop1.ancestralGens()
        #else:
        #    self.assertEqual(pop2.ancestralGens(), ancGen)
        #    topGen = ancGen
        #for gen in range(topGen):
        #    self.compare(pop1, pop2, gen=gen)
        

    def testExtract(self):
        'Testing population::Extract(field=None, loci=None, infoFields=None, ancGen =-1)'
        # If subpoulation size is too small, the last subpopulation
        # may not have any individual.
        pop = population(size=[30, 50], loci=[2, 3],  infoFields=['x', 'y'])
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
                if AlleleType() == 'binary':
                    self.assertEqual(ind.genotype(), [1]*(pop.totNumLoci()*pop.ploidy()))
                else:
                    self.assertEqual(ind.genotype(), [sp+1]*(pop.totNumLoci()*pop.ploidy()))
        # FIXME: test extraction of loci and ancestral generations

    def testMergeSubPops(self):
        'Testing population::mergeSubPops(subpops=[])'
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

    def testRemoveSubPops(self):
        'Testing population::removeEmptySubPops(), removeSubPops()'
        pop = self.getPop(size=[0, 100, 0, 20, 30, 0, 50], subPopNames=['A', 'B', 'C', 'D', 'E', 'F', 'G'])
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
        self.assertRaises(exceptions.IndexError, pop.removeSubPops, [8])
        # accept single input
        pop.removeSubPops(0)
        
    def testRemoveIndividuals(self):
        'Testing population::removeIndividuals(inds)'
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
        self.assertRaises(exceptions.IndexError, pop.removeIndividuals, 500)
        # 2) pop.removeIndividuals([]) should not change anything (self.assertEqual(pop, pop1))
        pop = pop1.clone()
        pop.removeIndividuals([])
        self.assertEqual(pop, pop1)
        # 3) pop.removeIndividuals(range(15, 25)) ...
        pop = pop1.clone()
        inds = range(15, 25)
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
        inds = range(15, 125)
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
        inds = range(pop.subPopBegin(1), pop.subPopEnd(1))
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
        inds = range(0, 150)
        random.shuffle(inds)
        pop.removeIndividuals(inds)
        self.assertEqual(pop.subPopSizes(), (0, 0, 0))
        self.assertEqual(pop.subPopNames(), pop1.subPopNames())

    def testRemoveLoci(self):
        'Testing population::removeLoci(loci=[], keep=[])'
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

    def testResize(self):
        'Testing population::resize(newSubPopSizes, propagate=false)'
        pop = self.getPop(size=[100, 20, 30], loci=[4, 5, 1])
        InitByFreq(pop, [.2, .3, .5])
        pop1 = pop.clone()
        pop2 = pop.clone()
        # number of subpop mismatch
        self.assertRaises(exceptions.ValueError, pop1.resize, [50, 50])
        # resize without propagation
        pop1.resize([50, 50, 80], propagate=False)
        for sp in range(pop1.numSubPop()):
            for i in range(min(pop1.subPopSize(sp), pop.subPopSize(sp))):
                self.assertEqual(pop1.individual(i, sp), pop.individual(i, sp))
            for i in range(min(pop1.subPopSize(sp), pop.subPopSize(sp)), pop1.subPopSize(sp)):
                self.assertEqual(pop1.individual(i, sp).genotype(), [0]*20)
        # resize with population
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
        'Testing population::splitSubPop(subPop, sizes)'
        pop = population(size=[100, 80, 50], subPopNames=['A', 'B', 'C'])
        pop1 = pop.clone()
        self.assertRaises(exceptions.ValueError, pop.splitSubPop, 1, [20, 70])
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
        pop = population(size=[100, 80, 50])
        self.assertRaises(exceptions.ValueError, pop.splitSubPop, 1, [20, 70], names=['A1'])
        ids = pop.splitSubPop(1, [20, 60], names=['A1', 'A2'])
        self.assertEqual(ids, (1, 2))
        self.assertEqual(pop.subPopName(1), 'A1')
        self.assertEqual(pop.subPopName(2), 'A2')
        self.assertEqual(pop.subPopNames(), ('', 'A1', 'A2', ''))

    def testSetSubPopByIndInfo(self):
        'Testing population::setSubPopByIndInfo(field)'
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

    def testAddInfoFields(self):
        'Testing population::addInfoFields(fields, init=0)'
        pop = self.getPop()
        pop1 = pop.addInfoFields('fitness', 5.0)
        self.assertEqual(pop.infoSize(), 2)
        pop1 = pop.addInfoFields('misc', 6.0)
        self.assertEqual(pop.infoSize(), 3)
        self.assertEqual(pop.indInfo('fitness'), tuple([5.0]*pop.popSize()))
        self.assertEqual(pop.indInfo('misc'), tuple([6.0]*pop.popSize()))
        pop = self.getPop()
        pop1 = pop.addInfoFields(['x', 'fitness', 'misc'],  2.0)
        self.assertEqual(pop.infoSize(), 3)
        self.assertEqual(pop.indInfo('fitness'), tuple([2.0]*pop.popSize()))
        self.assertEqual(pop.indInfo('misc'), tuple([2.0]*pop.popSize()))
        # info field x is re-initialized
        self.assertEqual(pop.indInfo('x'), tuple([2.0]*pop.popSize()))


    def testIndInfo(self):
        'Testing population::indInfo(idx), indInfo(name), indInfo(idx, subPop)'
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
            self.assertEqual(pop.indInfo('x', 0), tuple([1, 2]*(pop.subPopSize(0)/2)))
            self.assertEqual(pop.indInfo(0, 1), tuple([3, 4]*(pop.subPopSize(1)/2)))
        #
        testSetAndRead(self.getPop())
        testSetAndRead(self.getPop(True))
        # test for virtual subpopulation
        def testVSPSetAndRead(pop):
            pop.setIndInfo([1, 2], 'x', [1, 0])
            pop.setIndInfo([3], 0, [1, 1])
            for idx, ind in enumerate(pop.individuals([1, 0])):
                self.assertEqual(ind.sex(), Male)
                if idx % 2 == 0:
                    self.assertEqual(ind.info('x'), 1)
                else:
                    self.assertEqual(ind.info('x'), 2)
            for idx, ind in enumerate(pop.individuals([1, 1])):
                self.assertEqual(ind.sex(), Female)
                self.assertEqual(ind.info('x'), 3)

            self.assertEqual(pop.indInfo('x', [1, 0]), tuple(([1, 2]*pop.subPopSize(1))[:pop.subPopSize([1, 0])]))
            self.assertEqual(pop.indInfo(0, [1, 1]), tuple([3]*pop.subPopSize([1, 1])))
        #
        self.assertRaises(exceptions.IndexError, testVSPSetAndRead, self.getPop())
        testVSPSetAndRead(self.getPop(True))

    def testSetInfoFields(self):
        'Testing population::setInfoFields(fields, init=0)'
        pop = self.getPop()
        pop1 = pop.setInfoFields(['fitness', 'misc'],  3)
        self.assertEqual(pop.infoSize(), 2)
        # info field x is removed
        self.assertEqual(pop.indInfo('fitness'), tuple([3]*pop.popSize()))
        self.assertEqual(pop.indInfo('misc'), tuple([3]*pop.popSize()))

    def testUpdateInfoFieldsFrom(self):
        'Testing population::updateInfoFieldsFrom(fields, pop, fromFields=[], ancGen=-1)'
        pop = self.getPop(size = 100, ancGen = 5)
        pop1 = self.getPop(size = 100, ancGen = 5)
        for gen in range(0, 6):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            self.assertNotEqual(pop.indInfo('x'), pop1.indInfo('x'))
        pop.updateInfoFieldsFrom('x', pop1)
        for gen in range(0, 6):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            self.assertEqual(pop.indInfo('x'), pop1.indInfo('x'))
        # do not update all ancestral generations
        pop1 = self.getPop(size = 100, ancGen = 5)
        pop.updateInfoFieldsFrom('x', pop1, ancGen=2)
        for gen in range(0, 6):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            if gen <= 2:
                self.assertEqual(pop.indInfo('x'), pop1.indInfo('x'))
            else:
                self.assertNotEqual(pop.indInfo('x'), pop1.indInfo('x'))
        pop3 = self.getPop(size = 200)
        self.assertRaises(exceptions.ValueError, pop.updateInfoFieldsFrom, 'x', pop3)

    def testClone(self):
        'Testing population::clone()'
        pop = self.getPop(ancGen = 5)
        pop1 = pop.clone()
        for gen in range(pop.ancestralGens(), -1, -1):
            pop.useAncestralGen(gen)
            pop1.useAncestralGen(gen)
            self.assertEqual(pop, pop1)

    def testSave(self):
        'Testing population::save(filename)'
        pop = self.getPop()
        pop.save("popout");
        pop1 = LoadPopulation("popout")
        self.assertEqual(pop, pop1)
        os.remove('popout')

    def testVars(self):
        'Testing population::vars(), vars(subPop), dvars(), dvars(subPop)'
        pop = self.getPop(size=1000, loci=[2, 4])
        InitByFreq(pop, [.2, .3, .5])
        Stat(pop, alleleFreq=range(0, 6))
        pop1 = pop.clone()
        self.assertEqual(len(pop.vars()["alleleFreq"]), 6)
        self.assertEqual(len(pop.dvars().alleleFreq), 6)
        self.assertEqual(len(pop1.vars()["alleleFreq"]), 6)
        self.assertEqual(len(pop1.dvars().alleleFreq), 6)
        # with subPop
        pop = self.getPop(size=[20, 80], loci=[2, 4])
        InitByFreq(pop, [.2, .3, .5])
        Stat(pop, alleleFreq=range(0, 6), vars='alleleFreq_sp')
        pop1 = pop.clone()
        self.assertEqual(len(pop.vars(0)["alleleFreq"]), 6)
        self.assertEqual(len(pop.dvars(1).alleleFreq), 6)
        self.assertEqual(len(pop1.vars(0)["alleleFreq"]), 6)
        self.assertEqual(len(pop1.dvars(1).alleleFreq), 6)

    def testSexSplitter(self):
        'Testing sexSplitter::sexSplitter()'
        pop = population(size=[20, 80])
        InitByFreq(pop, [0.4, 0.6])
        Stat(pop, numOfMale=True, vars=['numOfMale_sp', 'numOfFemale_sp'])
        pop.setVirtualSplitter(sexSplitter())
        self.assertEqual(pop.subPopSize([1, 0]), pop.dvars(1).numOfMale)
        self.assertEqual(pop.subPopSize([1, 1]), pop.dvars(1).numOfFemale)
        self.assertEqual(pop.subPopName([1, 0]), 'Male')
        self.assertEqual(pop.subPopName([1, 1]), 'Female')
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.sex(), Male)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.sex(), Female)
        numMale = 0
        numFemale = 0
        for ind in pop.individuals(1):
            if ind.sex() == Male:
                numMale += 1
            else:
                numFemale += 1
        #print numMale, numFemale
        self.assertEqual(numMale == 0, False)
        self.assertEqual(numFemale == 0, False)

    def testAffectionSplitter(self):
        'Testing affectionSplitter::affectionSplitter()'
        pop = population(size=[20, 80], loci=[1, 2])
        InitByFreq(pop, [0.4, 0.6])
        MaPenetrance(pop, loci=0, wildtype=0, penetrance=[0.2, 0.4, 0.8])
        Stat(pop, numOfAffected=True, vars=['numOfAffected_sp', 'numOfUnaffected_sp'])
        pop.setVirtualSplitter(affectionSplitter())
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
        'Testing infoSplitter::infoSplitter(field, values=[], cutoff=[])'
        pop = population(1000, infoFields=['x'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(infoSplitter('x', values=range(10, 15)))
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
                self.assertEqual(ind.intInfo('x'), 10+i)
        # test cutoff
        pop.setVirtualSplitter(infoSplitter('x', cutoff=[11.5, 13.5]))
        self.assertEqual(pop.subPopName([0, 0]), "x < 11.5")
        self.assertEqual(pop.subPopName([0, 1]), "11.5 <= x < 13.5")
        self.assertEqual(pop.subPopName([0, 2]), "x >= 13.5")
        self.assertEqual(pop.subPopSize([0, 0]), infos.count(10) + infos.count(11))
        self.assertEqual(pop.subPopSize([0, 1]), infos.count(12) + infos.count(13))
        self.assertEqual(pop.subPopSize([0, 2]),
            sum([infos.count(x) for x in range(14, 21)]))
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.intInfo('x') < 11.5, True)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(11.5 <= ind.intInfo('x') < 13.5, True)
        for ind in pop.individuals([0, 2]):
            self.assertEqual(ind.intInfo('x') >=13.5, True)

    def testProportionSplitter(self):
        'Testing proportionSplitter::proportionSplitter(proportions=[])'
        pop = population(10)
        pop.setVirtualSplitter(proportionSplitter([0.01]*100))
        for i in range(100):
            self.assertEqual(pop.subPopName([0, i]), "Prop 0.01")
            if i != 99:
                self.assertEqual(pop.subPopSize([0, i]), 0)
            else:
                # the last vsp is specially treated to avoid such problem.
                self.assertEqual(pop.subPopSize([0, i]), 10)
        #
        pop = population(1000)
        pop.setVirtualSplitter(proportionSplitter([0.4, 0.6]))
        self.assertEqual(pop.subPopSize([0, 0]), 400)
        self.assertEqual(pop.subPopSize([0, 1]), 600)

    def testRangeSplitter(self):
        'Testing rangeSplitter::rangeSplitter(ranges)'
        pop = population(100)
        pop.setVirtualSplitter(rangeSplitter(ranges=[[10, 20], [80, 200]]))
        self.assertEqual(pop.subPopName([0, 0]), "Range [10, 20)")
        self.assertEqual(pop.subPopName([0, 1]), "Range [80, 200)")
        self.assertEqual(pop.subPopSize([0, 0]), 10)
        self.assertEqual(pop.subPopSize([0, 1]), 20)

    def testGenotypeSplitter(self):
        'Testing genotypeSplitter::genotypeSplitter(loci(or locus), alleles, phase=False)'
        pop = population(1000, loci=[2, 3])
        InitByFreq(pop, [0.3, 0.7])
        pop.setVirtualSplitter(genotypeSplitter(loci=1, alleles=[[0, 0], [1, 0]], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 0")
        self.assertEqual(pop.subPopName([0, 1]), "Genotype 1: 1 0")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0), 0)
            self.assertEqual(ind.allele(1, 1), 0)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0), 1)
            self.assertEqual(ind.allele(1, 1), 0)
        Stat(pop, genoFreq=[1])
        self.assertEqual(pop.subPopSize([0, 0]), pop.dvars().genoNum[1][(0,0)])
        self.assertEqual(pop.subPopSize([0, 1]), pop.dvars().genoNum[1][(1,0)])
        # non-phased case
        pop.setVirtualSplitter(genotypeSplitter(loci=1, alleles=[[0, 0], [1, 0]], phase=False))
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0), 0)
            self.assertEqual(ind.allele(1, 1), 0)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0)==0 or ind.allele(1, 0)==1, True)
            self.assertEqual(ind.allele(1, 1)==0 or ind.allele(1, 1)==1, True)
        # multiple loci
        pop.setVirtualSplitter(genotypeSplitter(loci=[0, 1], alleles=[0, 1, 1, 1], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 0, 1: 0 1 1 1")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(0, 0), 0)
            self.assertEqual(ind.allele(0, 1), 1)
            self.assertEqual(ind.allele(1, 0), 1)
            self.assertEqual(ind.allele(1, 1), 1)
        # multiple genotype at the same locus
        pop.setVirtualSplitter(genotypeSplitter(loci=1, alleles=[0, 1, 1, 1], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 1 1 1")
        for ind in pop.individuals([0, 0]):
             self.assertEqual(ind.allele(1, 0)==1  or ind.allele(1, 0)==0, True)
             self.assertEqual(ind.allele(1, 1), 1)
        # haploid case
        pop = population(1000, ploidy = 1, loci=[2, 3])
        InitByFreq(pop, [0.3, 0.7])
        pop.setVirtualSplitter(genotypeSplitter(loci=1, alleles=[[0, 1], [2]], phase=True))
        self.assertEqual(pop.subPopName([0, 0]), "Genotype 1: 0 1")
        self.assertEqual(pop.subPopName([0, 1]), "Genotype 1: 2")
        for ind in pop.individuals([0, 0]):
            self.assertEqual(ind.allele(1, 0)==1 or ind.allele(1, 0)==0, True)
        for ind in pop.individuals([0, 1]):
            self.assertEqual(ind.allele(1, 0), 2)

    def testCombinedSplitter(self):
        'Testing combinedSplitter:: combinedSplitter(splitters=[])'
        pop = population(1000, loci=[2, 3])
        InitByFreq(pop, [0.3, 0.7])
        pop.setVirtualSplitter(combinedSplitter([
            genotypeSplitter(loci=1, alleles=[[0, 0], [1, 0]], phase=True),
            sexSplitter()]))
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
            self.assertEqual(ind.sex(), Male)
        for ind in pop.individuals([0, 3]):
            self.assertEqual(ind.sex(), Female)
        Stat(pop, numOfMale=True, vars='numOfFemale_sp')
        self.assertEqual(pop.subPopSize([0, 3]), pop.dvars(0).numOfFemale)


if __name__ == '__main__':
    unittest.main()


