#!/usr/bin/env python
#
# test initiailization operators
#
# # Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#


import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions, random

class TestInitialization(unittest.TestCase):

    def clearGenotype(self, pop):
        pop.setGenotype([0])

    def getGenotype(self, pop, atLoci=[], subPop=[], atPloidy=[]):
        geno = []
        if type(atPloidy) == type(1):
            ploidy = [atPloidy]
        elif len(atPloidy) > 0:
            ploidy = atPloidy
        else:
            ploidy = range(0, pop.ploidy())
        if len(atLoci) > 0:
            loci = atLoci
        else:
            loci = range(pop.totNumLoci())
        gs = pop.genoSize()
        tl = pop.totNumLoci()
        if len(subPop) > 0:
            for sp in subPop:
                for ind in pop.individuals(sp):
                        for p in ploidy:
                            for loc in loci:
                                geno.append(ind.allele(loc, p))
        else:
            arr = pop.genotype()
            if len(ploidy) == 0 and len(atLoci) == 0:
                geno = pop.arrGenotype(True)
            else:
                for i in range(pop.popSize()):
                    for p in ploidy:
                        for loc in loci:
                            geno.append( arr[ gs*i + p*tl +loc] )
        return geno

    def assertGenotype(self, pop, genotype,loci=[], subPop=[], atPloidy=[]):
        'Assert if the genotype of subPop of pop is genotype '
        geno = self.getGenotype(pop, loci, subPop, atPloidy)
        if AlleleType() == 'binary':
            if type(genotype) == type(1):
                self.assertEqual(geno, [genotype>0]*len(geno))
            else:
                self.assertEqual(geno, [x>0 for x in genotype])
        else:
            if type(genotype) == type(1):
                self.assertEqual(geno, [genotype]*len(geno))
            else:
                self.assertEqual(geno, genotype)

    def assertGenotypeFreq(self, pop, freqLow, freqHigh,loci=[], subPop=[], atPloidy=[]):
        'Assert if the genotype has the correct allele frequency'
        geno = self.getGenotype(pop, loci, subPop, atPloidy)
        if AlleleType() == 'binary':
            if len(freqLow) == 1:    # only one
                freq0 = geno.count(0)*1.0 / len(geno)
                assert freq0 >= freqLow[0] and freq0 <= freqHigh[0]
            else:    # 0 and 1, but group all other freq.
                f0 = [freqLow[0], sum(freqLow[1:])]
                f1 = [freqHigh[0], sum(freqHigh[1:])]
                freq0 = geno.count(0)*1.0 / len(geno)
                freq1 = geno.count(1)*1.0 / len(geno)
                assert freq0 >= f0[0] and freq0 <= f1[0]
                assert freq1 >= f0[1] and freq1 <= f1[1]
        else:    # all loci
            for i in range(len(freqLow)):
                freq = geno.count(i)*1.0 / len(geno)
                assert freq >= freqLow[i] and freq <= freqHigh[i]

    
    def testInitSex(self):
        'Testing operator initSex'
        pop = population(size=[500, 1000], loci=[1])
        InitSex(pop, sex=[Male, Female, Female])
        for sp in range(2):
            for idx, ind in enumerate(pop.individuals(sp)):
                if idx % 3 == 0:
                    self.assertEqual(ind.sex(), Male)
                else:
                    self.assertEqual(ind.sex(), Female)
        # maleFreq
        InitSex(pop, maleFreq=0.3)
        count = 0
        for ind in pop.individuals():
            if ind.sex() == Male:
                count += 1
        assert count * 1.0 / 1500 > 0.25 and count * 1.0 /1500 < 0.35
        # suBPop, virtual subPop
        pop = population(size=[500, 1000], loci=[1], infoFields=['x'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(infoSplitter('x', values=range(10, 15)))
        InitSex(pop, sex=[Male, Female, Female], subPops=[[0,0],[1,0]])
        for sp in range(2):
            for idx, ind in enumerate(pop.individuals([sp,0])):
                if idx % 3 == 0:
                    self.assertEqual(ind.sex(), Male)
                else:
                    self.assertEqual(ind.sex(), Female)

    def testInitByFreq(self):
        'Testing operator initByFreq '
        pop = population(size=[500, 1000, 500], loci=[2,4,2])
        # initialize all
        InitByFreq(pop, [.2, .3, .5])
        self.assertGenotypeFreq(pop, [.15, .25, .45],
            [.25, .35, .55])
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .3, .4, .1], loci=[2,4,6])
        self.assertGenotypeFreq(pop, [.15, .25, .35, .05],
            [.25, .35, .45, .15], loci=[2,4,6])
        self.assertGenotype(pop, 0, loci=[0,1,3,5,7])
        # use maleFreq=1 to avoid problem when comparing individuals
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .3, .4, .1], identicalInds=True,
            maleFreq=1)
        self.assertEqual(pop.individual(0), pop.individual(1))
        self.assertEqual(pop.individual(10), pop.individual(20))       
        # subPop
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .8], identicalInds=1, subPops=[0, 1],
            maleFreq=1)
        self.assertEqual(pop.individual(0), pop.individual(1))
        self.assertNotEqual(pop.individual(2), pop.individual(500))
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2],[.5,.5]])
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85], subPop=[0])
        self.assertGenotypeFreq(pop, [.75, .15], [.85, .25], subPop=[1])
        self.assertGenotypeFreq(pop, [.45, .45], [.55, .55], subPop=[2])
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .8], identicalInds=1, subPops=[0],
            maleFreq=1)
        self.assertEqual(pop.individual(6), pop.individual(7))
        self.assertNotEqual(pop.individual(0), pop.individual(500))
        self.assertGenotype(pop, 0, subPop=[1,2] )
        self.assertRaises(exceptions.ValueError,
            InitByFreq, pop, alleleFreq=[[.2, .8],[.8,.2]], identicalInds=1)
        # ploidy in initByFreq'
        pop = population(size=[500, 1000, 500], loci=[2,4,2])
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .3, .5], loci=[2,4,6], ploidy=[0])
        self.assertGenotypeFreq(pop, [.15, .25, .45], [.25, .35, .55],
            loci=[2,4,6], atPloidy=0)
        self.assertGenotype(pop, 0, loci=[0,3,5,7])
        self.assertGenotype(pop, 0, atPloidy=1)
        # virtual subPop
        pop = population(size=[500, 1000, 500], loci=[2,4,2], infoFields=['x'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(infoSplitter('x', values=range(10, 15)))
        InitByFreq(pop, [[0.2, 0.3, 0.5], [0.2,0.8], [0.5, 0.5]], 
            subPops=[[0,0],[1,1], [2,0]], loci=[2,4,6])
        self.assertGenotypeFreq(pop, [.15, .25, .45], [.25, .35, .55], 
            subPop=[[0,0]], loci=[2,4,6])
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85], 
            subPop=[[1,1]], loci=[2,4,6])
        self.assertGenotypeFreq(pop, [.45, .45], [.55, .55], 
            subPop=[[2,0]], loci=[2,4,6])
        # corner case
        self.clearGenotype(pop)
        InitByFreq(pop, [[0, 0, 1], [0, 1], [1]], subPops=[[0,0],[1,1], [2,1]])
        for ind in pop.individuals([0,0]):
            for allele in ind.genotype():
                 self.assertEqual(allele, 2) 
        for ind in pop.individuals([1,1]):
            for allele in ind.genotype():
                 self.assertEqual(allele, 1)
        for ind in pop.individuals([2,1]):
            for allele in ind.genotype():
                 self.assertEqual(allele, 0)
        self.clearGenotype(pop)
        self.assertRaises(exceptions.ValueError,InitByFreq, pop, alleleFreq=[-1,2])
 
    def testInitByValue(self):
        'Testing operator initByValue'
        pop = population(size=[500,1000, 500], loci=[2,4,2], infoFields=['x'])
        for ind in pop.individuals():
            ind.setInfo(random.randint(10, 20), 'x')
        pop.setVirtualSplitter(infoSplitter('x', values=range(10, 15)))
        # can initialize an invidiausl
        InitByValue(pop, [0]*5 + [2]*3 + [3]*5 +[4]*3)
        self.assertGenotype(pop, ([0]*5 + [2]*3 + [3]*5 +[4]*3)*pop.popSize())
        # one copy of chromosomes
        self.clearGenotype(pop)
        InitByValue(pop, [0]*5 + [7]*3)
        self.assertGenotype(pop, ([0]*5 + [7]*3)*(pop.popSize()*pop.ploidy()))
        # by proportion
        InitByValue(pop, value= [ [0]*8, [1]*8 ],
            proportions=[.3,.7])
        self.assertGenotypeFreq(pop, [0.25, 0.65], [0.35, 0.75])
        # ploidy
        self.clearGenotype(pop)
        InitByValue(pop, value=[0]*5 + [1]*3 , ploidy=[1])
        self.assertGenotype( pop, ([0]*5 + [1]*3)*pop.popSize(), atPloidy=1)
        self.assertGenotype(pop, 0, atPloidy=0)
        # subPop, virtual subPop
        self.clearGenotype(pop)
        InitByValue(pop, [0]*5 + [2]*3 + [3]*5 +[4]*3, subPops=[[0,1], [1,0]])
        self.assertGenotype(pop, ([0]*5 + [2]*3 + [3]*5 +[4]*3)*
            (pop.subPopSize([0,1])+pop.subPopSize([1,0])), subPop=[[0,1], [1,0]])
        # loci
        self.clearGenotype(pop)
        InitByValue(pop, value=[0,1,5], loci=[2,4,5], subPops=[[0,0], [2,1]])
        self.assertGenotype(pop, [0,1,5]*2*(pop.subPopSize([0,0])+
            pop.subPopSize([2,1])), loci=[2,4,5], subPop=[[0,0], [2,1]])
        self.assertGenotype(pop, 0, loci=[0,1,3,6,7])
        # error
        self.clearGenotype(pop)
        self.assertRaises(exceptions.ValueError,InitByValue, pop, [0]*16, ploidy=[0])    
        self.assertRaises(exceptions.ValueError,InitByValue, pop, [])  
        self.assertRaises(exceptions.ValueError,InitByValue, pop, [], subPops=[[0,1], [1,0]])  


if __name__ == '__main__':
    unittest.main()
