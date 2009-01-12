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
import unittest, os, sys, exceptions

def getGenotype(pop, atLoci=[], subPop=[], indRange=[], atPloidy=[]):
    '''HIDDEN
    Obtain genotype as specified by parameters

        atLoci
            subset of loci, default to all

        subPop
            subset of subpopulations, default ao all

        indRange
            individual ranges

    This is mostly used for testing purposes because the returned
    array can be large for large populations.
    '''
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
    if len(indRange) > 0:
        if type(indRange[0]) not in [type([]), type(())]:
            indRange = [indRange]
        arr = pop.genotype()
        for r in indRange:
            for i in range(r[0], r[1]):
                for p in ploidy:
                    for loc in loci:
                        geno.append( arr[ gs*i + p*tl + loc] )
    elif len(subPop) > 0:
        for sp in subPop:
            arr = pop.genotype(sp)
            for i in range(pop.subPopSize(sp)):
                for p in ploidy:
                    for loc in loci:
                        geno.append(arr[ gs*i + p*tl +loc])
    else:
        arr = pop.genotype()
        if len(ploidy) == 0 and len(atLoci) == 0:
            geno = pop.genotype()
        else:
            for i in range(pop.popSize()):
                for p in ploidy:
                    for loc in loci:
                        geno.append( arr[ gs*i + p*tl +loc] )
    return geno


class TestInitialization(unittest.TestCase):

    def clearGenotype(self, pop):
        pop.genotype()[:] = 0

    def assertGenotype(self, pop, genotype,
        loci=[], subPops=[], indRange=[], atPloidy=[]):
        'Assert if the genotype of subPop of pop is genotype '
        geno = getGenotype(pop, loci, subPops, indRange, atPloidy)
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

    def assertGenotypeFreq(self, pop, freqLow, freqHigh,
        loci=[], subPops=[], indRange=[], atPloidy=[]):
        'Assert if the genotype has the correct allele frequency'
        geno = getGenotype(pop, loci, subPops, indRange, atPloidy)
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
        #
        self.clearGenotype(pop)
        # use maleFreq=1 to avoid problem when comparing individuals
        InitByFreq(pop, [.2, .3, .4, .1], identicalInds=True,
            maleFreq=1)
        self.assertEqual(pop.individual(0), pop.individual(1))
        self.assertEqual(pop.individual(10), pop.individual(20))
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .8], identicalInds=1, subPops=[0],
            maleFreq=1)
        self.assertEqual(pop.individual(0), pop.individual(1))
        self.assertNotEqual(pop.individual(2), pop.individual(500))
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .8], subPops=[0])
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85], subPops=[0])
        self.assertGenotype(pop, 0, subPops=[1,2])
        #
        return
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .8], identicalInds=1, indRange=[6,9],
            maleFreq=1)
        self.assertEqual(pop.individual(6), pop.individual(7))
        self.assertNotEqual(pop.individual(0), pop.individual(7))
        #
        self.assertRaises(exceptions.ValueError,
            InitByFreq, pop, alleleFreq=[[.2, .8],[.8,.2]])
        #
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2],[.5,.5]])
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85], subPops=[0])
        self.assertGenotypeFreq(pop, [.75, .15], [.85, .25], subPops=[1])
        self.assertGenotypeFreq(pop, [.45, .45], [.55, .55], subPops=[2])
        #
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]], indRange=[[0,300],[500,1300]])
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85], indRange=[0,300] )
        self.assertGenotypeFreq(pop, [.75, .15], [.85, .25], indRange=[500, 1300] )
        self.assertGenotype(pop, 0, indRange=[[300,500], [1300, 2000]] )
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .8], identicalInds=1, subPops=[0],
            maleFreq=1)
        self.assertEqual(pop.individual(6), pop.individual(7))
        self.assertNotEqual(pop.individual(0), pop.individual(500))
        self.assertGenotype(pop, 0, subPops=[1,2] )
        #
        self.assertRaises(exceptions.ValueError,
            InitByFreq, pop, alleleFreq=[[.2, .8],[.8,.2]], identicalInds=1)
        #
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2],[1]], identicalInds=1,
            maleFreq=1)
        self.assertEqual(pop.individual(6), pop.individual(7))
        self.assertEqual(pop.individual(6,1), pop.individual(7,1))
        self.assertEqual(pop.individual(6,2), pop.individual(7,2))
        self.assertNotEqual(pop.individual(0), pop.individual(500))
        #
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
            indRange=[[0,300],[500,1300]], loci=[2,3,5,6])
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85],
            indRange=[0,300], loci=[2,3,5,6] )
        self.assertGenotypeFreq(pop, [.75, .15], [.85, .25],
            indRange=[500, 1300], loci=[2,3,5,6] )
        self.assertGenotype(pop, 0, indRange=[[300,500], [1300, 2000]] )
        self.assertGenotype(pop, 0, indRange=[[300,500],[1300,2000]], loci=[2,3,5,6])
        #
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
            indRange=[[0,4],[8,15]], loci=[2,3,5,6],
            identicalInds=1, maleFreq=1)
        self.assertEqual(pop.individual(0), pop.individual(3))
        self.assertEqual(pop.individual(8), pop.individual(14))
        self.assertNotEqual(pop.individual(0), pop.individual(8))
        #
        #Testing atPloidy in initByFreq'
        pop = population(size=[500, 1000, 500], loci=[2,4,2])
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .3, .5], loci=[2,4,6], atPloidy=0)
        self.assertGenotypeFreq(pop, [.15, .25, .45], [.25, .35, .55],
            loci=[2,4,6], atPloidy=0)
        self.assertGenotype(pop, 0, loci=[0,3,5,7])
        self.assertGenotype(pop, 0, atPloidy=1)
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .3, .5], atPloidy=1)
        self.assertGenotypeFreq(pop, [.15, .25, .45], [.25, .35, .55],
            atPloidy=1)
        self.assertGenotype(pop, 0, atPloidy=0)
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .3, .4, .1], identicalInds=1, atPloidy=0,
            maleFreq=1)
        self.assertEqual(pop.individual(0), pop.individual(1))
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .8], subPops=[0], atPloidy=1)
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85],
            atPloidy=1, subPops=[0])
        self.assertGenotype(pop, 0, atPloidy=0)
        self.assertGenotype(pop, 0, subPops=[1])
        #
        self.clearGenotype(pop)
        InitByFreq(pop, [.2, .8], identicalInds=1, indRange=[0,1000], atPloidy=0,
            maleFreq=1)
        self.assertEqual(pop.individual(0), pop.individual(1))
        self.assertNotEqual(pop.individual(0), pop.individual(1001))
        #
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2],[.5,.5]], atPloidy=1)
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85], subPops=[0], atPloidy=1)
        self.assertGenotypeFreq(pop, [.75, .15], [.85, .25], subPops=[1], atPloidy=1)
        self.assertGenotypeFreq(pop, [.45, .45], [.55, .55], subPops=[2], atPloidy=1)
        self.assertGenotype(pop, 0, atPloidy=0)
        #
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
                             indRange=[[0,300],[500,1300]], atPloidy=1)
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85], indRange=[0,300], atPloidy=1 )
        self.assertGenotypeFreq(pop, [.75, .15], [.85, .25], indRange=[500, 1300], atPloidy=1 )
        self.assertGenotype(pop, 0, indRange=[[300,500], [1300, 2000]] )
        self.assertGenotype(pop, 0, atPloidy=0)
        #
        self.clearGenotype(pop)
        InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
            indRange=[[0,300],[500,1300]], loci=[2,3,5,6], atPloidy=0)
        self.assertGenotypeFreq(pop, [.15, .75], [.25, .85], loci=[2,3,5,6],
            indRange=[0,300], atPloidy=0 )
        self.assertGenotypeFreq(pop, [.75, .15], [.85, .25], loci=[2,3,5,6],
            indRange=[500, 1300], atPloidy=0 )
        self.assertGenotype(pop, 0, indRange=[[300,500], [1300, 2000]] )
        self.assertGenotype(pop, 0, atPloidy=1)

    def testInitByValue(self):
        'Testing operator initByValue'
        pop = population(size=[500,1000, 500], loci=[2,4,2])
        # can initialize an invidiausl
        InitByValue(pop, [0]*5 + [2]*3 + [3]*5 +[4]*3)
        self.assertGenotype(pop, ([0]*5 + [2]*3 + [3]*5 +[4]*3)*pop.popSize())
        #
        self.clearGenotype(pop)
        # or one copy of chromosomes
        InitByValue(pop, [0]*5 + [7]*3)
        self.assertGenotype(pop, ([0]*5 + [7]*3)*(pop.popSize()*pop.ploidy()))
        #
        self.clearGenotype(pop)
        return 
        InitByValue(pop, [0]*5 + [2]*3 + [3]*5 +[4]*3,
            indRange=[[2],[5]])
        self.assertGenotype(pop, ([0]*5 + [2]*3 + [3]*5 +[4]*3)*2,
            indRange=[[2,3],[5,6]])
        self.assertGenotype(pop, 0, indRange=[3,4])
        #
        self.clearGenotype(pop)
        InitByValue(pop, value=[0,1,5], loci=[2,4,5],
            indRange=[3,6])
        self.assertGenotype(pop, [0,1,5]*6, loci=[2,4,5],
            indRange=[3,6])
        self.assertGenotype(pop, 0, loci=[0,1,3,6,7])
        self.assertGenotype(pop, 0, indRange=[6, pop.popSize()-1] )
        # by proportion
        InitByValue(pop, value= [ [0]*8, [1]*8 ],
            proportions=[.3,.7])
        self.assertGenotypeFreq(pop, [0.25, 0.65], [0.35, 0.75])
        # atPloidy
        self.clearGenotype(pop)
        InitByValue(pop, value=[0]*5 + [1]*3 , atPloidy=1)
        self.assertGenotype( pop, ([0]*5 + [1]*3)*pop.popSize(), atPloidy=1)
        self.assertGenotype(pop, 0, atPloidy=0)
        # error if.
        self.clearGenotype(pop)
        self.assertRaises(exceptions.ValueError,
            InitByValue, pop, [0]*16, atPloidy=0)
        #
        self.clearGenotype(pop)
        InitByValue(pop, [0]*5 + [2]*3,
            indRange=[[2,4],[5,7]], atPloidy=1)
        self.assertGenotype(pop, ([0]*5 + [2]*3)*4, atPloidy=1, indRange=[[2,4],[5,7]])
        # whole ind
        self.clearGenotype(pop)
        InitByValue(pop, value=[[0]*3, [1]*3], loci=[2,4,5],
            proportions=[.3,.7], indRange=[[300,600],[700,1000]] )
        self.assertGenotype(pop, 0, loci=[0,1,3,6,7])
        self.assertGenotype(pop, 0, indRange=[[0,300],[600,700]])
        self.assertGenotypeFreq(pop, [0.25, 0.65], [0.35, 0.75],
            loci=[2,4,5], indRange=[[300,600],[700,1000]])

    def testInitSex(self):
        'Testing operator initSex'
        #
        pop = population(size=[500, 1000], loci=[1])
        InitSex(pop, sex=[Male, Female, Female])
        for sp in range(2):
            for idx, ind in enumerate(pop.individuals(sp)):
                if idx % 3 == 0:
                    self.assertEqual(ind.sex(), Male)
                else:
                    self.assertEqual(ind.sex(), Female)
        #
        InitSex(pop, maleFreq=0.3)
        count = 0
        for ind in pop.individuals():
            if ind.sex() == Male:
                count += 1
        assert count * 1.0 / 1500 > 0.25 and count * 1.0 /1500 < 0.35

if __name__ == '__main__':
    unittest.main()
