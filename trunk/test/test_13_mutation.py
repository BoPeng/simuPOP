#!/usr/bin/env python
#
# Purpose:
#    testing of interfaces of mutators    of simuPOP
#
# Author:
#    Bo Peng (bpeng@rice.edu)
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


class TestMutator(unittest.TestCase):

    def assertGenotype(self, pop, genotype,
        loci=[], subPops=[], indRange=[], atPloidy=[]):
        'Assert if the genotype of subPop of pop is genotype '
        geno = getGenotype(pop, loci, subPops, indRange, atPloidy)
        if moduleInfo()['alleleType'] == 'binary':
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
        'Testing if the genotype has the correct allele frequency'
        geno = getGenotype(pop, loci, subPops, indRange, atPloidy)
        if moduleInfo()['alleleType'] == 'binary':
            if len(freqLow) == 1:    # only one
                freq0 = geno.count(0)*1.0 / len(geno)
                self.assertTrue(freq0 >= freqLow[0] and freq0 <= freqHigh[0])
            else:
                f0 = [freqLow[0], sum(freqLow[1:])]
                f1 = [freqHigh[0], sum(freqHigh[1:])]
                freq0 = geno.count(0)*1.0 / len(geno)
                freq1 = geno.count(1)*1.0 / len(geno)
                #print f0,f1,freq0,freq1
                self.assertTrue(freq0 >= f0[0] and freq0 <= f1[0])
                self.assertTrue(freq1 >= f0[1] and freq1 <= f1[1])
        else:
            for i in range(len(freqLow)):
                freq = geno.count(i)*1.0 / len(geno)
                self.assertTrue(freq >= freqLow[i])
                self.assertTrue(freq <= freqHigh[i])

    def testSP(self):
        'Testing the subpop support of mutators'
        pop = population(size=[1000, 2000], ploidy=2, loci=[2, 3])
        simu = simulator(pop )
        simu.evolve(initOps = [initSex()],
            matingScheme = randomMating(),
            postOps = [ KamMutator(k=2, rates=0.5, loci=[1,4],
                subPops=1)],
            gen = 100)
        pop = simu.extract(0)
        Stat(pop, alleleFreq=range(5), vars='alleleFreq_sp')
        for loc in [0, 2, 3]:
            self.assertEqual(pop.dvars(0).alleleFreq[loc][0], 1.0)
            self.assertEqual(pop.dvars(1).alleleFreq[loc][0], 1.0)
        for loc in [1, 4]:
            self.assertEqual(pop.dvars(0).alleleFreq[loc][0], 1.0)
            self.assertNotEqual(pop.dvars(1).alleleFreq[loc][0], 1.0)

    def testVSP(self):
        'Testing the subpops parameter of mutator'
        pop = population(size=2000, ploidy=2, loci=[2, 3])
        pop.setVirtualSplitter(SexSplitter())
        simu = simulator(pop )
        simu.evolve(initOps = [initSex()],
            matingScheme = randomMating(),
            postOps = [ KamMutator(k=2, rates=0.5, loci=[1,4], subPops=[(0, 0)])],
            gen = 1)
        pop = simu.extract(0)
        Stat(pop, alleleFreq=range(5))
        for loc in [0, 2, 3]:
            self.assertEqual(pop.dvars().alleleFreq[loc][0], 1.0)
        for loc in [1, 4]:
            self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 1.0)
        for ind in pop.individuals():
            for loc in [0, 2, 3]:
                self.assertEqual(ind.allele(loc, 0), 0)
                self.assertEqual(ind.allele(loc, 1), 0)
            if ind.sex() == FEMALE:
                self.assertEqual(ind.genotype(), [0]*10)

    def testUntouchedLoci(self):
        'Testing if mutator would mutate irrelevant locus'
        simu = simulator( population(size=1000, ploidy=2, loci=[2, 3]))
        simu.evolve(initOps = [initSex()],
            matingScheme = randomMating(),
            postOps = [ KamMutator(k=2, rates=0.5, loci=[1,4])], gen=200)
        self.assertGenotype(simu.population(0), 0,
            loci=[0,2,3])

    def testSnpMutator(self):
        'Testing diallelic mutator (SNP mutator)'
        simu = simulator( population(size=1000, ploidy=2, loci=[2, 3]), rep=5)
        simu.evolve(
                initOps = [ initSex(), initByFreq([.5, .5], loci=[0, 4])],
            matingScheme = randomMating(),
                postOps = [SnpMutator(u=0.1, loci=[0, 4]),
                    #stat(alleleFreq=[0, 4]),
                    #pyEval(r'"%.3f %.3f\n" % (alleleFreq[0][0], alleleFreq[4][0])')
                ],
                gen=100)
        self.assertGenotype(simu.population(0), 0,
            loci=[1, 2, 3])
        # fewer and fewer allele 0
        self.assertGenotypeFreq(simu.population(0),
            [0.], [0.05], loci=[0, 4])

    def testAlleleMapping(self):
        'Testing the allele mapping feature'
        simu = simulator(population(size=1000, ploidy=2, loci=[2, 3]),
            rep=5)
        simu.evolve(
                initOps = [initSex(), initByFreq([0, 0, 0, 0, 0, .5, .5], loci=[0, 4])],
            matingScheme = randomMating(),
                postOps = [SnpMutator(u=0.1, loci=[0, 4],
                    mapIn=[0, 0, 0, 0, 0, 0, 1],
                    mapOut=[5, 6]),
                    #stat(alleleFreq=[0, 4]),
                    #pyEval(r'"%.3f %.3f\n" % (alleleFreq[0][5], alleleFreq[4][5])')
                ],
                gen=100)
        self.assertGenotype(simu.population(0), 0,
            loci=[1, 2, 3])
        # fewer and fewer allele 0
        self.assertGenotypeFreq(simu.population(0),
            [0, 0, 0, 0, 0, 0, 0.95], [0, 0, 0, 0, 0, 0.05, 1], loci=[0, 4])
        #
        def mapIn(allele):
            return allele - 5
        def mapOut(allele):
            return allele + 5
        simu = simulator(population(size=1000, ploidy=2, loci=[2, 3]),
            rep=5)
        simu.evolve(
                initOps = [initSex(), initByFreq([0, 0, 0, 0, 0, .5, .5], loci=[0, 4])],
            matingScheme = randomMating(),
                postOps = [SnpMutator(u=0.1, loci=[0, 4],
                    mapIn=mapIn, mapOut=mapOut),
                    #stat(alleleFreq=[0, 4]),
                    #pyEval(r'"%.3f %.3f\n" % (alleleFreq[0][5], alleleFreq[4][5])')
                ],
                gen=100)
        self.assertGenotype(simu.population(0), 0,
            loci=[1, 2, 3])
        # fewer and fewer allele 0
        self.assertGenotypeFreq(simu.population(0),
            [0, 0, 0, 0, 0, 0, 0.95], [0, 0, 0, 0, 0, 0.05, 1], loci=[0, 4])



    def testKamMutator(self):
        'Testing k-allele mutator'
        simu = simulator( population(size=1000, ploidy=2, loci=[2, 3]),
            rep=5)
        # simu.apply( [ initSex(), initByFreq([.2,.8])])
        simu.evolve(
                initOps = [ initSex(), initByFreq([.2,.8])],
            matingScheme = randomMating(),
                postOps = [ KamMutator(k=2, rates=0.1)],
                gen=200)
        # at loci
        simu = simulator( population(size=10000, ploidy=2, loci=[2, 3]),
            rep=5)
        simu.evolve(
            initOps = [initSex()],
            matingScheme = randomMating(),
            postOps = [ KamMutator(k=2, rates=0.1, loci=[0,4])],
            gen = 1)
        # frequency seems to be OK.
        self.assertGenotypeFreq(simu.population(0),
            [0.85],[0.95], loci=[0,4])
        self.assertGenotype(simu.population(0), 0,
            loci=[1,2,3])

    def testSmmMutator(self):
        'Testing generalized step-wise mutation mutator'
        if moduleInfo()['alleleType'] == 'binary':
            return
        simu = simulator( population(size=1000, ploidy=2, loci=[2, 3]),
            rep=5)
        # simu.apply( [ initSex(), initByFreq([.2,.8])])
        simu.evolve(initOps=[initSex(), initByFreq([.2,.8])],
            matingScheme = randomMating(),
             postOps = [ SmmMutator(rates=0.2)], gen=200)
        # at loci
        simu = simulator( population(size=10000, ploidy=2, loci=[2, 3]),
            rep=5)
        simu.evolve(initOps = [initSex()],
            matingScheme = randomMating(),
            postOps = [ SmmMutator(rates=0.2, loci=[0,4])],
            gen = 1)
        # frequency seems to be OK.
        self.assertGenotypeFreq(simu.population(0),
            [0.85],[0.95], loci=[0,4])
        self.assertGenotype(simu.population(0), 0,
            loci=[1,2,3])


    def testPyMutator(self):
        'Testing the hybrid PyMutator'
        pop = population(size=10, loci=[2])
        # cutom mutator
        def mut(x):
            return 1
        PyMutate(pop, rates=1, func=mut)
        self.assertEqual(pop.individual(0).allele(0), 1)

    def testMixedMutator(self):
        'Testing mixed mutator'
        pop = population(1000, loci=[1])
        simu = simulator(pop)
        self.assertRaises(exceptions.ValueError,
            MixedMutator, mutators=[KamMutator(k=10), KamMutator(k=10)],
            prob=[0.2, 0.4])
        self.assertRaises(exceptions.ValueError,
            MixedMutator, mutators=[KamMutator(k=10), KamMutator(k=10)],
            prob=[0.2, 0.4, 0.4])

    def testContextMutator(self):
        'Testing context mutator'
        simu = simulator(population(50000, loci=[3, 3]))
        simu.evolve(
            # initialize locus by 0, 0, 0, 1, 0, 1
            initOps = [initSex(), initByValue([1, 1], loci=[3, 5])],
            matingScheme = randomMating(),
            postOps = [
                ContextMutator(mutators=[
                    SnpMutator(u=0.1),
                    SnpMutator(u=1),
                    ],
                    contexts=[(0, 0), (1, 1)],
                    loci=[1, 4],
                    rates=0.1
                ),
            ], 
            gen = 1
        )
        pop = simu.extract(0)
        Stat(pop, alleleFreq=(1, 4))
        # test combined allele frequency
        self.assertAlmostEqual(pop.dvars().alleleFreq[1][1], 0.01, places=2)
        self.assertAlmostEqual(pop.dvars().alleleFreq[4][1], 0.1, places=2)
        # test parameters
        self.assertRaises(exceptions.ValueError, ContextMutate, pop, 
            contexts=[(0, 0, 0)])
        self.assertRaises(exceptions.ValueError, ContextMutate, pop, 
            mutators=[SnpMutator(u=0.1)], contexts=[(0, 0), (1, 1)])
        self.assertRaises(exceptions.ValueError, ContextMutate, pop, 
            mutators=[SnpMutator(u=0.1), SnpMutator(u=0.01)],
            contexts=[(0, 0), (1, 1, 2, 2)])

    def testPointMutator(self):
        'Testing point mutator'
        # test point mutator
        pop = population(size=10, ploidy=2, loci=[5])
        InitByValue(pop, value=[[1]*5, [2]*5], proportions=[.3,.7])
        PointMutate(pop, inds=[1,2,3], allele=0, loci=[1,3])
        self.assertEqual(pop.individual(1).allele(1,0), 0)
        self.assertNotEqual(pop.individual(1).allele(1,1), 0)
        PointMutate(pop, inds=[1,2,3], ploidy=[1],
            allele=0, loci=[1,2])
        self.assertEqual(pop.individual(1).allele(2,1), 0)
        self.assertNotEqual(pop.individual(1).allele(2,0), 0)

    def testSnpMutator(self):
        'Testing SNP mutator'
        cnt0 = 0
        cnt1 = 0
        for i in range(50):
            pop = population(size=10000, loci=[1])
            InitByFreq(pop, [0.6, 0.4])
            SnpMutate(pop, u=0.2, v=0.1, loci=0)
            Stat(pop, alleleFreq=[0])
            # u = 10000*2*(0.6-0.12+0.04), v = 10000*2*(0.4-0.04+0.12)
            cnt0 += pop.dvars().alleleNum[0][0]
            cnt1 += pop.dvars().alleleNum[0][1]
        self.assertEqual( abs(cnt0/50. - 20000*(0.6-0.12+0.04)) < 20, True)
        self.assertEqual( abs(cnt1/50. - 20000*(0.4-0.04+0.12)) < 20, True)

    def testMutationRate(self):
        'Testing mutation rate'
        cnt = 0
        for i in range(50):
            pop = population(size=10000, loci=[1])
            # Mutate autosome
            SnpMutate(pop, u=0.01, loci=0)
            Stat(pop, alleleFreq=[0])
            # 10000 x 2 x 0.01 = 200
            cnt += pop.dvars().alleleNum[0][1]
        self.assertEqual( abs(cnt/50. - 200) < 5, True)

    def testMutationSexChromosomes(self):
        'Testing mutation on chromosome X'
        cnt = 0
        for i in range(50):
            pop = population(size=10000, loci=[1, 1], chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
            InitSex(pop, sex=[MALE, FEMALE])
            # Mutate X chromosomes
            SnpMutate(pop, u=0.01, loci=0)
            Stat(pop, alleleFreq=[0])
            # MALE: 5000 x 0.01, FEMALE: 5000 x 2 x 0.01 = 50 + 100
            cnt += pop.dvars().alleleNum[0][1]
        self.assertEqual( abs(cnt/50. - (5000*0.01*2 + 5000*0.01)) < 5, True)
        #
        # Mutate Y chromosomes
        cnt = 0
        for i in range(50):
            pop = population(size=10000, loci=[1, 1], chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
            InitSex(pop, sex=[MALE, FEMALE])
            SnpMutate(pop, u=0.01, loci=1)
            Stat(pop, alleleFreq=[1])
            # MALE: 5000 x 0.01 = 50
            cnt += pop.dvars().alleleNum[1][1]
        self.assertEqual( abs(cnt/50. - 5000*0.01) < 2, True)

if __name__ == '__main__':
    unittest.main()
