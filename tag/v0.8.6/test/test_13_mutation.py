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
from simuUtil import getGenotype
import unittest, os, sys, exceptions

class TestMutator(unittest.TestCase):

    def assertGenotype(self, pop, genotype,
        loci=[], subPop=[], indRange=[], atPloidy=[]):
        'Assert if the genotype of subPop of pop is genotype '
        geno = getGenotype(pop, loci, subPop, indRange, atPloidy)
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
        loci=[], subPop=[], indRange=[], atPloidy=[]):
        'Assert if the genotype has the correct allele frequency'
        geno = getGenotype(pop, loci, subPop, indRange, atPloidy)
        if AlleleType() == 'binary':
            if len(freqLow) == 1:    # only one
                freq0 = geno.count(0)*1.0 / len(geno)
                assert freq0 >= freqLow[0] and freq0 <= freqHigh[0]
            else:
                f0 = [freqLow[0], sum(freqLow[1:])]
                f1 = [freqHigh[0], sum(freqHigh[1:])]
                freq0 = geno.count(0)*1.0 / len(geno)
                freq1 = geno.count(1)*1.0 / len(geno)
                #print f0,f1,freq0,freq1
                assert freq0 >= f0[0] and freq0 <= f1[0]
                assert freq1 >= f0[1] and freq1 <= f1[1]
        else:
            for i in range(len(freqLow)):
                freq = geno.count(i)*1.0 / len(geno)
                assert freq >= freqLow[i] and freq <= freqHigh[i]

    def testMemberFunctions(self):
        'Testing common mutator functions'
        # need a parameter
        self.assertRaises(exceptions.ValueError,
            kamMutator)
        #
        m = kamMutator(rate=0.001)
        self.assertEqual(m.rate(), (0.001,))
        #
        m.setRate(0.1)
        self.assertEqual(m.rate(), (0.1,))
        #
        m.setRate(0.1, loci=[0,2,4])
        self.assertEqual(m.rate(), (0.1,))
        #
        m.setRate([0.1, 0.001], loci=[2,4])
        self.assertEqual(m.rate(), (0.1,0.001))
        #
        self.assertRaises(exceptions.ValueError,
            m.setRate, [0.1, 0.001], loci=[2,4,6])

    def testUntouchedLoci(self):
        'Testing if mutator would mutate irrelevant locus'
        simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
            randomMating() )
        simu.evolve([ kamMutator(rate=0.5, loci=[1,4])], gen=200)
        self.assertGenotype(simu.population(0), 0,
            loci=[0,2,3])

    def testKamMutator(self):
        'Testing k-allele mutator'
        simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
            randomMating(), rep=5)
        # simu.apply( [ initByFreq([.2,.8])])
        simu.evolve(
                preOps =    [ initByFreq([.2,.8])],
                ops = [ kamMutator(rate=0.1)], gen=200)
        # at loci
        simu = simulator( population(size=10000, ploidy=2, loci=[2, 3]),
            randomMating(), rep=5)
        simu.step([ kamMutator(rate=0.1, loci=[0,4])])
        # frequency seems to be OK.
        self.assertGenotypeFreq(simu.population(0),
            [0.85],[0.95], loci=[0,4])
        self.assertGenotype(simu.population(0), 0,
            loci=[1,2,3])

    def testSmmMutator(self):
        'Testing generalized step-wise mutation mutator'
        if AlleleType() == 'binary':
            return
        simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
            randomMating(), rep=5)
        # simu.apply( [ initByFreq([.2,.8])])
        simu.evolve(preOps=[initByFreq([.2,.8])],
             ops = [ smmMutator(rate=0.2)], gen=200)
        # at loci
        simu = simulator( population(size=10000, ploidy=2, loci=[2, 3]),
            randomMating(), rep=5)
        simu.step([ smmMutator(rate=0.2, loci=[0,4])])
        # frequency seems to be OK.
        self.assertGenotypeFreq(simu.population(0),
            [0.85],[0.95], loci=[0,4])
        self.assertGenotype(simu.population(0), 0,
            loci=[1,2,3])

    def testGsmMutator(self):
        'Testing generalized step-wise mutation mutator (imcomplete)'
        if AlleleType() == 'binary':
            return
        simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
            randomMating(), rep=5)
        simu.evolve( preOps = [ initByFreq([.2,.8])],
                ops = [ gsmMutator(rate=0.2)], gen=200)
        # at loci
        simu = simulator( population(size=10000, ploidy=2, loci=[2, 3]),
            randomMating(), rep=5)
        simu.population(0).arrGenotype(True)[:] = 1
        simu.step([ gsmMutator(rate=0.2, loci=[0,4])])
        # frequency? Genometric distribution of step
        #self.assertGenotypeFreq(simu.population(0),
        #    [0.85],[0.95], loci=[0,4])
        self.assertGenotype(simu.population(0), 1,
            loci=[1,2,3])

    def testPyMutator(self):
        pop = population(size=10, loci=[2])
        # cutom mutator
        def mut(x):
            return 1
        m = pyMutator(rate=1, func=mut)
        m.apply(pop)
        assert pop.individual(0).allele(0) == 1, \
            "PyMutator failed"

    def testMutationCount(self):
        N = 10000
        r = [0.001, 0.002, 0.003]
        G = 100
        pop = population(size=N, ploidy=2, loci=[5])
        simu = simulator(pop, randomMating())
        mut = kamMutator(rate=r, maxAllele=10, loci=[0,2,4])
        simu.evolve( [mut], gen=G)
        assert abs( mut.mutationCount(0) - 2*N*r[0]*G) < 200, \
            "Number of mutation event is not as expected."
        assert abs( mut.mutationCount(2) - 2*N*r[1]*G) < 200, \
            "Number of mutation event is not as expected."
        assert abs( mut.mutationCount(4) - 2*N*r[2]*G) < 200, \
            "Number of mutation event is not as expected."
        self.assertEqual( mut.mutationCount(1), 0)
        self.assertEqual( mut.mutationCount(3), 0)

    def testPointMutator(self):
        # test point mutator
        pop = population(size=10, ploidy=2, loci=[5])
        InitByValue(pop, value=[[1]*5, [2]*5], proportions=[.3,.7])
        PointMutate(pop, inds=[1,2,3], toAllele=0, loci=[1,3])
        assert pop.individual(1).allele(1,0) == 0
        assert pop.individual(1).allele(1,1) != 0
        PointMutate(pop, inds=[1,2,3], atPloidy=[1],
            toAllele=0, loci=[1,2])
        assert pop.individual(1).allele(2,1) == 0
        assert pop.individual(1).allele(2,0) != 0

if __name__ == '__main__':
    unittest.main()
