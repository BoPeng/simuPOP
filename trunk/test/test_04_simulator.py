#!/usr/bin/env python
#
# testing for simulator
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#
import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

class TestSimulator(unittest.TestCase):

    def testClone(self):
        'Testing simulator::clone() of simulator'
        pop = population(size=[100, 40], loci=[2, 5])
        simu = simulator(pop, rep = 3)
        simu.evolve(
            initOps = [InitSex(), InitByFreq([0.3, .7])],
            matingScheme=randomMating(),
            postOps = [stat(alleleFreq=range(pop.totNumLoci()))],
            gen = 10
        )
        simu1 = simu.clone()
        for i in range(3):
            self.assertEqual(simu.population(i), simu1.population(i))
        self.assertEqual(simu1.dvars(0).gen, simu.dvars(0).gen)
        # this test should be enough
        self.assertEqual(simu, simu1)
        # test if a cloned simulator can evolve again
        simu1.evolve(
            initOps = [InitSex()],
            matingScheme=randomMating(),
            postOps = stat(alleleFreq=range(pop.totNumLoci())),
            gen = 20
        )

    def testEvolve(self):
        'Testing simulator:: evolve()'
        pop = population(size=1, loci=[1])
        simu = simulator(pop, rep=3)
        # no terminator, no ending generation is specified
        self.assertRaises(exceptions.ValueError, simu.evolve)
        # sample
        pop = population(size=[200, 80], loci=[3])
        pop.evolve(
            initOps = [InitSex(), InitByFreq([0.2, 0.8])],
            matingScheme=randomMating(ops = Recombinator(rates=0.001)),
            postOps = stat(alleleFreq=[1]),
            finalOps = stat(),
            gen=10
        )

    def testCreateSimulator(self):
        'Testing the construction of simulator'
        pop = population(size=[20, 80], loci=1)
        pop1 = population(size=[20, 40], loci=2)
        simu = simulator([pop, pop1], steal=False)
        self.assertEqual(pop.popSize(), 100)
        self.assertEqual(pop1.popSize(), 60)
        self.assertEqual(simu.population(0).popSize(), 100)
        self.assertEqual(simu.population(1).popSize(), 60)
        # steal
        simu = simulator([pop, pop1])
        self.assertEqual(pop.popSize(), 0)
        self.assertEqual(pop1.popSize(), 0)
        self.assertEqual(simu.population(0).popSize(), 100)
        self.assertEqual(simu.population(1).popSize(), 60)
        # rep
        pop = population(size=[20, 80], loci=1)
        pop1 = population(size=[20, 40], loci=2)
        simu = simulator([pop, pop1], rep=3, steal=False)
        self.assertEqual(pop.popSize(), 100)
        self.assertEqual(pop1.popSize(), 60)
        self.assertEqual(simu.population(2).popSize(), 100)
        self.assertEqual(simu.population(3).popSize(), 60)
        self.assertEqual(simu.numRep(), 6)

    def testExtract(self):
        'Testing simulator::extract(rep), numRep()'
        pop = population(size=[20, 80], loci=[3])
        simu = simulator(pop, rep=5)
        repnum = simu.numRep()
        simu.extract(2)
        self.assertEqual(simu.numRep(), repnum-1)
        self.assertRaises(exceptions.IndexError, simu.population, 5)
        self.assertRaises(exceptions.IndexError, simu.extract, 5)

    def testAdd(self):
        'Testing simulator::add(pop)'
        pop = population(size=[20, 80], loci=[3])
        simu = simulator(pop, rep=5)
        repnum = simu.numRep()
        pop1 = population(size=[20, 50], loci=[2])
        simu.add(pop1)
        self.assertEqual(simu.numRep(), repnum+1)
        self.assertEqual(simu.population(repnum).subPopSizes(), (20, 50))
        simu.population(repnum).removeSubPops(1)
        self.assertEqual(simu.population(repnum).subPopSizes(), (20,))
        self.assertEqual(pop1.subPopSizes(), (0,))
        # add without strealing
        pop = population(size=[300, 500], loci=1)
        repnum = simu.numRep()
        simu.add(pop, steal=False)
        self.assertEqual(simu.numRep(), repnum + 1)
        self.assertEqual(pop.subPopSizes(), (300, 500))



    def testPopulation(self):
        'Testing simulator::population(rep), populations()'
        pop = population(size=1000, loci=[1])
        simu = simulator(pop, rep=5)
        self.assertEqual(pop.popSize(), 0)
        # pop is not affected if simu changes
        for rep in range(5):
            for idx in range(simu.population(rep).popSize()):
                simu.population(rep).individual(idx).setAllele(1, 0)
                self.assertEqual(simu.population(rep).individual(idx).allele(0), 1)
        # reference to the rep-th population
        for rep in range(5):
            pop = simu.population(rep)
            for idx in range(pop.popSize()):
                pop.individual(rep).setAllele(1, 0)
                self.assertEqual(simu.population(rep).individual(idx).allele(0), 1)
        # independent copy of the population
        for rep in range(5):
            pop1 = simu.population(rep).clone()
            for idx in range(pop.popSize()):
                simu.population(rep).individual(idx).setAllele(0, 0)
                self.assertEqual(pop1.individual(idx).allele(0), 1)
                self.assertEqual(simu.population(rep).individual(idx).allele(0), 0)
        # populations
        for onepop in simu.populations():
            for idx in range(onepop.popSize()):
                self.assertEqual(onepop.individual(idx).allele(0), 0)
        self.assertRaises(exceptions.IndexError, simu.population, 5)

    def testAddInfoField(self):
        'Testing setMatingScheme(matingScheme)'
        simu = simulator(population(100, infoFields=['a']), rep=3)
        simu.evolve(initOps=[InitSex()], 
            matingScheme=CloneMating(), gen=1)
        simu.evolve(initOps=[InitSex()],
            matingScheme=randomMating(), gen=1)

    def testVars(self):
        'Testing simulator::vars(rep), vars(rep, subPop), dvars(rep), dvars(rep, subPop)'
        pop = population(size=100, loci=[2, 4])
        initByFreq(pop, [.2, .3, .5])
        Stat(pop, alleleFreq=range(0, 6))
        simu = simulator(pop, rep=5)
        for rep in range(5):
            self.assertEqual(len(simu.vars(rep)["alleleFreq"]), 6)
            self.assertEqual(len(simu.dvars(rep).alleleFreq), 6)
        # with subPop
        pop = population(size=[20, 80], loci=[2, 4])
        initByFreq(pop, [.2, .3, .5])
        Stat(pop, alleleFreq=range(0, 6), vars=['alleleFreq', 'alleleFreq_sp'])
        simu = simulator(pop, rep=5)
        for rep in range(5):
            self.assertEqual(len(simu.vars(rep)["alleleFreq"]), 6)
            self.assertEqual(len(simu.dvars(rep, 1).alleleFreq), 6)

if __name__ == '__main__':
    unittest.main()
