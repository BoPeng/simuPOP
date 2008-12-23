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
import unittest, os, sys

class TestSimulator(unittest.TestCase):

    def testClone(self):
        'Testing simulator::clone() of simulator'
        pop = population(size=[100, 40], loci=[2, 5])
        simu = simulator(pop, randomMating(), rep = 3)
        simu.evolve(
            preOps = [initByFreq([0.3, .7])],
            ops = [stat(alleleFreq=range(pop.totNumLoci()))],
            gen = 10
        )
        simu1 = simu.clone()
        for i in range(3):
            self.assertEqual(simu.population(i), simu1.population(i))
        self.assertEqual(simu1.gen(), simu.gen())
        # test if a cloned simulator can evolve again
        simu1.evolve(
            preOps = [initSex()],
            ops = [stat(alleleFreq=range(pop.totNumLoci()))],
            gen = 20
        )

    def testSave(self):
        'Testing simulator::save(filename)'
        pop = population(size=[100, 40], loci=[2, 5])
        simu = simulator(pop, randomMating(), rep = 3)
        simu.evolve(
            preOps = [initByFreq([0.3, .7])],
            ops = [stat(alleleFreq=range(pop.totNumLoci()))],
            gen = 10
        )

        simu.save("simuout");
        simu1 = LoadSimulator("simuout", randomMating())
        self.assertEqual(simu, simu1)

    def testSetGen(self):
        'Testing simulator::gen(), setGen(gen)'
        pop = population(size=1, loci=[1])
        simu = simulator(pop, randomMating(), rep=3)
        self.assertEqual(simu.numRep(), 3)
        self.assertEqual(simu.gen(), 0)
        simu.setGen(10)
        self.assertEqual(simu.gen(), 10)

    def testEvolve(self):
        'Testing simulator:: evolve()'
        pop = population(size=1, loci=[1])
        simu = simulator(pop, randomMating(), rep=3)
        self.assertEqual(simu.gen(), 0)
        # no terminator, no ending generation is specified
        self.assertRaises(exceptions.ValueError,
            simu.evolve, ops=[] )
        # sample
        simu = simulator(
            population(size=[20, 80], loci=[3]), randomMating())
        simu.evolve(
            preOps = [initByFreq([0.2, 0.8])],
            ops = [
                recombinator(rate=0.001),
                stat(alleleFreq=[1]),
            ],
            postOps = [stat()],
            #dryrun=True,
            gen=10
        )

    def testExtract(self):
        'Testing simulator::extract(rep), numRep()'
        pop = population(size=[20, 80], loci=[3])
        simu = simulator(pop, randomMating(), rep=5)
        repnum = simu.numRep()
        simu.extract(2)
        self.assertEqual(simu.numRep(), repnum-1)
        self.assertRaises(exceptions.IndexError, simu.extract, 5)

    def testPopulation(self):
        'Testing simulator::population(rep), populations()'
        pop = population(size=1000, loci=[1])
        simu = simulator(pop, randomMating(), rep=5)
        # pop is not affected if simu changes
        for rep in range(5):
            for idx in range(pop.popSize()):
                simu.population(rep).individual(idx).setAllele(1, 0)
                self.assertEqual(simu.population(rep).individual(idx).allele(0), 1)
                self.assertEqual(pop.individual(idx).allele(0), 0)
        # reference to the rep-th population
        for rep in range(5):
            pop = simu.population(rep)
            for idx in range(pop.popSize()):
                pop.individual(rep).setAllele(1, 0)
                self.assertEqual(simu.population(rep).individual(idx).allele(0), 1)
                self.assertEqual(pop.individual(idx).allele(0), 1)
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
        'Testing simulator::addInfoFields(field, init=0), addInfoFields(fields, init=0)'
        'setAncestralDepth(depth), setMatingScheme(matingScheme)'
        simu = simulator(population(100, infoFields=['a']), cloneMating(), rep=3)
        simu.setMatingScheme(randomMating())
        simu.addInfoField('b')
        self.assertEqual(simu.infoFields(), ('a', 'b'))
        simu.addInfoFields(['c', 'd'])
        simu.setAncestralDepth(2)
        simu.evolve(preOps = [initSex()], ops=[], gen=10)
        self.assertEqual(simu.population(0).ancestralGens(), 2)
        self.assertEqual(simu.infoFields(), ('a', 'b', 'c', 'd'))
        self.assertEqual(simu.population(0).infoFields(), ('a', 'b', 'c', 'd'))
        simu.population(0).addInfoField('l')
        self.assertRaises(exceptions.ValueError, simu.addInfoField, 'j')
        self.assertRaises(exceptions.ValueError, simu.evolve, ops=[])

    def testVars(self):
        'Testing simulator::vars(rep), vars(rep, subPop), dvars(rep), dvars(rep, subPop)'
        pop = population(size=100, loci=[2, 4])
        InitByFreq(pop, [.2, .3, .5])
        Stat(pop, alleleFreq=range(0, 6))
        simu = simulator(pop, randomMating(), rep=5)
        for rep in range(5):
            self.assertEqual(len(simu.population(rep).vars()["alleleFreq"]), 6)
            self.assertEqual(len(simu.population(rep).dvars().alleleFreq), 6)
        # with subPop
        pop = population(size=[20, 80], loci=[2, 4])
        InitByFreq(pop, [.2, .3, .5])
        Stat(pop, alleleFreq=range(0, 6))
        simu = simulator(pop, randomMating(), rep=5)
        for rep in range(5):
            self.assertEqual(len(simu.population(rep).vars(0)["alleleFreq"]), 6)
            self.assertEqual(len(simu.population(rep).dvars(1).alleleFreq), 6)

    def testIntegrity(self):
        'Testing checking of simulator integrity'
        simu = simulator(population(1), cloneMating())
        simu.population(0).addInfoField('something')
        # one can not change the genotype structure of the populations
        # in a simulator
        self.assertRaises(exceptions.ValueError, simu.evolve, ops=[])

if __name__ == '__main__':
    unittest.main()
