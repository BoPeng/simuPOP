#!/usr/bin/env python
#
# testing for Simulator
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import unittest, os, sys
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

class TestSimulator(unittest.TestCase):

    def testClone(self):
        'Testing Simulator::clone() of Simulator'
        pop = Population(size=[100, 40], loci=[2, 5])
        simu = Simulator(pop, rep = 3)
        simu.evolve(
            initOps = [InitSex(), InitGenotype(freq=[0.3, .7])],
            matingScheme=RandomMating(),
            postOps = [Stat(alleleFreq=list(range(pop.totNumLoci())))],
            gen = 10
        )
        simu1 = simu.clone()
        for i in range(3):
            self.assertEqual(simu.population(i), simu1.population(i))
        self.assertEqual(simu1.dvars(0).gen, simu.dvars(0).gen)
        # this test should be enough
        self.assertEqual(simu, simu1)
        # test if a cloned Simulator can evolve again
        simu1.evolve(
            initOps = [InitSex()],
            matingScheme=RandomMating(),
            postOps = Stat(alleleFreq=list(range(pop.totNumLoci()))),
            gen = 20
        )

    def testEvolve(self):
        'Testing Simulator:: evolve()'
        pop = Population(size=1, loci=[1])
        simu = Simulator(pop, rep=3)
        # no terminator, no ending generation is specified
        # but a mating scheme can also terminate evolution
        #self.assertRaises(ValueError, simu.evolve)
        # sample
        pop = Population(size=[200, 80], loci=[3])
        pop.evolve(
            initOps = [InitSex(), InitGenotype(freq=[0.2, 0.8])],
            matingScheme=RandomMating(ops = Recombinator(rates=0.001)),
            postOps = Stat(alleleFreq=[1]),
            finalOps = Stat(),
            gen=10
        )

    def testCreateSimulator(self):
        'Testing the construction of Simulator'
        pop = Population(size=[20, 80], loci=1)
        pop1 = Population(size=[20, 40], loci=2)
        simu = Simulator([pop, pop1], stealPops=False)
        self.assertEqual(pop.popSize(), 100)
        self.assertEqual(pop1.popSize(), 60)
        self.assertEqual(simu.population(0).popSize(), 100)
        self.assertEqual(simu.population(1).popSize(), 60)
        # steal
        simu = Simulator([pop, pop1])
        self.assertEqual(pop.popSize(), 0)
        self.assertEqual(pop1.popSize(), 0)
        self.assertEqual(simu.population(0).popSize(), 100)
        self.assertEqual(simu.population(1).popSize(), 60)
        # rep
        pop = Population(size=[20, 80], loci=1)
        pop1 = Population(size=[20, 40], loci=2)
        simu = Simulator([pop, pop1], rep=3, stealPops=False)
        self.assertEqual(pop.popSize(), 100)
        self.assertEqual(pop1.popSize(), 60)
        self.assertEqual(simu.population(2).popSize(), 100)
        self.assertEqual(simu.population(3).popSize(), 60)
        self.assertEqual(simu.numRep(), 6)

    def testExtract(self):
        'Testing Simulator::extract(rep), numRep()'
        pop = Population(size=[20, 80], loci=[3])
        simu = Simulator(pop, rep=5)
        repnum = simu.numRep()
        simu.extract(2)
        self.assertEqual(simu.numRep(), repnum-1)
        self.assertRaises(IndexError, simu.population, 5)
        self.assertRaises(IndexError, simu.extract, 5)

    def testAdd(self):
        'Testing Simulator::add(pop)'
        pop = Population(size=[20, 80], loci=[3])
        simu = Simulator(pop, rep=5)
        repnum = simu.numRep()
        pop1 = Population(size=[20, 50], loci=[2])
        simu.add(pop1)
        self.assertEqual(simu.numRep(), repnum+1)
        self.assertEqual(simu.population(repnum).subPopSizes(), (20, 50))
        simu.population(repnum).removeSubPops(1)
        self.assertEqual(simu.population(repnum).subPopSizes(), (20,))
        self.assertEqual(pop1.subPopSizes(), (0,))
        # add without strealing
        pop = Population(size=[300, 500], loci=1)
        repnum = simu.numRep()
        simu.add(pop, stealPop=False)
        self.assertEqual(simu.numRep(), repnum + 1)
        self.assertEqual(pop.subPopSizes(), (300, 500))



    def testPopulation(self):
        'Testing Simulator::Population(rep), populations()'
        pop = Population(size=1000, loci=[1])
        simu = Simulator(pop, rep=5)
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
        self.assertRaises(IndexError, simu.population, 5)

    def testAddInfoField(self):
        'Testing setMatingScheme(matingScheme)'
        simu = Simulator(Population(100, infoFields=['a']), rep=3)
        simu.evolve(initOps=[InitSex()], 
            matingScheme=CloneMating(), gen=1)
        simu.evolve(initOps=[InitSex()],
            matingScheme=RandomMating(), gen=1)

    def testVars(self):
        'Testing Simulator::vars(rep), vars(rep, subPop), dvars(rep), dvars(rep, subPop)'
        pop = Population(size=100, loci=[2, 4])
        initGenotype(pop, freq=[.2, .3, .5])
        stat(pop, alleleFreq=list(range(0, 6)))
        simu = Simulator(pop, rep=5)
        for rep in range(5):
            self.assertEqual(len(simu.vars(rep)["alleleFreq"]), 6)
            self.assertEqual(len(simu.dvars(rep).alleleFreq), 6)
        # with subPop
        pop = Population(size=[20, 80], loci=[2, 4])
        initGenotype(pop, freq=[.2, .3, .5])
        stat(pop, alleleFreq=list(range(0, 6)), vars=['alleleFreq', 'alleleFreq_sp'])
        simu = Simulator(pop, rep=5)
        for rep in range(5):
            self.assertEqual(len(simu.vars(rep)["alleleFreq"]), 6)
            self.assertEqual(len(simu.dvars(rep, 1).alleleFreq), 6)

    def demo(self, pop):
        if pop.dvars().gen == 5:
            return [-1]
        else:
            return pop.subPopSizes()

    def testNegSize(self):
        'Testing negative population size returned by demographic function'
        pop = Population(size=[100], loci=[2])
        self.assertRaises(ValueError, pop.evolve,
            initOps=InitSex(),
            matingScheme=RandomMating(subPopSize=self.demo),
            gen=10)

if __name__ == '__main__':
    unittest.main()
