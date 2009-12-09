#!/usr/bin/env python
#
# unittests for mating schemes
#
# Author:
#   Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, random, math, sets, exceptions

class TestMatingSchemes(unittest.TestCase):

    def getFamSize(self, ms, gen=1, N=1000):
        '''Check the number of offspring for each family using
           information field father_idx'''
        simu = simulator(
            population(size=[N], infoFields=['father_idx', 'mother_idx']),
            matingScheme=ms)
        simu.evolve(
            initOps = initSex(),
            duringOps = parentsTagger(),
            gen=gen)
        # get the parents of each offspring
        parents = [(x, y) for x, y in zip(simu.population(0).indInfo('mother_idx'),
            simu.population(0).indInfo('father_idx'))]
        # Individuals with identical parents are considered as siblings.
        famSize = []
        lastParent = (-1, -1)
        for parent in parents:
            if parent == lastParent:
                famSize[-1] += 1
            else:
                lastParent = parent
                famSize.append(1)
        return famSize

    def testSubPopSizes(self):
        'Testing parameter subPopSize of mating schemes'
        def demo(gen, oldSize=[]):
            return (500 + gen*10, 1000 + gen*10)
        def demoSize(pop):
            gen = pop.dvars().gen
            intended_size = demo(gen)
            self.assertEqual(pop.subPopSizes(), intended_size)

        pop = population(size=[500, 1000], infoFields=['migrate_to'])
        simu = simulator(pop, randomMating(subPopSize=demo))
        simu.evolve(
            initOps = [initSex()],
            postOps = pyOperator(func=demoSize),
            gen = 100
        )


    def testNumOffspring(self):
        'Testing number of offspring'
        # single number
        self.assertEqual(
             self.getFamSize(randomMating(numOffspring=5), N=50),
             [5]*10)

        # Python function
        def nos(gen):
            return gen%2+1
        self.assertEqual(
            self.getFamSize(randomMating(numOffspring=nos), gen=2),
            [2]*500)
        self.assertEqual(
            self.getFamSize(randomMating(numOffspring=nos), gen=3),
            [1]*1000)
        # randomnumber
        def nos(gen):
            return random.randint(1, 3)
        cnt = self.getFamSize(randomMating(numOffspring=nos), N=1000)
        self.assertEqual(sum(cnt), 1000)
        num = [ cnt.count(i) for i in range(1, 4) ]
        #
        mean = sum(num)/3.
        for i in range(3):
            assert num[i] < mean + 50 and num[i] > mean - 50
        # GeometricDistribution
        p = 0.33
        cnt = self.getFamSize(randomMating(
            numOffspring=(GeometricDistribution, p)), N=10000)
        # mean should be 1/p, variance (1-p)/(p*p)
        mean = sum(cnt)*1.0/len(cnt)
        var = sum([x*x for x in cnt])*1.0/len(cnt) - mean*mean
        self.assertEqual(abs(mean - 1./p) < 0.1, True)
        self.assertEqual(abs(var - (1-p)/(p*p)) < 1, True)
        # PoissonDistribution
        p = 3
        cnt = self.getFamSize(randomMating(
            numOffspring=(PoissonDistribution, p)), N=100000)
        mean = sum(cnt)*1.0/len(cnt)
        var = sum([x*x for x in cnt])*1.0/len(cnt) - mean*mean
        self.assertEqual(abs(mean - (p+1)) < 0.1, True)
        self.assertEqual(abs(var - p) < 0.2, True)
        # BinomialDistribution
        p = 0.3
        n = 10
        cnt = self.getFamSize(randomMating(
            numOffspring=(BinomialDistribution, p, n)), N=10000)
        mean = sum(cnt)*1.0/len(cnt)
        var = sum([x*x for x in cnt])*1.0/len(cnt) - mean*mean
        self.assertEqual(abs(mean - ((n-1)*p+1)) < 0.1, True)
        self.assertEqual(abs(var - (n-1)*p*(1-p)) < 0.2, True)
        # UniformDistribution
        a = 3
        b = 6
        cnt = self.getFamSize(randomMating(
            numOffspring=(UniformDistribution, a, b)), N=10000)
        mean = sum(cnt)*1.0/len(cnt)
        var = sum([x*x for x in cnt])*1.0/len(cnt) - mean*mean
        self.assertEqual(abs(mean - (a + b)/2.) < 0.1, True)

    def checkSexMode(self, ms):
        simu = simulator(
            population(size=[40]),
            matingScheme=ms)
        simu.evolve(initOps = initSex(), gen=1)
        # return individual sex as a string
        return ''.join([ind.sexChar() for ind in simu.population(0).individuals()])

    def testSexMode(self):
        'Testing parameter sexMode of mating schemes'
        # noSex
        self.assertEqual(
            self.checkSexMode(randomMating(sexMode=NoSex)),
            'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM')
        # NumOfMales
        self.assertEqual(
            self.checkSexMode(randomMating(numOffspring=3,
            sexMode=(NumOfMales, 1))),
            'MFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFM')
        # NumOfFemales
        self.assertEqual(
            self.checkSexMode(randomMating(numOffspring=4,
            sexMode=(NumOfFemales, 2))),
            'FFMMFFMMFFMMFFMMFFMMFFMMFFMMFFMMFFMMFFMM')
        # ProbOfMales
        pop = population(10000)
        simu = simulator(pop, randomMating(sexMode=(ProbOfMales, 0.3)))
        simu.evolve(
            initOps = [initSex(), initByFreq([0.5, 0.5])],
            postOps = [
                stat(numOfMales=True),
                # number of male should be variable, but not too much
                terminateIf('numOfMales < 2500 or numOfMales > 3500'),
            ],
            gen = 10
        )
        self.assertEqual(simu.gen(), 10)

    def testMonoMating(self):
        'Testing monogemous mating scheme'
        pop = population(size=[2000], loci=[3,5], infoFields=['father_idx', 'mother_idx'])
        InitByFreq(pop, [0.2, 0.3, 0.5])
        simu = simulator(pop, monogamousMating(numOffspring=2, sexMode=(NumOfMales, 1)))
        simu.evolve(
            initOps = initSex(sex=(Male, Female)), 
            duringOps = parentsTagger(),
            gen = 5)
        self.assertEqual(len(sets.Set(simu.population(0).indInfo('father_idx'))), 1000)
        self.assertEqual(len(sets.Set(simu.population(0).indInfo('mother_idx'))), 1000)
        pop = simu.extract(0)
        self.assertEqual([ind.sex() for ind in pop.individuals()], [1,2]*1000)
               
    def testHeteroMating(self):
        'Testing heterogeneous mating schemes'
        pop = population(size=[10000, 10000], loci=[2], infoFields=['father_idx', 'mother_idx'])
        simu = simulator(pop,
            heteroMating(
                [randomMating(numOffspring=2, subPops=0),
                randomMating(numOffspring=4, subPops=1)])
        )
        simu.evolve(
            initOps = initSex(),
            duringOps = parentsTagger(),
            gen=10)      
        parents = [(x, y) for x, y in zip(simu.population(0).indInfo('mother_idx'),
            simu.population(0).indInfo('father_idx'))]
        # Individuals with identical parents are considered as siblings.
        famSize = []
        lastParent = (-1, -1)
        for parent in parents:
            if parent == lastParent:
                famSize[-1] += 1
            else:
                lastParent = parent
                famSize.append(1)
        self.assertEqual(famSize, [2]*5000+[4]*2500)

        # virtual subpopulation
        pop = population(size =[20000, 20000], loci=[2], infoFields=['father_idx', 'mother_idx'])
        pop.setVirtualSplitter(proportionSplitter([0.2, 0.8]))
        simu = simulator(pop, heteroMating(
            matingSchemes = [
                randomMating(numOffspring=1, subPops=[(0,0)]),
                randomMating(numOffspring=2, subPops=[(1,1)])
            ])
        )
        simu.evolve(
            initOps = initSex(),
            duringOps = parentsTagger(),
            gen =10
        )
        parents = [(x, y) for x, y in zip(simu.population(0).indInfo('mother_idx'),
            simu.population(0).indInfo('father_idx'))]
        # Individuals with identical parents are considered as siblings.
        famSize = []
        lastParent = (-1, -1)
        for parent in parents:
            if parent == lastParent:
                famSize[-1] += 1
            else:
                lastParent = parent
                famSize.append(1)
        self.assertEqual(famSize, [1]*20000+[2]*10000)
         
    def testPolygamousMating(self):
        'Testing polygamous mating scheme'
        pop = population(size=[200], loci=[3,5], infoFields=['father_idx', 'mother_idx'])
        InitByFreq(pop, [0.2, 0.3, 0.5])
        # exactly 100 males and 100 females
        for i in range(100):
            pop.individual(i).setSex(Male)
            pop.individual(100+i).setSex(Female)
        simu = simulator(pop, polygamousMating(polySex=Male, polyNum=3, numOffspring=2))
        simu.evolve(
            initOps = [],
            duringOps = parentsTagger(),
            gen = 1)
        # there is only one Male...
        fi = simu.population(0).indInfo('father_idx')
        self.assertEqual(fi[0], fi[1])
        self.assertEqual(fi[0], fi[5])
        self.assertNotEqual(fi[0], fi[6])
        mi = simu.population(0).indInfo('mother_idx')
        self.assertEqual(mi[0], mi[1])
        self.assertNotEqual(mi[0], mi[2])              
       
    def testPedigreeMating(self):
        'Testing pedigree mating using a population object'
        pop = population(size=[100, 100], loci=[2, 5], ancGen=-1,
            infoFields=['father_idx', 'mother_idx'])
        simu = simulator(pop, randomMating())
        simu.evolve(
            initOps = initSex(),
            duringOps = parentsTagger(),
            gen = 20
        )
        return
        ped = pedigree(simu.extract(0), infoFields=['father_idx', 'mother_idx'])
        simu = simulator(pop, pedigreeMating(ped,
            offspringGenerator(mendelianGenoTransmitter())))
        simu.evolve(
            duringOps = parentsTagger(),
            gen = 100
        )
        ped1 = simu.extract(0)
        self.assertEqual(simu.gen(), 21)
        # compare two populations!
        self.assertEqual(ped.ancestralGens(), ped1.ancestralGens())
        self.assertEqual(ped.ancestralGens(), 20)
        for gen in range(21):
            ped.useAncestralGen(gen)
            ped1.useAncestralGen(gen)
            self.assertEqual(ped.indInfo('father_idx'), ped1.indInfo('father_idx'))
            self.assertEqual(ped.indInfo('mother_idx'), ped1.indInfo('mother_idx'))

    def testSequentialParentsChooser(self):
        'Testing sequential parent chooser'
        pop = population(size=[100, 200], infoFields=['parent_idx'])
        InitByFreq(pop, [.3, .7])
        simu = simulator(pop, homoMating(
            sequentialParentsChooser(),
            offspringGenerator(selfingGenoTransmitter())))
        simu.evolve(
            duringOps = parentsTagger(infoFields='parent_idx'), 
            gen=1)

    def testRandomParentsChooser(self):
        'Testing sequential parent chooser'
        def traj(gen):
            return [0.5 + gen*0.01]
        pop = population(size=[1000, 2000], infoFields=['parent_idx'])
        InitByFreq(pop, [.2, .8])
        simu = simulator(pop, homoMating(
            randomParentChooser(),
            offspringGenerator(selfingGenoTransmitter())))
        simu.evolve(
            duringOps = parentsTagger(infoFields='parent_idx'), 
            gen=1)

    def testPyParentsChooserRetValue(self):
        'Testing the return value of Python parents chooser'
        import random
        def retIndex(pop, sp):
            while True:
                yield random.randint(0, pop.subPopSize(sp) - 1)
        def retIndexes(pop, sp):
            while True:
                yield random.randint(0, pop.subPopSize(sp) - 1), random.randint(0, pop.subPopSize(sp) -1)
        def retInd(pop, sp):
            while True:
                yield pop.individual(random.randint(0, pop.subPopSize(sp) - 1))
        def retInds(pop, sp):
            while True:
                yield pop.individual(random.randint(0, pop.subPopSize(sp) - 1)), \
                     pop.individual(random.randint(0, pop.subPopSize(sp) - 1))
        def retPop(pop, sp):
            while True:
                yield pop
        def retWrongIndex(pop, sp):
            while True:
                yield pop.subPopSize(sp)
        def retWrongIndexes(pop, sp):
            while True:
                yield 0, pop.subPopSize(sp)
        def testPyRetValue(func):
            pop = population([200]*5)
            simu = simulator(pop,
                homoMating(
                    pyParentsChooser(func),
                    offspringGenerator(cloneGenoTransmitter())
                )
            )
            simu.evolve(
                gen = 5
            )
        testPyRetValue(retIndex)
        testPyRetValue(retIndexes)
        testPyRetValue(retInd)
        testPyRetValue(retInds)
        self.assertRaises(exceptions.ValueError, testPyRetValue, retPop)
        self.assertRaises(exceptions.ValueError, testPyRetValue, retWrongIndex)
        self.assertRaises(exceptions.ValueError, testPyRetValue, retWrongIndexes)

  
    def testHaploidRandomMating(self):
        'Testing random mating in haploid populations'
        pop = population(size=[50, 100], loci=[5]*5, ploidy=1,
            chromTypes=[Customized]*5)
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, 
            randomMating(ops=[mitochondrialGenoTransmitter()]))
        simu.evolve(
            initOps = [initSex(),
                # female has [1]
                initByValue([1]*25, subPops=[(0, 1), (1, 1)]),
                ],
            gen = 1
        )
        self.assertEqual(simu.population(0).genotype(), [1]*(150*25))



if __name__ == '__main__':
  unittest.main()
  sys.exit(0)


