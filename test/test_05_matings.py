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

import unittest, os, sys, random, math
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

class TestMatingSchemes(unittest.TestCase):

    def getFamSize(self, numOffspring, gen=1, N=1000):
        '''Check the number of offspring for each family using
           information field father_idx'''
        simu = Simulator(
            Population(size=[N], infoFields=['father_idx', 'mother_idx']))
        simu.evolve(
            initOps = InitSex(),
            matingScheme=RandomMating(numOffspring=numOffspring,
                ops=[MendelianGenoTransmitter(), ParentsTagger()]),
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
        def demo(gen):
            return (500 + gen*10, 1000 + gen*10)
        def demoSize(pop):
            gen = pop.dvars().gen
            intended_size = demo(gen)
            self.assertEqual(pop.subPopSizes(), intended_size)
            return True

        pop = Population(size=[500, 1000], infoFields=['migrate_to'])
        pop.evolve(
            initOps = [InitSex()],
            matingScheme = RandomMating(subPopSize=demo),
            postOps = PyOperator(func=demoSize),
            gen = 100
        )


    def testNumOffspring(self):
        'Testing number of offspring'
        # single number
        self.assertEqual(
             self.getFamSize(numOffspring=5, N=50),
             [5]*10)

        # Python function
        def nos(gen):
            return gen%2+1
        self.assertEqual(
            self.getFamSize(numOffspring=nos, gen=2),
            [2]*500)
        self.assertEqual(
            self.getFamSize(numOffspring=nos, gen=3),
            [1]*1000)
        # a generator
        def nos_gen():
            while True:
                yield random.randint(5, 10)
        cnt = self.getFamSize(numOffspring=nos_gen, N=1000)
        self.assertEqual(sum(cnt), 1000)
        num = [ cnt.count(i) for i in range(1, 4) ]
        #
        mean = sum(num)/3.
        # for i in range(3):
        #     self.assertTrue(num[i] < mean + 50, 
        #     "Expression num[i] (test value %f) be less than mean + 50. This test may occasionally fail due to the randomness of outcome." % (num[i]))
        #     self.assertTrue(num[i] > mean - 50, 
        #     "Expression num[i] (test value %f) be greater than to mean - 50. This test may occasionally fail due to the randomness of outcome." % (num[i]))
        # randomnumber
        def nos():
            return random.randint(1, 3)
        cnt = self.getFamSize(numOffspring=nos, N=1000)
        self.assertEqual(sum(cnt), 1000)
        num = [ cnt.count(i) for i in range(1, 4) ]
        #
        mean = sum(num)/3.
        # for i in range(3):
        #     self.assertTrue(num[i] < mean + 50, 
        #     "Expression num[i] (test value %f) be less than mean + 50. This test may occasionally fail due to the randomness of outcome." % (num[i]))
        #     self.assertTrue(num[i] > mean - 50, 
        #     "Expression num[i] (test value %f) be greater than to mean - 50. This test may occasionally fail due to the randomness of outcome." % (num[i]))
        # GEOMETRIC_DISTRIBUTION
        p = 0.33
        cnt = self.getFamSize( numOffspring=(GEOMETRIC_DISTRIBUTION, p), N=10000)
        # mean should be 1/p, variance (1-p)/(p*p)
        mean = sum(cnt)*1.0/len(cnt)
        var = sum([x*x for x in cnt])*1.0/len(cnt) - mean*mean
        self.assertEqual(abs(mean - 1./p) < 0.1, True)
        self.assertEqual(abs(var - (1-p)/(p*p)) < 1, True)
        # POISSON_DISTRIBUTION
        p = 3
        cnt = self.getFamSize( numOffspring=(POISSON_DISTRIBUTION, p), N=100000)
        mean = sum(cnt)*1.0/len(cnt)
        self.assertEqual(abs(mean - p/(1-math.exp(-p))) < 0.1, True)
        # BINOMIAL_DISTRIBUTION
        p = 0.3
        n = 10
        cnt = self.getFamSize( numOffspring=(BINOMIAL_DISTRIBUTION, p, n), N=10000)
        mean = sum(cnt)*1.0/len(cnt)
        self.assertEqual(abs(mean - ((n)*p/(1-(1-p)**n))) < 0.1, True)
        # UNIFORM_DISTRIBUTION
        a = 3
        b = 6
        cnt = self.getFamSize( numOffspring=(UNIFORM_DISTRIBUTION, a, b), N=10000)
        mean = sum(cnt)*1.0/len(cnt)
        var = sum([x*x for x in cnt])*1.0/len(cnt) - mean*mean
        self.assertEqual(abs(mean - (a + b)/2.) < 0.1, True)


    def checkSexMode(self, ms):
        pop = Population(size=[40])
        pop.evolve(initOps = InitSex(), matingScheme=ms, gen=1)
        # return individual sex as a string
        def sexChar(sex):
            if sex == MALE:
                return 'M'
            else:
                return 'F'
        return ''.join([sexChar(ind.sex()) for ind in pop.individuals()])

    def testSexMode(self):
        'Testing parameter sexMode of mating schemes'
        # noSex
        self.assertEqual(
            self.checkSexMode(RandomMating(sexMode=NO_SEX)),
            'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM')
        # NUM_OF_MALES
        self.assertEqual(
            self.checkSexMode(RandomMating(numOffspring=3,
            sexMode=(NUM_OF_MALES, 1))),
            'MFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFM')
        # NUM_OF_FEMALES
        self.assertEqual(
            self.checkSexMode(RandomMating(numOffspring=4,
            sexMode=(NUM_OF_FEMALES, 2))),
            'FFMMFFMMFFMMFFMMFFMMFFMMFFMMFFMMFFMMFFMM')
        # SEQUENCE_OF_SEX
        self.assertEqual(
            self.checkSexMode(RandomMating(numOffspring=4,
            sexMode=(SEQUENCE_OF_SEX, MALE, FEMALE, FEMALE))),
            'MFFMMFFMMFFMMFFMMFFMMFFMMFFMMFFMMFFMMFFM')
        # SEQUENCE_OF_SEX
        self.assertEqual(
            self.checkSexMode(RandomMating(numOffspring=4,
            sexMode=(GLOBAL_SEQUENCE_OF_SEX, MALE, FEMALE, FEMALE))),
            'MFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFM')
        # PROB_OF_MALES
        pop = Population(10000)
        simu = Simulator(pop)
        simu.evolve(
            initOps = [InitSex(), InitGenotype(freq=[0.5, 0.5])],
            matingScheme = RandomMating(sexMode=(PROB_OF_MALES, 0.3)),
            postOps = [
                Stat(numOfMales=True),
                # number of male should be variable, but not too much
                TerminateIf('numOfMales < 2500 or numOfMales > 3500'),
            ],
            gen = 10
        )
        self.assertEqual(simu.dvars(0).gen, 10)
        # Using a function.
        def sexFunc():
            return random.randint(1, 2)
        self.assertNotEqual(
            self.checkSexMode(RandomMating(numOffspring=(UNIFORM_DISTRIBUTION, 2, 6),
                sexMode=sexFunc)),
            'FMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFM')
        # Using a generator function
        def sexFunc():
            i = 0
            while True:
                random.random()
                i += 1
                if i % 2 == 0:
                    yield MALE
                else:
                    yield FEMALE
        self.assertEqual(
            self.checkSexMode(RandomMating(numOffspring=(UNIFORM_DISTRIBUTION, 2, 6),
                sexMode=sexFunc)),
            'FMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFM')
        

    def testMonoMating(self):
        'Testing monogemous mating scheme'
        pop = Population(size=[2000], loci=[3,5], infoFields=['father_idx', 'mother_idx'])
        initGenotype(pop, freq=[0.2, 0.3, 0.5])
        simu = Simulator(pop)
        simu.evolve(
            initOps = InitSex(sex=(MALE, FEMALE)), 
            matingScheme =  MonogamousMating(numOffspring=2, sexMode=(NUM_OF_MALES, 1),
                ops=[MendelianGenoTransmitter(), ParentsTagger()]),
            gen = 5)
        self.assertEqual(len(set(simu.population(0).indInfo('father_idx'))), 1000)
        self.assertEqual(len(set(simu.population(0).indInfo('mother_idx'))), 1000)
        pop = simu.extract(0)
        self.assertEqual([ind.sex() for ind in pop.individuals()], [1,2]*1000)
               
    def testHeteroMating(self):
        'Testing heterogeneous mating schemes'
        pop = Population(size=[10000, 10000], loci=[2], infoFields=['father_idx', 'mother_idx'])
        pop.evolve(
            initOps = InitSex(),
            matingScheme=HeteroMating(
                [RandomMating(numOffspring=2, subPops=0, ops=[MendelianGenoTransmitter(), ParentsTagger()]),
                RandomMating(numOffspring=4, subPops=1, ops=[MendelianGenoTransmitter(), ParentsTagger()])]),
            gen=10)      
        parents = [(x, y) for x, y in zip(pop.indInfo('mother_idx'),
            pop.indInfo('father_idx'))]
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
        pop = Population(size =[20000, 20000], loci=[2], infoFields=['father_idx', 'mother_idx'])
        pop.setVirtualSplitter(ProportionSplitter([0.2, 0.8]))
        pop.evolve(
            initOps = InitSex(),
            matingScheme = HeteroMating([
                RandomMating(numOffspring=1, subPops=[(0,0)], ops=[MendelianGenoTransmitter(), ParentsTagger()]),
                RandomMating(numOffspring=2, subPops=[(1,1)], ops=[MendelianGenoTransmitter(), ParentsTagger()]),
                ]),
            gen =10
        )
        parents = [(x, y) for x, y in zip(pop.indInfo('mother_idx'),
            pop.indInfo('father_idx'))]
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
         
    def testWeightingScheme(self):
        'Testing weighting schemes of heterogeneous mating schemes'
        pop = Population(size=[1000], loci=2, infoFields='mark')
        pop.setVirtualSplitter(RangeSplitter([[0, 500], [200, 1000]]))
        # weighting scheme of -0.5, 2, 3
        pop.evolve(
            initOps=InitSex(),
            matingScheme=HeteroMating([
                RandomMating(subPops=0, weight=-0.5,
                    ops=[InfoExec('mark=0'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 0)], weight=2,
                    ops=[InfoExec('mark=1'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 1)], weight=3,
                    ops=[InfoExec('mark=2'), MendelianGenoTransmitter()])
            ]),
            gen = 10
        )
        marks = list(pop.indInfo('mark'))
        self.assertEqual(marks.count(0.), 500)
        self.assertEqual(marks.count(1.), 200)
        self.assertEqual(marks.count(2.), 300)
        # weighting scheme of 1, 2, 7
        pop.evolve(
            initOps=InitSex(),
            matingScheme=HeteroMating([
                RandomMating(subPops=0, weight=1,
                    ops=[InfoExec('mark=0'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 0)], weight=2,
                    ops=[InfoExec('mark=1'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 1)], weight=7,
                    ops=[InfoExec('mark=2'), MendelianGenoTransmitter()])
            ]),
            gen = 10
        )
        marks = list(pop.indInfo('mark'))
        self.assertEqual(marks.count(0.), 100)
        self.assertEqual(marks.count(1.), 200)
        self.assertEqual(marks.count(2.), 700)
        # weighting scheme of 0, 2, 3
        pop.evolve(
            initOps=InitSex(),
            matingScheme=HeteroMating([
                RandomMating(subPops=0, weight=0,
                    ops=[InfoExec('mark=0'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 0)], weight=2,
                    ops=[InfoExec('mark=1'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 1)], weight=3,
                    ops=[InfoExec('mark=2'), MendelianGenoTransmitter()])
            ]),
            gen = 10
        )
        marks = list(pop.indInfo('mark'))
        self.assertEqual(marks.count(0.), 0)
        self.assertEqual(marks.count(1.), 400)
        self.assertEqual(marks.count(2.), 600)
        # weighting scheme of -1, 0
        pop.evolve(
            initOps=InitSex(),
            matingScheme=HeteroMating([
                RandomMating(subPops=[(0, 0)], weight=-1,
                    ops=[InfoExec('mark=0'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 1)],
                    ops=[InfoExec('mark=1'), MendelianGenoTransmitter()])
            ]),
            gen = 10
        )
        marks = list(pop.indInfo('mark'))
        self.assertEqual(marks.count(0.), 500)
        self.assertEqual(marks.count(1.), 500)

    def testWeightingSchemeBySex(self):
        'Testing weighting schemes of heterogeneous mating schemes'
        pop = Population(size=[1000], loci=2, infoFields='mark')
        pop.setVirtualSplitter(RangeSplitter([[0, 500], [200, 1000]]))
        # weighting scheme of -0.5, 2, 3
        pop.evolve(
            initOps=InitSex(sex=[MALE, FEMALE]),
            matingScheme=HeteroMating([
                RandomMating(subPops=0, weight=-0.5,
                    ops=[InfoExec('mark=0'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 0)], weight=1,
                    ops=[InfoExec('mark=1'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 1)], weight=2,
                    ops=[InfoExec('mark=2'), MendelianGenoTransmitter()])
            ], weightBy=MALE_ONLY),
            gen = 1
        )
        marks = list(pop.indInfo('mark'))
        self.assertEqual(marks.count(0.), 250)
        self.assertEqual(marks.count(1.), 250)
        self.assertEqual(marks.count(2.), 500)
        #
        pop = Population(size=[1000], loci=2, infoFields='mark')
        pop.setVirtualSplitter(SexSplitter())
        # weighting scheme of -0.5, 2, 3
        #turnOnDebug('DBG_DEVEL,DBG_WARNING')
        pop.evolve(
            initOps=InitSex(sex=[MALE, FEMALE]),
            matingScheme=HeteroMating([
                RandomMating(subPops=0, weight=-0.5,
                    ops=[InfoExec('mark=0'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 0)], weight=1,
                    ops=[InfoExec('mark=1'), MendelianGenoTransmitter()]),
                RandomMating(subPops=0, weight=2,
                    ops=[InfoExec('mark=2'), MendelianGenoTransmitter()])
            ], weightBy=PAIR_ONLY),
            gen = 1
        )
        marks = list(pop.indInfo('mark'))
        self.assertEqual(marks.count(0.), 250)
        self.assertEqual(marks.count(1.), 0)
        self.assertEqual(marks.count(2.), 750)
        #
        pop.evolve(
            initOps=InitSex(sex=[MALE, FEMALE]),
            matingScheme=HeteroMating([
                RandomMating(subPops=0, weight=2,
                    ops=[InfoExec('mark=0'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 0)], weight=1,
                    ops=[InfoExec('mark=1'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(0, 1)], weight=2,
                    ops=[InfoExec('mark=2'), MendelianGenoTransmitter()])
            ], weightBy=PAIR_ONLY),
            gen = 1
        )
        marks = list(pop.indInfo('mark'))
        self.assertEqual(marks.count(0.), 1000)
        self.assertEqual(marks.count(1.), 0)
        self.assertEqual(marks.count(2.), 0)

    def testWeightByWithEmptySubpop(self):
        '''Test weightBy option when there are empty subpopulations'''
        pop = Population(size=[1000, 0], loci=2, infoFields='mark')
        pop.setVirtualSplitter(RangeSplitter([[0, 500], [200, 1000]]))
        # weighting scheme of -0.5, 2, 3
        pop.evolve(
            initOps=InitSex(sex=[MALE, FEMALE]),
            matingScheme=HeteroMating([
                RandomMating(weight=-0.5,
                    ops=[InfoExec('mark=0'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(ALL_AVAIL, 0)], weight=1,
                    ops=[InfoExec('mark=1'), MendelianGenoTransmitter()]),
                RandomMating(subPops=[(ALL_AVAIL, 1)], weight=2,
                    ops=[InfoExec('mark=2'), MendelianGenoTransmitter()])
            ], weightBy=MALE_ONLY),
            gen = 1
        )

    def testPolygamousMating(self):
        'Testing polygamous mating scheme'
        pop = Population(size=[200], loci=[3,5], infoFields=['father_idx', 'mother_idx'])
        initGenotype(pop, freq=[0.2, 0.3, 0.5])
        # exactly 100 males and 100 females
        for i in range(100):
            pop.individual(i).setSex(MALE)
            pop.individual(100+i).setSex(FEMALE)
        simu = Simulator(pop)
        simu.evolve(
            initOps = [],
            matingScheme = PolygamousMating(polySex=MALE, polyNum=3, numOffspring=2,ops=[MendelianGenoTransmitter(), ParentsTagger()]),
            gen = 1)
        # there is only one MALE...
        fi = simu.population(0).indInfo('father_idx')
        self.assertEqual(fi[0], fi[1])
        self.assertEqual(fi[0], fi[5])
        self.assertNotEqual(fi[0], fi[6])
        mi = simu.population(0).indInfo('mother_idx')
        self.assertEqual(mi[0], mi[1])
        self.assertNotEqual(mi[0], mi[2])              
       
    def testPedigreeMating(self):
        'Testing pedigree mating using a population object'
        pop = Population(size=[100, 100], loci=[2, 5], ancGen=-1,
            infoFields=['father_idx', 'mother_idx'])
        pop.evolve(
            initOps = InitSex(),
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(), ParentsTagger()]),
            gen = 20
        )

    def testSequentialParentsChooser(self):
        'Testing sequential parent chooser'
        pop = Population(size=[100, 200], infoFields=['parent_idx'])
        initSex(pop)
        initGenotype(pop, freq=[.3, .7])
        pop.evolve(
            matingScheme = HomoMating(
                SequentialParentsChooser(),
                OffspringGenerator(ops=[
                    SelfingGenoTransmitter(),
                    ParentsTagger(infoFields='parent_idx'), 
                    ])),
            gen=1)

    def testRandomParentChooser(self):
        'Testing sequential parent chooser'
        def traj(gen):
            return [0.5 + gen*0.01]
        pop = Population(size=[1000, 2000], infoFields=['parent_idx'])
        initGenotype(pop, freq=[.2, .8])
        pop.evolve(
            matingScheme= HomoMating(
                RandomParentChooser(),
                OffspringGenerator([SelfingGenoTransmitter(),
                    ParentsTagger(infoFields='parent_idx'), 
                ])),
            gen=1)

    def testPyParentsChooserRetValue(self):
        'Testing the return value of Python parents chooser'
        import random
        def retIndex(pop, subPop):
            while True:
                yield random.randint(0, pop.subPopSize(subPop) - 1)
        def retIndexes(pop, subPop):
            while True:
                yield random.randint(0, pop.subPopSize(subPop) - 1), random.randint(0, pop.subPopSize(subPop) -1)
        def retInd(pop, subPop):
            while True:
                yield pop.individual(random.randint(0, pop.subPopSize(subPop) - 1))
        def retInds(pop, subPop):
            while True:
                yield pop.individual(random.randint(0, pop.subPopSize(subPop) - 1)), \
                     pop.individual(random.randint(0, pop.subPopSize(subPop) - 1))
        def retVarInds(pop, subPop):
            while True:
                male = pop.individual(random.randint(0, pop.subPopSize(subPop) - 1))
                female = pop.individual(random.randint(0, pop.subPopSize(subPop) - 1))
                yield male, female
        def retPop(pop, subPop):
            while True:
                yield pop
        def retWrongIndex(pop, subPop):
            while True:
                yield pop.subPopSize(subPop)
        def retWrongIndexes(pop, subPop):
            while True:
                yield 0, pop.subPopSize(subPop)
        def testPyRetValue(func):
            pop = Population([200]*5)
            pop.evolve(
                matingScheme= HomoMating(
                    PyParentsChooser(func),
                    OffspringGenerator(CloneGenoTransmitter())
                ),
                gen = 5
            )
        testPyRetValue(retIndex)
        testPyRetValue(retIndexes)
        testPyRetValue(retInd)
        testPyRetValue(retInds)
        testPyRetValue(retVarInds)
        self.assertRaises(ValueError, testPyRetValue, retPop)
        self.assertRaises(ValueError, testPyRetValue, retWrongIndex)
        self.assertRaises(ValueError, testPyRetValue, retWrongIndexes)

  
    def testHaploidRandomMating(self):
        'Testing random mating in haploid populations'
        pop = Population(size=[50, 100], loci=[5]*5, ploidy=1,
            chromTypes=[CUSTOMIZED]*5)
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            initOps = [InitSex(),
                # female has [1]
                InitGenotype(genotype=[1]*25, subPops=[(0, 1), (1, 1)]),
                ],
            matingScheme=RandomMating(ops=[MitochondrialGenoTransmitter()]),
            gen = 1
        )
        self.assertEqual(pop.genotype(), [1]*(150*25))

    def testFixedCondMating(self):
        'Testing ConditionalMating with constant condition'
        pop = Population(size=100)
        pop.evolve(
            matingScheme=ConditionalMating(False, RandomMating(),
                RandomSelection()),
            gen=5)
        # forgot init sex, RandomMating will fail
        self.assertRaises(RuntimeError, pop.evolve,
            matingScheme=ConditionalMating(True, RandomMating(),
                RandomSelection()),
            gen=5)

    def testExprCondMating(self):
        'Testing ConditionalMating with constant condition'
        simu = Simulator(Population(size=1000), rep=2)
        simu.evolve(
            preOps=InitSex(),
            matingScheme=ConditionalMating('rep==0',
                RandomMating(),
                RandomMating(sexMode=(PROB_OF_MALES, 0.7))),
            finalOps=Stat(numOfMales=True),
            gen=5)
        self.assertLess(simu.dvars(0).numOfMales, 550)
        self.assertGreater(simu.dvars(1).numOfMales, 650)

    def cmFunc(self, pop):
        return pop.dvars().rep == 0

    def testFuncCondMating(self):
        'Testing ConditionalMating with constant condition'
        
        simu = Simulator(Population(size=1000), rep=2)
        simu.evolve(
            preOps=InitSex(),
            matingScheme=ConditionalMating(self.cmFunc,
                RandomMating(),
                RandomMating(sexMode=(PROB_OF_MALES, 0.7))),
            finalOps=Stat(numOfMales=True),
            gen=5)
        self.assertLess(simu.dvars(0).numOfMales, 550)
        self.assertGreater(simu.dvars(1).numOfMales, 650)

    def testSequentialParentChooser(self):
        'Testing SequentialParentChooser'
        pop = Population([10, 10], loci=[1], infoFields='a')
        initInfo(pop, range(10, 20), subPops=1, infoFields='a')
        c = SequentialParentChooser()
        c.initialize(pop=pop, subPop=1)
        idx = 0
        while idx < 10:
            ind, tmp = c.chooseParents()
            self.assertEqual(ind.a, 10 + idx)
            idx += 1

    def testRandomParentChooser(self):
        'Test random parent chooser'
        pop = Population([10, 10], loci=[1], infoFields='a')
        initInfo(pop, range(10), subPops=0, infoFields='a')
        initInfo(pop, range(10, 20), subPops=1, infoFields='a')
        c = RandomParentChooser()
        c.initialize(pop=pop, subPop=1)
        for idx in range(100):
            ind, tmp = c.chooseParents()
            self.assertEqual(tmp, None)
            self.assertLess(ind.a, 20)
            self.assertGreaterEqual(ind.a, 10)
        #
        c = RandomParentChooser()
        c.initialize(pop=pop, subPop=0)
        for idx in range(100):
            ind, tmp = c.chooseParents()
            self.assertEqual(tmp, None)
            self.assertLess(ind.a, 10)
            self.assertGreaterEqual(ind.a, 0)
        # fixed sex?
        initSex(pop, sex=[MALE, FEMALE])
        c = RandomParentChooser(sexChoice=MALE_ONLY)
        c.initialize(pop=pop, subPop=0)
        for idx in range(100):
            f, tmp = c.chooseParents()
            self.assertEqual(f.sex(), MALE)
            # only males are chosen
            self.assertEqual(int(f.a) % 2, 0)
            self.assertLess(f.a, 10)
            self.assertGreaterEqual(f.a, 0)

    def testRandomParentChooserWithFitness(self):
        'Test random parent chooser'
        pop = Population([10, 10], loci=[1], infoFields=['a', 'fitness'])
        initInfo(pop, range(10), subPops=0, infoFields='a')
        initInfo(pop, [0.5]*5 + [0]*5, subPops=0, infoFields='fitness')
        initInfo(pop, range(10, 20), subPops=1, infoFields='a')
        initInfo(pop, [1]*5 + [0]*5, subPops=1, infoFields='fitness')
        c = RandomParentChooser()
        c.initialize(pop=pop, subPop=1)
        for idx in range(500):
            ind, tmp = c.chooseParents()
            self.assertEqual(tmp, None)
            self.assertLess(ind.a, 15)
            self.assertGreaterEqual(ind.a, 10)
        #
        c = RandomParentChooser()
        c.initialize(pop=pop, subPop=0)
        for idx in range(500):
            ind, tmp = c.chooseParents()
            self.assertEqual(tmp, None)
            self.assertLess(ind.a, 5)
            self.assertGreaterEqual(ind.a, 0)

    def testRandomParentChooserWithoutReplacement(self):
        'Test random parent chooser'
        pop = Population([10, 10], loci=[1], infoFields='a')
        initInfo(pop, range(10), subPops=0, infoFields='a')
        initInfo(pop, range(10, 20), subPops=1, infoFields='a')
        c = RandomParentChooser(replacement=False)
        c.initialize(pop=pop, subPop=1)
        for idx in range(10):
            ind, tmp = c.chooseParents()
            self.assertEqual(tmp, None)
            self.assertLess(ind.a, 20)
            self.assertGreaterEqual(ind.a, 10)
        self.assertRaises(RuntimeError, c.chooseParents)
        #
        c = RandomParentChooser(replacement=False)
        c.initialize(pop=pop, subPop=0)
        for idx in range(10):
            ind, tmp = c.chooseParents()
            self.assertEqual(tmp, None)
            self.assertLess(ind.a, 10)
            self.assertGreaterEqual(ind.a, 0)
        self.assertRaises(RuntimeError, c.chooseParents)

    def testRandomParentsChooser(self):
        'Test random parent chooser'
        pop = Population([10, 10], loci=[1], infoFields='a')
        initSex(pop)
        initInfo(pop, range(10), subPops=0, infoFields='a')
        initInfo(pop, range(10, 20), subPops=1, infoFields='a')
        c = RandomParentsChooser()
        c.initialize(pop=pop, subPop=1)
        for idx in range(100):
            f, m = c.chooseParents()
            self.assertEqual(f.sex(), MALE)
            self.assertEqual(m.sex(), FEMALE)
            self.assertLess(f.a, 20)
            self.assertGreaterEqual(f.a, 10)
            self.assertLess(m.a, 20)
            self.assertGreaterEqual(m.a, 10)
        #
        c = RandomParentsChooser()
        c.initialize(pop=pop, subPop=0)
        for idx in range(100):
            f, m = c.chooseParents()
            self.assertEqual(f.sex(), MALE)
            self.assertEqual(m.sex(), FEMALE)
            self.assertLess(f.a, 10)
            self.assertGreaterEqual(f.a, 0)
            self.assertLess(m.a, 10)
            self.assertGreaterEqual(m.a, 0)


    def testRandomParentsChooserWithFitness(self):
        'Test random parent chooser with Fitness'
        pop = Population([10, 10], loci=[1], infoFields=['a', 'fitness'])
        #
        # Note: if we ignore sex=[MALE, FEMALE]
        # there is a case when there is no female in the first 5 individual, therefore
        # all females have fitness 0, and a female out of boundary could be selected
        initSex(pop, sex=[MALE, FEMALE])
        initInfo(pop, range(10), subPops=0, infoFields='a')
        initInfo(pop, [0.5]*5 + [0]*5, subPops=0, infoFields='fitness')
        initInfo(pop, range(10, 20), subPops=1, infoFields='a')
        initInfo(pop, [1]*5 + [0]*5, subPops=1, infoFields='fitness')
        #print([ind.sex() for ind in pop.individuals()])
        c = RandomParentsChooser()
        c.initialize(pop=pop, subPop=1)
        for idx in range(5000):
            f, m = c.chooseParents()
            self.assertEqual(f.sex(), MALE)
            self.assertEqual(m.sex(), FEMALE)
            self.assertLess(f.a, 15)
            self.assertGreaterEqual(f.a, 10)
            self.assertLess(m.a, 15)
            self.assertGreaterEqual(m.a, 10)
        #
        c = RandomParentsChooser()
        c.initialize(pop=pop, subPop=0)
        for idx in range(5000):
            f, m = c.chooseParents()
            self.assertEqual(f.sex(), MALE)
            self.assertEqual(m.sex(), FEMALE)
            self.assertLess(f.a, 5)
            self.assertGreaterEqual(f.a, 0)
            self.assertLess(m.a, 5)
            self.assertGreaterEqual(m.a, 0)

    def testRandomParentsChooserWithoutReplacement(self):
        'Test random parents chooser without replacement'
        pop = Population([10, 10], loci=[1], infoFields='a')
        initSex(pop, sex=[MALE, FEMALE])
        initInfo(pop, range(10), subPops=0, infoFields='a')
        initInfo(pop, range(10, 20), subPops=1, infoFields='a')
        c = RandomParentsChooser(replacement=False)
        c.initialize(pop=pop, subPop=1)
        for idx in range(5):
            f, m = c.chooseParents()
            self.assertEqual(f.sex(), MALE)
            self.assertEqual(m.sex(), FEMALE)
            self.assertLess(f.a, 20)
            self.assertGreaterEqual(f.a, 10)
            self.assertLess(m.a, 20)
            self.assertGreaterEqual(m.a, 10)
        self.assertRaises(RuntimeError, c.chooseParents)
        #
        c = RandomParentsChooser(replacement=False)
        c.initialize(pop=pop, subPop=0)
        for idx in range(5):
            f, m = c.chooseParents()
            self.assertEqual(f.sex(), MALE)
            self.assertEqual(m.sex(), FEMALE)
            self.assertLess(f.a, 10)
            self.assertGreaterEqual(f.a, 0)
            self.assertLess(m.a, 10)
            self.assertGreaterEqual(m.a, 0)
        self.assertRaises(RuntimeError, c.chooseParents)


    def testPolyParentsChooser(self):
        'Test PolyParentsChooser'
        pop = Population([10, 10], loci=[1], infoFields='a')
        initSex(pop)
        initInfo(pop, range(10), subPops=0, infoFields='a')
        initInfo(pop, range(10, 20), subPops=1, infoFields='a')
        c = PolyParentsChooser(polySex=FEMALE, polyNum=2)
        c.initialize(pop=pop, subPop=1)
        for idx in range(100):
            f, m = c.chooseParents()
            if idx % 2 == 0:
                lastInfo = m.a
            else:
                self.assertEqual(lastInfo, m.a)
            self.assertEqual(f.sex(), MALE)
            self.assertEqual(m.sex(), FEMALE)
            self.assertLess(f.a, 20)
            self.assertGreaterEqual(f.a, 10)
            self.assertLess(m.a, 20)
            self.assertGreaterEqual(m.a, 10)
        #
        c = PolyParentsChooser(polySex=MALE, polyNum=2)
        c.initialize(pop=pop, subPop=0)
        for idx in range(100):
            f, m = c.chooseParents()
            if idx % 2 == 0:
                lastInfo = f.a
            else:
                self.assertEqual(lastInfo, f.a)
            self.assertEqual(f.sex(), MALE)
            self.assertEqual(m.sex(), FEMALE)
            self.assertLess(f.a, 10)
            self.assertGreaterEqual(f.a, 0)
            self.assertLess(m.a, 10)
            self.assertGreaterEqual(m.a, 0)

    def testHermaphroditicMating(self):
        'Test HermaphroditicMating'
        pop = Population(100, loci=10)
        pop.evolve(
            matingScheme=HermaphroditicMating(),
            gen=10)
        #
        def theSame(dad, mom):
            if dad.ind_id == mom.ind_id:
                raise RuntimeError()
            return True
        #
        pop = Population(100, loci=10, infoFields='ind_id')
        pop.evolve(
            preOps=IdTagger(),
            matingScheme=HermaphroditicMating(allowSelfing=False,
                ops=[MendelianGenoTransmitter(),
                    IdTagger(),
                    PyOperator(func=theSame),
                    ]),
            gen=100)

if __name__ == '__main__':
    unittest.main()
    sys.exit(0)


