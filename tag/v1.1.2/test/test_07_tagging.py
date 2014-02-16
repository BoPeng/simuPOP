#!/usr/bin/env python
#
#    This is a unittest file for taggers
#
#    Bo Peng (bpeng@rice.edu)
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

class TestTagger(unittest.TestCase):

    def testParentsTagger(self):
        'Testing parents tagger.'
        simu = Simulator(
            Population(size=[50,150], ploidy=2, loci=[2,4],
                    infoFields=['father_idx', 'mother_idx']))
        simu.evolve(
            initOps = [InitSex()],
            matingScheme = RandomMating(numOffspring=2, ops=[MendelianGenoTransmitter(), ParentsTagger()]),
            gen = 1
        )
        pop = simu.population(0)
        # check if all siblings have the same parents
        for sp in range(pop.numSubPop()):
            for i in range(int(pop.subPopSize(sp)/2)):
                self.assertEqual(pop.individual(i*2,sp).info(0),
                    pop.individual(i*2+1,sp).info(0) )
                self.assertEqual(pop.individual(i*2,sp).info(1),
                    pop.individual(i*2+1,sp).info(1) )
            # note that the last one may be left alone

    def testOffspringTagger(self):
        'Testing offspring tagger.'
        simu = Simulator(
            Population(size=[50,150], ploidy=2, loci=[2,4],
                    infoFields='offspring_idx'))
        simu.evolve(
            initOps = [InitSex()],
            matingScheme = RandomMating(numOffspring=3, ops=[MendelianGenoTransmitter(), OffspringTagger()]),
            gen = 1
        )
        pop = simu.population(0)
        self.assertEqual(pop.indInfo('offspring_idx', subPop=0), tuple([0., 1., 2.]*16 + [0., 1.]))
        # varying family size?
        pop = Population(size=[500,1000], ploidy=2, loci=[2,4],
                    infoFields=['offspring_idx', 'father_id', 'ind_id', 'mother_id'])
        pop.evolve(
            initOps = [InitSex(), IdTagger(), ],
            matingScheme = RandomMating(numOffspring=(UNIFORM_DISTRIBUTION, 2, 5),
                ops=[MendelianGenoTransmitter(), IdTagger(),
                PedigreeTagger(), OffspringTagger()]),
            gen = 1
        )
        lf = 0
        lm = 0
        li = 0
        for ind in pop.individuals():
            if ind.father_id == lf and ind.mother_id == lm:
                li += 1
            else:
                lf = ind.father_id
                lm = ind.mother_id
                li = 0
            self.assertEqual(ind.offspring_idx, li)


    def testInheritTagger(self):
        'Testing inherit tagger (pass info from parents to offspring'
        # this operator pass tag from one or both parents to offspring
        # the game is not:
        # who is the offspring of one parent?
        pop = Population(size=[50,150], ploidy=2, loci=[2,4],
                infoFields=['paternal_tag', 'maternal_tag'])
        pop.individual(0).setInfo(1, 'paternal_tag')
        pop.individual(50).setInfo(2, 'paternal_tag')
        simu = Simulator(pop)
        # other mode include mode=MATERNAL, TAG_Both
        simu.evolve(
            initOps = [InitSex()],
            matingScheme = RandomMating(ops=[MendelianGenoTransmitter(), InheritTagger(mode=PATERNAL)]),
            gen = 1)
        pop = simu.population(0)
        # we only know subpopulation 0 can not have tag 2
        # we only know subpopulation 1 can not have tag 1
        for i in range(pop.subPopSize(0)):
            self.assertNotEqual(pop.individual(i,0).info('paternal_tag'), 2)
        for i in range(pop.subPopSize(1)):
            self.assertNotEqual(pop.individual(i,1).info('paternal_tag'), 1)
        # from this test, we can see that genetic drift
        # can easily remove a signal (tag) from population.

    def testInheritTaggerToFile(self):
        'Testing inherit tagger that record indexes to a file'
        pop = Population(size=[50,150], ploidy=2, loci=[2,4],
                infoFields=['paternal_tag', 'maternal_tag'])
        for ind in pop.individuals(0):
            ind.setInfo(1, 'paternal_tag')
        for ind in pop.individuals(1):
            ind.setInfo(2, 'paternal_tag')
        simu = Simulator( pop)
        # other mode include mode=MATERNAL, TAG_Both
        simu.evolve(
            initOps = [InitSex()],
            matingScheme = RandomMating(ops=[MendelianGenoTransmitter(), InheritTagger(mode=PATERNAL)]),
            gen = 1
        )
        pop = simu.extract(0)
        # we only know subpopulation 0 can not have tag 2
        # we only know subpopulation 1 can not have tag 1
        for i in range(pop.subPopSize(0)):
            self.assertNotEqual( pop.individual(i,0).info('paternal_tag'), 2 )
        for i in range(pop.subPopSize(1)):
            self.assertNotEqual( pop.individual(i,1).info('paternal_tag'), 1 )

 
    def testPyTagger(self):
        'Testing python tagger (pass trait from parents to offspring)'
        pop = Population(size=[50,150], ploidy=2, loci=[2,4],
                infoFields=['trait1', 'trait2'])
        pop.setIndInfo([1], 'trait1')
        pop.setIndInfo([2], 'trait2')
        def myfunc(trait1, trait2):
            return [trait1[0]+trait1[1], trait2[0]*trait2[1]]
        #
        simu = Simulator(pop)
        simu.evolve(
            initOps = [InitSex()],
            matingScheme = RandomMating(ops=[MendelianGenoTransmitter(), 
                PyTagger(func=myfunc),
            ]),
            gen = 4)
        pop = simu.population(0)
        for ind in pop.individuals():
            # 1 + 1 = 2, 2 + 2 = 4, ...
            self.assertEqual(ind.info('trait1'), 16)
            # 2 * 2 = 4, 4 * 4 = 16, ...
            self.assertEqual(ind.info('trait2'), 65536)


    def TestPedigree(self):
        'Testing the handling of Pedigrees (FIXME)'
        pop = Population(size=[100, 100], loci=[2,5], infoFields=['x', 'y', 'z'])
        initGenotype(pop, freq=[0.2, 0.8])
        def addToZ(val):
            return [val[0]+1]
        simu = Simulator(pop)
        simu.evolve(
            matingScheme = RandomMating(ops=[MendelianGenoTransmitter(), 
                ParentsTagger(output='>>Pedigree.dat', infoFields=[]),
                PyTagger(output='>>z.dat', func=addToZ, infoFields=['z'])
                ]),
            end=10
        )
        return
        ped = Pedigree(pedfile='Pedigree.dat')
        ped.loadInfo('affection.dat', 'affection')
        ped.saveInfo('aff1.dat', 'affection')
        ped.loadInfo('info.dat', ['x', 'y'])
        ped.loadInfo('z.dat', 'z')
        ped.loadInfo('sex.dat', 'sex')
        ped.saveInfo('all.dat', ['affection', 'sex', 'x', 'y', 'z'])
        # reduce pedigree
        ped.markUnrelated()
        ped.removeUnrelated()
        ped.saveInfo('z1.dat', 'z')
        # add and set info?
        ped.addInfo('a')
        self.assertEqual(ped.info(0, 0, 'a'), 0)
        ped.addInfo('b', -1)
        self.assertEqual(ped.info(0, 0, 'b'), -1)
        #
        ped.addGen([20, 20])
        self.assertEqual(ped.subPopSizes(ped.gen()-1), (20, 20))
        ped.setFather(5, 11, 0)
        ped.setMother(5, 11, 0)
        self.assertRaises(IndexError, ped.setMother, 500, 11, 0)
        self.assertEqual(ped.father(11, 0), 5)
        self.assertEqual(ped.mother(11, 0), 5)


if __name__ == '__main__':
    unittest.main()
