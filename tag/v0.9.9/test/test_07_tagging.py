#!/usr/bin/env python
#
#    This is a unittest file for taggers
#
#    Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

class TestTagger(unittest.TestCase):

    def testParentsTagger(self):
        'Testing parents tagger.'
        simu = simulator(
            population(size=[50,150], ploidy=2, loci=[2,4],
                    infoFields=['father_idx', 'mother_idx']),
            randomMating(numOffspring=2))
        simu.evolve(
            initOps = [initSex()],
            duringOps =[parentsTagger()],
            gen = 1
        )
        pop = simu.population(0)
        # check if all siblings have the same parents
        for sp in range(pop.numSubPop()):
            for i in range(pop.subPopSize(sp)/2):
                self.assertEqual(pop.individual(i*2,sp).info(0),
                    pop.individual(i*2+1,sp).info(0) )
                self.assertEqual(pop.individual(i*2,sp).info(1),
                    pop.individual(i*2+1,sp).info(1) )
            # note that the last one may be left alone

    def testInheritTagger(self):
        'Testing inherit tagger (pass info from parents to offspring'
        # this operator pass tag from one or both parents to offspring
        # the game is not:
        # who is the offspring of one parent?
        pop = population(size=[50,150], ploidy=2, loci=[2,4],
                infoFields=['paternal_tag', 'maternal_tag'])
        pop.individual(0).setInfo(1, 'paternal_tag')
        pop.individual(50).setInfo(2, 'paternal_tag')
        simu = simulator(pop, randomMating())
        # other mode include mode=MATERNAL, TAG_Both
        simu.evolve(
            initOps = [initSex()],
            duringOps =[inheritTagger(mode=PATERNAL)],
            gen = 1)
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
        pop = population(size=[50,150], ploidy=2, loci=[2,4],
                infoFields=['paternal_tag', 'maternal_tag'])
        for ind in pop.individuals(0):
            ind.setInfo(1, 'paternal_tag')
        for ind in pop.individuals(1):
            ind.setInfo(2, 'paternal_tag')
        simu = simulator( pop, randomMating())
        # other mode include mode=MATERNAL, TAG_Both
        simu.evolve(
            initOps = [initSex()],
            duringOps =[inheritTagger(mode=PATERNAL)],
            gen = 1
        )
        # we only know subpopulation 0 can not have tag 2
        # we only know subpopulation 1 can not have tag 1
        for i in range(pop.subPopSize(0)):
            self.assertNotEqual( pop.individual(i,0).info('paternal_tag'), 2 )
        for i in range(pop.subPopSize(1)):
            self.assertNotEqual( pop.individual(i,1).info('paternal_tag'), 1 )

 
    def testPyTagger(self):
        'Testing python tagger (pass trait from parents to offspring)'
        pop = population(size=[50,150], ploidy=2, loci=[2,4],
                infoFields=['trait1', 'trait2'])
        pop.setIndInfo([1], 'trait1')
        pop.setIndInfo([2], 'trait2')
        def myfunc(values):
            'values are t1_pa, t2_pa, t1_mo, t2_mo'
            return [values[0]+values[2], values[1]*values[3]]
        #
        simu = simulator(pop, randomMating())
        simu.evolve(
            initOps = [initSex()],
            duringOps =[
                pyTagger(infoFields=['trait1', 'trait2'], func=myfunc),
            ],
            gen = 4)
        pop = simu.population(0)
        for ind in pop.individuals():
            # 1 + 1 = 2, 2 + 2 = 4, ...
            self.assertEqual(ind.info('trait1'), 16)
            # 2 * 2 = 4, 4 * 4 = 16, ...
            self.assertEqual(ind.info('trait2'), 65536)


    def TestPedigree(self):
        'Testing the handling of pedigrees (FIXME)'
        pop = population(size=[100, 100], loci=[2,5], infoFields=['x', 'y', 'z'])
        InitByFreq(pop, [0.2, 0.8])
        def addToZ(val):
            return [val[0]+1]
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps =[
                parentsTagger(output='>>pedigree.dat', infoFields=[]),
                pyTagger(output='>>z.dat', func=addToZ, infoFields=['z'])
                ],
            end=10
        )
        return
        ped = pedigree(pedfile='pedigree.dat')
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
        self.assertRaises(exceptions.IndexError, ped.setMother, 500, 11, 0)
        self.assertEqual(ped.father(11, 0), 5)
        self.assertEqual(ped.mother(11, 0), 5)


if __name__ == '__main__':
    unittest.main()