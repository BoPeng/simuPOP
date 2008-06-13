#!/usr/bin/env python
#
# This is a unittest file for migration operators
#
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision$
# $LastChangedDate$
# 

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, time, exceptions

class TestMigrator(unittest.TestCase):
    
    def testMigrateByCounts(self):
        'Testing migrate by counts'
        pop = population(size=[2000,4000,4000], loci=[2])
        Migrate(pop, mode=MigrByCounts, 
            rate = [ [0, 50, 50],
                             [0, 0, 0],
                             [0, 0, 0] ])
        assert pop.subPopSizes() == (1900, 4050, 4050)
        Migrate(pop, mode=MigrByCounts, 
            rate = [ [0, 50, 50],
                             [50, 0, 0],
                             [50, 0, 0] ])
        assert pop.subPopSizes() == (1900, 4050, 4050)
        # rate should be a matrix
        self.assertRaises(exceptions.ValueError, 
            Migrate, pop, [ [0, 50],
                             [50, 0, 0],
                             [50, 0, 0] ], MigrByCounts)
        
    def testMigrateByProportion(self):
        'Testing migrate by proportion'
        pop = population(size=[2000,4000,4000], loci=[2])
        # now if we want to inject a mutation whenever fixation happens
        Migrate(pop, mode=MigrByProportion, 
            rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ])
        assert pop.subPopSizes() == (2000, 4000, 4000)
        Migrate(pop, mode=MigrByProportion, 
            rate = [ [0, .25, .25],
                             [0.25, 0, 0],
                             [0, 0.25, 0] ])
        assert pop.subPopSizes() == (2000, 4500, 3500)
        
    def testMigrateByProbability(self):
        'Testing migrate by probability'
        pop = population(size=[2000,4000,4000], loci=[2])
        # now if we want to inject a mutation whenever fixation happens
        Migrate(pop, mode=MigrByProbability, 
            rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ])
        # print pop.subPopSizes()
        assert abs(pop.subPopSize(0) - 2000) < 100
        assert abs(pop.subPopSize(1) - 4000) < 100
        assert abs(pop.subPopSize(2) - 4000) < 100
        Migrate(pop, mode=MigrByProbability, 
            rate = [ [0, .25, .25],
                             [0.25, 0, 0],
                             [0, 0.25, 0] ])
        # print pop.subPopSizes()
        assert abs(pop.subPopSize(0) - 2000) < 100
        assert abs(pop.subPopSize(1) - 4500) < 100
        assert abs(pop.subPopSize(2) - 3500) < 100

    def testMigrateFromTo(self):
        'Testing parameter from and to of migrators'
        pop = population(size=[2000,4000,4000], loci=[2])
        # now if we want to inject a mutation whenever fixation happens
        Migrate(pop, mode=MigrByProbability, 
            fromSubPop = [0], toSubPop = [1,2], 
            rate = [.05, .05] )
        # print pop.subPopSizes()
        assert abs(pop.subPopSize(0) - 1800) < 50
        assert abs(pop.subPopSize(1) - 4100) < 50
        assert abs(pop.subPopSize(2) - 4100) < 50
        # other parameter form can be used as well
        pop = population(size=[2000,4000,4000], loci=[2])
        Migrate(pop, mode=MigrByProbability, 
            fromSubPop = 0, toSubPop = [1,2], 
            rate = [[.05, .05]] )
        assert abs(pop.subPopSize(0) - 1800) < 50
        assert abs(pop.subPopSize(1) - 4100) < 50
        assert abs(pop.subPopSize(2) - 4100) < 50

    def testMigrateBySexAndCounts(self):
        'Testing migrate by sex and counts'
        # everyone is Male
        pop = population(size=[2000,4000,4000], loci=[2])
        Migrate(pop, mode=MigrByCounts, 
            rate = [ [0, 50, 50],
                             [0, 0, 0],
                             [0, 0, 0] ])
        assert pop.subPopSizes() == (1900, 4050, 4050)
        Migrate(pop, mode=MigrByCounts, 
            rate = [ [0, 50, 50],
                             [50, 0, 0],
                             [50, 0, 0] ])
        assert pop.subPopSizes() == (1900, 4050, 4050)
        # rate should be a matrix
        self.assertRaises(exceptions.ValueError, 
            Migrate, pop, [ [0, 50],
                             [50, 0, 0],
                             [50, 0, 0] ], MigrByCounts)
        
    def testMigrateBySexAndProportion(self):
        'Testing migrate by sex and proportion'
        pop = population(size=[2000,4000,4000], loci=[2])
        # now if we want to inject a mutation whenever fixation happens
        Migrate(pop, mode=MigrByProportion, 
            rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ])
        assert pop.subPopSizes() == (2000, 4000, 4000)
        Migrate(pop, mode=MigrByProportion, 
            rate = [ [0, .25, .25],
                             [0.25, 0, 0],
                             [0, 0.25, 0] ])
        assert pop.subPopSizes() == (2000, 4500, 3500)
        
    def testMigrateBySexAndProbability(self):
        'Testing migrate by sex and probability'
        pop = population(size=[2000,4000,4000], loci=[2])
        # now if we want to inject a mutation whenever fixation happens
        Migrate(pop, mode=MigrByProbability, 
            rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ])
        # print pop.subPopSizes()
        assert abs(pop.subPopSize(0) - 2000) < 100
        assert abs(pop.subPopSize(1) - 4000) < 100
        assert abs(pop.subPopSize(2) - 4000) < 100
        Migrate(pop, mode=MigrByProbability, 
            rate = [ [0, .25, .25],
                             [0.25, 0, 0],
                             [0, 0.25, 0] ])
        # print pop.subPopSizes()
        assert abs(pop.subPopSize(0) - 2000) < 100
        assert abs(pop.subPopSize(1) - 4500) < 100
        assert abs(pop.subPopSize(2) - 3500) < 100

    
    def testMigrConstAlleleFreq(self):
        'Testing that migration does not change allele frequency'
        pop = population(size=[2000,4000,4000], loci=[2])
        InitByFreq(pop, [.2, .8])
        Stat(pop, alleleFreq=[0])
        af = pop.dvars().alleleFreq[0][1]    # ~.2
        # migrate and check if allele frequency changes
        Migrate(pop, mode=MigrByProbability, 
            fromSubPop = [0], toSubPop = [1,2], 
            rate = [.05, .05] )
        Stat(pop, alleleFreq=[0])
        assert pop.dvars().alleleFreq[0][1] == af
        Migrate(pop, mode=MigrByProbability, 
            rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ])
        assert pop.dvars().alleleFreq[0][1] == af
        Migrate(pop, mode=MigrByProbability, 
            rate = [ [0, .25, .25],
                             [0.25, 0, 0],
                             [0, 0.25, 0] ])
        assert pop.dvars().alleleFreq[0][1] == af
             
    def testPyMigrator(self):
        'Testing operator pyMigrator'
        #pop = population(size=[2,4,4], loci=[2,6])
        #InitByFreq(pop, [.2,.4,.4])
        #ind1 = pop.individual(0).arrGenotype(True)
        #PyMigrate(pop, subPopID=[2,2,0,0,2,2,1,1,1,1])
        #assert ind1 == pop.individual(6).arrGenotype(True)
        #PyMigrate(pop, subPopID=[0,0,2,2,2,2,1,1,1,1])
        #assert ind1 == pop.individual(2).arrGenotype(True)
        #self.assertRaises(exceptions.ValueError, 
        #    PyMigrate, pop, [0,0,2,2,2,2,1,1])
        #TurnOnDebug(DBG_MIGRATOR)
        def rFunc(gen, size = []):
            if gen == 0:
                return [ [0.1] ]
            elif gen == 1:
                return [ [0.5] ]
            else:
                return [[0]]
        pop = population(size=[2000,4000], loci=[2])
        InitByFreq(pop, [.2, .8])
        simu = simulator(pop, randomMating())
        simu.evolve(
            preOps = [],
            ops = [
                pyMigrator(rateFunc = rFunc, fromSubPop = [0],
                           toSubPop = [1], mode = MigrByProportion),
                stat( popSize = True ),
                # pyEval('subPopSize'),
                terminateIf( "subPopSize != [1800, 4200]", at = [0]),
                terminateIf( "subPopSize != [900, 5100]", at = [1])
            ],
            end=5
        )
        self.assertEqual(simu.gen(), 6)
        
    def testSplitSubPop(self):
        'Testing population split'
        pop = population(size=10, loci=[2,6])
        InitByFreq(pop, [.2,.4,.4])
        genotype = list(pop.arrGenotype(True))
        SplitSubPop(pop, which=0, sizes=[2,8], randomize=False)
        # individual untouched
        self.assertEqual(pop.arrGenotype(True), genotype)
        # split, with randomization
        SplitSubPop(pop, which=1, sizes=[6,2], randomize=True)
        self.assertNotEqual(pop.arrGenotype(True), genotype)


    def testRearrange(self):
        'Testing if info and genotype are migrated with individuals'
        if AlleleType() == 'binary':
            return
        #TurnOnDebug(DBG_POPULATION)
        pop = population(size=[4, 6], loci=[1], infoFields=['a','b'])
        pop.arrGenotype(True)[:] = range(20)
##         pop.arrIndInfo()[:] = range(20)
##         self.assertEqual([x.allele(0,0) for x in pop.individuals(0)],
##             [int(x) for x in pop.indInfo('a', 0)])
##         self.assertEqual([x.allele(0,0) for x in pop.individuals(1)],
##             [int(x) for x in pop.indInfo('a', 1)])
        #PyMigrate(pop, subPopID=[1,1,1,1,0,0,0,0,0,0])
        # from ind, this should be fine.
##         self.assertEqual([x.allele(0,0) for x in pop.individuals(0)],
##             [x.info('a') for x in pop.individuals(0)])
##         self.assertEqual([x.allele(0,0) for x in pop.individuals(1)],
##             [x.info('a') for x in pop.individuals(1)])
        # Is info rearranged?
        #self.assertequal([x.allele(0,0) for x in pop.individuals(0)],
        #    [int(x) for x in pop.indInfo('a', 0)])
        # if we force reorder
##         self.assertEqual([x.allele(0,0) for x in pop.individuals(0)],
##             [int(x) for x in pop.indInfo('a', 0)])
##         return

    def testMergeSubPop(self):
        'Testing population merge'
        pop = population(size=[2,4,4], loci=[2,6])
        MergeSubPops(pop, subPops=[0,2])
        self.assertEqual(pop.subPopSizes(), (6,4,0))
        MergeSubPops(pop, subPops=[0,1], removeEmptySubPops=True)
        self.assertEqual(pop.subPopSizes(), (10,))
    
    def testResizeSubPops(self):
        'Testing population resize'
        pop = population(size=[2, 4, 4], loci=[2,6])
        InitByFreq(pop, [0.3, 0.7])
        ResizeSubPops(pop, newSizes=[6], subPops=[0])
        self.assertEqual(pop.subPopSizes(), (6, 4, 4))
        for ind in (2, 4):
            self.assertEqual(pop.individual(ind, 0), pop.individual(0, 0))
            self.assertNotEqual(pop.individual(ind, 0), pop.individual(1, 0))
        for ind in (3, 5):
            self.assertEqual(pop.individual(ind, 0), pop.individual(1, 0))
            self.assertNotEqual(pop.individual(ind, 0), pop.individual(0, 0))
        # no propagate
        ResizeSubPops(pop, newSizes=[8, 7], subPops=[1,2], propagate=False)
        self.assertEqual(pop.subPopSizes(), (6, 8, 7))
        for ind in (10, 11, 12, 13, 18, 19, 20):
            self.assertEqual(pop.individual(ind, 0).arrGenotype(), [0]*16)
        
        
        
if __name__ == '__main__':
    unittest.main()
     



