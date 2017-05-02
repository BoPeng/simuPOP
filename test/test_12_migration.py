#!/usr/bin/env python
#
# This is a unittest file for migration operators
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import unittest, os, sys, time
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

import numpy as np
from numpy.linalg import inv


def repeated(func, times=100, *args, **kwargs):
    accu = []
    for i in range(times):
        res = func(*args, **kwargs)
        if isinstance(res, (int, float)):
            accu.append(res)
        else:
            if not accu:
                accu = [[x] for x in res]
            else:
                [x.append(y) for x,y in zip(accu, res)]
    if isinstance(accu[0], (int, float)):
        return sum(accu)/len(accu)
    else:
        return [sum(x)/len(x) for x in accu]



def Sp_From_Forward_and_S(FM, S):
    # calculate predicted subpopulation size from forward migration matrix (FM)
    # and current population size (S)  S'=F^T S
    #
    # We need to normalize FM(i,i) so that row sum is 1
    return np.array(
        [[1-sum(x)+k if j==i else k for j,k in enumerate(x)]
        for i,x in enumerate(FM)]).transpose().dot(np.array(S))

# for example: the following equals 2000, 4000, 4000
# Sp_From_Forward_and_S([ [0, .05, .05],
#                         [0.025, 0, 0],
#                         [0.025, 0, 0] ],
#         [2000, 4000, 4000]) 
    
def BM_From_Forward_and_S(FM, S):
    # calculate backward migration matrix from forward migration matrix
    # and current population size
    #
    # B = diag(Sp)^-1 F^T diag(S)
    Sp = Sp_From_Forward_and_S(FM, S)
    Ft = np.array([[1-sum(x)+k if j==i else k for j,k in enumerate(x)]
        for i,x in enumerate(FM)]).transpose()
    return np.diag(1. / Sp).dot(Ft).dot(np.diag(S))
        
# the following get
# [[ 0.9  ,  0.05 ,  0.05 ],
#       [ 0.025,  0.975,  0.   ],
#       [ 0.025,  0.   ,  0.975]]
#
#BM_From_Forward_and_S([ [0, .05, .05],
#                         [0.025, 0, 0],
#                         [0.025, 0, 0] ],
#         [2000, 4000, 4000]) 

def Sp_From_Backward_and_S(BM, S):
    # Calculate predicted subpopulation size from backward migration matrix
    # and current population size: S'=B^T^-1 * S
    #
    # normalize
    B = np.array([[1-sum(x)+k if j==i else k for j,k in enumerate(x)]
        for i,x in enumerate(BM)])
    return inv(np.array(B.transpose())).dot(np.array(S))

# the following yield [2000, 4000, 4000]
#Sp_From_Backward_and_S([[ 0.9  ,  0.05 ,  0.05 ],
#       [ 0.025,  0.975,  0.   ],
#       [ 0.025,  0.   ,  0.975]], 
#    [2000, 4000, 4000])

def FM_From_Backward_and_S(BM, S):
    #
    # F = diag(S)^-1 B^T diag(S')
    #
    Sp = Sp_From_Backward_and_S(BM, S)
    Bt = np.array([[1-sum(x)+k if j==i else k for j,k in enumerate(x)]
        for i,x in enumerate(BM)]).transpose()
    return np.diag(1. / np.array(S)).dot(Bt).dot(np.diag(Sp))

# the following yields
# [[ 0.9  ,  0.05 ,  0.05 ],
#  [ 0.025,  0.975,  0.   ],
#  [ 0.025,  0.   ,  0.975]]

#FM_From_Backward_and_S([[ 0.9  ,  0.05 ,  0.05 ],
#       [ 0.025,  0.975,  0.   ],
#       [ 0.025,  0.   ,  0.975]], 
#    [2000, 4000, 4000])



class TestMigrator(unittest.TestCase):

    def testmigrateByCounts(self):
        'Testing migrate by counts'
        pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
        migrate(pop, mode=BY_COUNTS,
            rate = [ [0, 50, 50],
                             [0, 0, 0],
                             [0, 0, 0] ])
        self.assertEqual(pop.subPopSizes(), (1900, 4050, 4050))
        migrate(pop, mode=BY_COUNTS,
            rate = [ [0, 50, 50],
                             [50, 0, 0],
                             [50, 0, 0] ])
        self.assertEqual(pop.subPopSizes(), (1900, 4050, 4050))
        # rate should be a matrix
        self.assertRaises(ValueError,
            migrate, pop, [ [0, 50],
                             [50, 0, 0],
                             [50, 0, 0] ], BY_COUNTS)

    def testMigrateByProportion(self):
        'Testing migrate by proportion'
        pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
        # now if we want to inject a mutation whenever fixation happens
        self.assertRaises(ValueError, migrate, pop, rate=[1, 0.1])
        migrate(pop, mode=BY_PROPORTION,
            rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ])
        self.assertEqual(pop.subPopSizes(), (2000, 4000, 4000))
        migrate(pop, mode=BY_PROPORTION,
            rate = [ [0, .25, .25],
                             [0.25, 0, 0],
                             [0, 0.25, 0] ])
        self.assertEqual(pop.subPopSizes(), (2000, 4500, 3500))



    def testmigrateByIndInfo(self):
        'Testing migrate by indinfo'
        pop = Population(size=[2000, 4000, 4000], loci=[2], infoFields=['migrate_to'])
        for sp in range(3):
            pop.setIndInfo([sp], 'migrate_to', sp)
        migrate(pop, mode=BY_IND_INFO)
        self.assertEqual(pop.subPopSizes(), (2000, 4000, 4000))
        # migrate a few individuals?
        for i in range(10):
            pop.individual(i, 0).setInfo(1, 'migrate_to')
            pop.individual(i, 1).setInfo(2, 'migrate_to')
        migrate(pop, mode=BY_IND_INFO)
        self.assertEqual(pop.subPopSizes(), (1990, 4000, 4010))
        # virtual subpopulations?
        pop = Population(size=[2000, 4000, 4000], loci=[2], infoFields=['migrate_to'])
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        # only get male out of the second subpopulation
        migrate(pop, mode=BY_IND_INFO, subPops=[(1, 0)])
        self.assertEqual(pop.subPopSizes(), (4000, 2000, 4000))
        for ind in pop.individuals(1):
            self.assertEqual(ind.sex(), FEMALE)
        pop.setIndInfo([1], 'migrate_to', [0])
        # only 1000 females will go
        migrate(pop, mode=BY_IND_INFO, subPops=[(0, 1)])
        self.assertEqual(pop.subPopSizes(), (3000, 3000, 4000))
        for ind in pop.individuals(1):
            self.assertEqual(ind.sex(), FEMALE)
        for ind in pop.individuals(0):
            self.assertEqual(ind.sex(), MALE)



    def testmigrateByProbability(self):
        'Testing migrate by probability'
        def migrateSize():
            pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
            # now if we want to inject a mutation whenever fixation happens
            migrate(pop, mode=BY_PROBABILITY,
                rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ])
            return pop.subPopSizes()
        
        tested = repeated(migrateSize, times=100)
        # print pop.subPopSizes()
        self.assertTrue(abs(tested[0] - 2000) < 20)
        self.assertTrue(abs(tested[1] - 4000) < 20)
        self.assertTrue(abs(tested[2] - 4000) < 20)

        def migrateSize():
            pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
            migrate(pop, mode=BY_PROBABILITY,
            rate = [ [0, .25, .25],
                             [0.25, 0, 0],
                             [0, 0.25, 0] ])
            return pop.subPopSizes()
        # print pop.subPopSizes()
        tested = repeated(migrateSize, times=100)
        self.assertTrue(abs(tested[0] - 2000) < 20)
        self.assertTrue(abs(tested[1] - 4500) < 20)
        self.assertTrue(abs(tested[2] - 3500) < 20)

    def testmigrateFromTo(self):
        'Testing parameter from and to of Migrators'
        def migrateSize():
            pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
            # now if we want to inject a mutation whenever fixation happens
            migrate(pop, mode=BY_PROBABILITY,
                subPops = [0], toSubPops = [1,2],
                rate = [.05, .05] )
            return pop.subPopSizes()

        tested = repeated(migrateSize, times=100)
        self.assertTrue(abs(tested[0] - 1800) < 10)
        self.assertTrue(abs(tested[1] - 4100) < 10)
        self.assertTrue(abs(tested[2] - 4100) < 10)
        # other parameter form can be used as well

        def migrateSize():
            pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
            migrate(pop, mode=BY_PROBABILITY,
                subPops = 0, toSubPops = [1,2],
                rate = [[.05, .05]] )
            return pop.subPopSizes()
        tested = repeated(migrateSize, times=100)
        self.assertTrue(abs(tested[0] - 1800) < 10)
        self.assertTrue(abs(tested[1] - 4100) < 10)
        self.assertTrue(abs(tested[2] - 4100) < 10)


    def testmigrateBySexAndCounts(self):
        'Testing migrate by sex and counts'
        # everyone is MALE
        pop = Population(size=[2000, 4000,4000], loci=[2], infoFields=['migrate_to'])
        initSex(pop, maleFreq=0, subPops=[0])
        initSex(pop, maleFreq=1, subPops=[1])
        initSex(pop, maleFreq=1, subPops=[2])
        pop.setVirtualSplitter(SexSplitter())
        migrate(pop, mode=BY_COUNTS,
            rate = [ [0, 50, 50],
                             [0, 0, 0],
                             [0, 0, 0] ],
            subPops = [(0, 0), 1, 2])
        self.assertEqual(pop.subPopSizes(), (2000, 4000, 4000))
        migrate(pop, mode=BY_COUNTS,
            rate = [ [0, 50, 50],
                             [0, 0, 0],
                             [0, 0, 0] ],
            subPops = [(0, 1), 1, 2])
        self.assertEqual(pop.subPopSizes(), (1900, 4050, 4050))
        self.assertEqual(pop.subPopSize([0, 0]), 0)
        self.assertEqual(pop.subPopSize([0, 1]), 1900)
        self.assertEqual(pop.subPopSize([1, 0]), 4000)
        self.assertEqual(pop.subPopSize([1, 1]), 50)
        self.assertEqual(pop.subPopSize([2, 0]), 4000)
        self.assertEqual(pop.subPopSize([2, 1]), 50)
        migrate(pop, mode=BY_COUNTS,
            rate = [ [0, 25, 25],
                     [0, 25, 25],
                             [50, 0, 0],
                             [50, 0, 0] ],
            subPops=[(0,0), (0,1), 1, 2])
        # 25F, 0 -> 1
        # 25F, 0 -> 2
        # 0M migrated.
        # 50, 1 -> 0
        # 50, 2 -> 0
        self.assertEqual(pop.subPopSizes(), (1950, 4025, 4025))


    def testmigrateBySexAndProportion(self):
        'Testing migrate by sex and proportion'
        pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
        initSex(pop, maleFreq=0, subPops=[0])
        initSex(pop, maleFreq=1, subPops=[1])
        initSex(pop, maleFreq=1, subPops=[2])
        pop.setVirtualSplitter(SexSplitter())
        # now if we want to inject a mutation whenever fixation happens
        migrate(pop, mode=BY_PROPORTION,
            rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ],
            subPops = [(0,0), (1,0), (2,0)])
        # no one move out from sp 1, move in 200
        # 200 male, 2000 female
        # 3900 male
        # 3900 male
        self.assertEqual(pop.subPopSizes(), (2200, 3900, 3900))
        stat(pop, numOfMales=True, vars=['numOfMales_sp', 'numOfFemales_sp'])
        self.assertEqual(pop.dvars(0).numOfMales, 200)
        self.assertEqual(pop.dvars(0).numOfFemales, 2000)
        self.assertEqual(pop.dvars(1).numOfMales, 3900)
        self.assertEqual(pop.dvars(1).numOfFemales, 0)
        self.assertEqual(pop.dvars(2).numOfMales, 3900)
        self.assertEqual(pop.dvars(2).numOfFemales, 0)
        # (200M, 2000F), (3900M, 0F), (3900M, 0F),
        # 20M, 200F ==> 1, 2
        # 390F, 1 => 0
        # 390F, 2 => 1
        migrate(pop, mode=BY_PROPORTION,
            rate = [ [0, .1, .1],
                     [0, .1, .1],
                             [0.1, 0, 0],
                             [0, 0.1, 0] ],
            subPops=[(0,0), (0,1), 1, 2]
        )
        # 220 to sp1 and sp2
        # noone goes from sp1 or sp2 to sp0.
        self.assertEqual(pop.subPopSizes(), (2200-440+390, 3900+220-390+390, 3900+220-390))


    def testmigrateBySexAndProbability(self):
        'Testing migrate by sex and probability'
        def migrateSize():
            pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
            initSex(pop, maleFreq=0, subPops=[0])
            initSex(pop, maleFreq=1, subPops=[1])
            initSex(pop, maleFreq=1, subPops=[2])
            pop.setVirtualSplitter(SexSplitter())
            # now if we want to inject a mutation whenever fixation happens
            migrate(pop, mode=BY_PROBABILITY,
                rate = [ [0, .1, .1],
                             [0.1, 0, 0],
                             [0.1, 0, 0]],
                subPops = [(0,0), (1,0), (2,0)])
            return pop.subPopSizes()

        tested = repeated(migrateSize, 100)
        # 2000 female -> 0
        # 4000 male -> 400 to 0
        # 4000 male -> 400 to 0
        self.assertTrue(abs(tested[0] - 2800) < 100) 
        self.assertTrue(abs(tested[1] - 3600) < 100) 
        self.assertTrue(abs(tested[2] - 3600) < 100)
        #
        pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
        initSex(pop, maleFreq=0, subPops=[0])
        initSex(pop, maleFreq=1, subPops=[1])
        initSex(pop, maleFreq=1, subPops=[2])
        pop.setVirtualSplitter(SexSplitter())
        v = pop.subPopSizes()
        vf = pop.subPopSize([0, 1])
        vm = pop.subPopSize([0, 0])
        def migrateSize(pop):
            pop = pop.clone()
            pop.setVirtualSplitter(SexSplitter())
            migrate(pop, mode=BY_PROBABILITY,
                rate = [ [0, 0.1, 0.1],
                         [0.1, 0, 0],
                         [0, 0.1, 0] ],
                subPops = [(0,1), (1,1), (2,1)])
            return pop.subPopSizes()

        tested = repeated(migrateSize, 100, pop)
        # 2000 female, 800 male -> 200 female to each
        # 3600 male no
        # 3600 male no
        self.assertTrue(abs(tested[0] - vm - vf*0.8) < 20, "{} not close to {}".format(tested[0], vm + vf*0.8))
        self.assertTrue(abs(tested[1] - v[1] - vf*0.1) < 20, "{} not close to {}".format(tested[1], v[1] + vf*0.1))
        self.assertTrue(abs(tested[2] - v[2] - vf*0.1) < 20, "{} not close to {}".format(tested[2], v[2] + vf*0.1))


    def testMigrConstAlleleFreq(self):
        'Testing that migration does not change allele frequency'
        pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
        initGenotype(pop, freq=[.2, .8])
        stat(pop, alleleFreq=[0])
        af = pop.dvars().alleleFreq[0][1]    # ~.2
        # migrate and check if allele frequency changes
        migrate(pop, mode=BY_PROBABILITY,
            subPops = [0], toSubPops = [1,2],
            rate = [.05, .05] )
        stat(pop, alleleFreq=[0])
        self.assertEqual(pop.dvars().alleleFreq[0][1], af)
        migrate(pop, mode=BY_PROBABILITY,
            rate = [ [0, .05, .05],
                             [0.025, 0, 0],
                             [0.025, 0, 0] ])
        self.assertEqual(pop.dvars().alleleFreq[0][1], af)
        migrate(pop, mode=BY_PROBABILITY,
            rate = [ [0, .25, .25],
                             [0.25, 0, 0],
                             [0, 0.25, 0] ])
        self.assertEqual(pop.dvars().alleleFreq[0][1], af)


    def testSplitSubPops(self):
        'Testing Population split'
        pop = Population(size=10, loci=[2,6], infoFields=['migrate_to'])
        initGenotype(pop, freq=[.2,.4,.4])
        genotype = list(pop.genotype())
        splitSubPops(pop, subPops=0, sizes=[2, 8], randomize=False)
        # individual untouched
        self.assertEqual(pop.genotype(), genotype)
        # split, with randomization
        splitSubPops(pop, subPops=1, sizes=[6,2], randomize=True)
        self.assertNotEqual(pop.genotype(), genotype)
        # test subpopulation names
        pop = Population(size=10, loci=[2,6], infoFields=['migrate_to'])
        self.assertRaises(ValueError, splitSubPops, pop, subPops=0,
            sizes=[2, 8], names='ab', randomize=False)
        splitSubPops(pop, subPops=0, sizes=[2, 8], names=['ab', 'cd'],
            randomize=False)
        self.assertEqual(pop.subPopName(0), 'ab')
        self.assertEqual(pop.subPopName(1), 'cd')
        splitSubPops(pop, subPops=1, proportions=[0.5, 0.5], randomize=False)
        self.assertEqual(pop.subPopName(2), 'cd')
        # names from noone...
        pop = Population(size=[10, 20], loci=[2,6], infoFields=['migrate_to'])
        splitSubPops(pop, subPops=1, sizes=[12, 8], names=['ab', 'cd'],
            randomize=False)
        self.assertEqual(pop.subPopName(0), '')
        self.assertEqual(pop.subPopName(1), 'ab')
        self.assertEqual(pop.subPopName(2), 'cd')

    def testRearrange(self):
        'Testing if info and genotype are migrated with individuals'
        if moduleInfo()['alleleType'] == 'binary':
            return
        #turnOnDebug('DBG_POPULATION')
        pop = Population(size=[4, 6], loci=[1], infoFields=['a', 'b', 'migrate_to'])
        pop.genotype()[:] = list(range(20))
##         pop.arrIndInfo()[:] = range(20)
##         self.assertEqual([x.allele(0,0) for x in pop.individuals(0)],
##             [int(x) for x in pop.indInfo('a', 0)])
##         self.assertEqual([x.allele(0,0) for x in pop.individuals(1)],
##             [int(x) for x in pop.indInfo('a', 1)])
        #Pymigrate(pop, subPopID=[1,1,1,1,0,0,0,0,0,0])
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
        'Testing Population merge'
        pop = Population(size=[2,4,4], loci=[2,6])
        mergeSubPops(pop, subPops=[0,2])
        self.assertEqual(pop.subPopSizes(), (6,4))
        mergeSubPops(pop, subPops=[0,1])
        self.assertEqual(pop.subPopSizes(), (10,))

    def testResizeSubPops(self):
        'Testing Population resize'
        pop = Population(size=[2, 4, 4], loci=[2,6])
        initGenotype(pop, freq=[0.3, 0.7])
        resizeSubPops(pop, sizes=[6], subPops=[0])
        self.assertEqual(pop.subPopSizes(), (6, 4, 4))
        for ind in (2, 4):
            self.assertEqual(pop.individual(ind, 0), pop.individual(0, 0))
            self.assertNotEqual(pop.individual(ind, 0), pop.individual(1, 0))
        for ind in (3, 5):
            self.assertEqual(pop.individual(ind, 0), pop.individual(1, 0))
            self.assertNotEqual(pop.individual(ind, 0), pop.individual(0, 0))
        # no propagate
        resizeSubPops(pop, sizes=[8, 7], subPops=[1,2], propagate=False)
        self.assertEqual(pop.subPopSizes(), (6, 8, 7))
        for ind in (10, 11, 12, 13, 18, 19, 20):
            self.assertEqual(pop.individual(ind).genotype(), [0]*16)

    def testMigrateByBackwardProportion(self):
        'Testing migrate by proportion'
        pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
        # now if we want to inject a mutation whenever fixation happens
        self.assertRaises(ValueError, migrate, pop, rate=[1, 0.1])
        backwardMigrate(pop, mode=BY_PROPORTION,
            rate =[[ 0.9  ,  0.05 ,  0.05 ],
                   [ 0.025,  0.975,  0.   ],
                   [ 0.025,  0.   ,  0.975]]
        )
        # is not exactly 2000, 4000, 4000 due to numeric issue.
        self.assertEqual(pop.subPopSizes(), (1998, 4001, 4001))
        # 
        pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
        backwardMigrate(pop, mode=BY_PROPORTION,
            rate = [[ 0.5       ,  0.5       ,  0.        ],
       [ 0.11111111,  0.66666667,  0.22222222],
       [ 0.14285714,  0.        ,  0.85714286]])
        self.assertEqual(pop.subPopSizes(), (2002, 4498, 3500))


    def testMigrateByBackwardProbability(self):
        'Testing migrate by probability'
        def migrateSize():
            pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
            # now if we want to inject a mutation whenever fixation happens
            backwardMigrate(pop, mode=BY_PROBABILITY,
                rate = [[ 0.9  ,  0.05 ,  0.05 ],
                       [ 0.025,  0.975,  0.   ],
                       [ 0.025,  0.   ,  0.975]])
            return pop.subPopSizes()
        #
        tested = repeated(migrateSize, times=100)
        # print pop.subPopSizes()
        self.assertTrue(abs(tested[0] - 2000) < 20) 
        self.assertTrue(abs(tested[1] - 4000) < 20) 
        self.assertTrue(abs(tested[2] - 4000) < 20) 
    
        def migrateSize():
            pop = Population(size=[2000,4000,4000], loci=[2], infoFields=['migrate_to'])
            backwardMigrate(pop, mode=BY_PROBABILITY,
                rate = [[ 0.5       ,  0.5       ,  0.        ],
               [ 0.11111111,  0.66666667,  0.22222222],
               [ 0.14285714,  0.        ,  0.85714286]])
            return pop.subPopSizes()

        tested = repeated(migrateSize, times=100)
        # print pop.subPopSizes()
        self.assertTrue(abs(tested[0] - 2000) < 20) 
        self.assertTrue(abs(tested[1] - 4500) < 20)
        self.assertTrue(abs(tested[2] - 3500) < 20)



if __name__ == '__main__':
    unittest.main()
