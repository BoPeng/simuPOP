#!/usr/bin/env python
#
# Testing plotting with R
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 149 $
# $LastChangedDate: 2006-02-03 15:51:04 -0600 (Fri, 03 Feb 2006) $
#

import simuOpt
simuOpt.setOptions(quiet=False)

from simuPOP import *
import unittest, sys, os, exceptions

from time import sleep

hasRPy = True
try:
    from simuRPy import *
    from rpy import r
except exceptions.ImportError:
    print "simuRPy can not be imported. Either rpy or r is not installed properly."
    hasRPy = False

class TestRPy(unittest.TestCase):
    def testVarPlotterBase(self):
        'Testing byRep parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize[0]',
                type='l', win=10, main='Three colorful lines, no legend, win=10',
                xlab='gen', ylim=[0, 100], 
                col_rep=['red', 'green', 'blue'])
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testVarPlotterByRep(self):
        'Testing byRep parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', byRep=True, lty_dim=[1, 2, 3],
                type='l', win=10, main='3 rep, 3 colorful lines, ylabs differ', xlab='gen', ylim=[0, 100], 
                col_dim=['red', 'green', 'blue'], 
                ylab_rep=['subPopSize (rep %d)' % x for x in range(3)])
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testVarPlotterByDim(self):
        'Testing byDim paramter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', ylab='sp', type='l',
                 col_rep=['red', 'green', 'blue'],
                 byDim=True, win=10, main='win=10, 3 dim, 3 colorful lines, legend',
                 legend=['a', 'b', 'c'])
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testVarPlotterSelectedRep(self):
        'Testing byDim paramter of varPlotter using selected replicates'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=5)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', ylab='sp', type='l',
                 col_rep=['red', 'green', 'blue'], byDim=True, 
                 main='3 out of 5 reps, 3 dim plots', rep=[0, 2, 3])
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testVarPlotterTogether(self):
        'Testing byDim paramter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', ylab='sp',
                 col_rep=['red', 'green', 'blue'],
                 main='9 lines, col rep, lty dim',
                 lty_dim=range(1, 4), ylim=[0, 10000],
                 legend=['rep1-sp1', 'rep1-sp2', 'rep1-sp3', 'rep2-sp1',
                    'rep2-sp2', 'rep2-sp3', 'rep3-sp1', 'rep3-sp2', 'rep3-sp3'])
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testVarPlotterByRepDim(self):
        'Testing byDim paramter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', ylab='sp', 
                 byDim=True, byRep=True,
                 win=10, main_rep_dim=['rep dim %d' % x for x in range(9)],
                 col_rep=['red', 'green', 'blue'], lty_dim=range(1, 4))
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testVarPlotterSaveAs(self):
        'Testing saveAs parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', byRep=True,
                 win=10, update=5, ylim=[0, 10000],
                 main='Save as, 3 rep, no color', ylab='sp',
                 saveAs='demo')
             ],
             gen = 31
        )
        for f in range(5,31,5):
            self.assertEqual(os.path.isfile('demo%d.eps' % f), True)
            os.remove('demo%d.eps' % f)
        sleep(5)
        r.dev_off()

    def testVarPlotterPar(self):
        'Testing ylim parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', byRep=True, ylim=[0,100],
                 ylab='sp', win=10, update=5, par_mfrow=[1, 3], main='mfrow=[1, 3]')
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testVarPlotterPar(self):
        'Testing ylim parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', byRep=True, ylim=[0,100],
                 ylab='sp', win=10, update=5, par_mfrow=[1, 3], main='mfrow=[1, 3]')
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testVarPlotterHook(self):
        'Testing ylim parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        def setPar(r):
            r.par(mar=[2]*4)
        def mtext(r):
            r.mtext('A marginal text', outer=True, side=3)
        def drawFrame(r, rep, dim):
            r.axis(1)
            r.axis(2)
            r.grid()

        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', byRep=True, ylim=[0,100], plot_axes=False,
                 ylab='sp', col_dim=['red', 'blue', 'black'],
                 update=5, main_rep=['dimension %d' % x for x in range(3)],
                 preHook=setPar,
                 postHook=mtext,
                 plotHook=drawFrame)
             ],
             gen = 30
        )
        sleep(5)
        r.dev_off()

    def testInfoPlotter(self):
        import random
        pop = population([100, 200], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            ops = [
                inheritTagger(TAG_Paternal, infoFields=['x']),
                inheritTagger(TAG_Maternal, infoFields=['y']),
                infoPlotter(['x', 'y'], main='B/W, 300 points'),
                #pause(),
            ],
            gen = 5,
        )
        sleep(5)
        r.dev_off()

    def testInfoPlotterSP(self):
        import random
        pop = population([100, 200], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            ops = [
                inheritTagger(TAG_Paternal, infoFields=['x']),
                inheritTagger(TAG_Maternal, infoFields=['y']),
                infoPlotter(['x', 'y'],
                    subPops = [(0, 0), (0, 1), (1, 0), (1, 1)],
                    col_sp = ['blue', 'red', 'green', 'purple'],
                    main='Color, 300 points, left right does not mix'),
                #pause(),
            ],
            gen = 5,
        )
        sleep(5)
        r.dev_off()

    def testInfoPlotterSubSet(self):
        import random
        pop = population([100, 200], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            ops = [
                inheritTagger(TAG_Paternal, infoFields=['x']),
                inheritTagger(TAG_Maternal, infoFields=['y']),
                infoPlotter(['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    col_sp = ['blue', 'red'],
                    pch_sp = [1, 2],
                    xlim = [0, 1],
                    main='Twoo colors, 100 points xlim=[0, 1]',
                    legend = ['Male', 'Female']),
                pause(),
            ],
            gen = 5,
        )
        sleep(5)
        r.dev_off()

if __name__ == '__main__':
    unittest.main()
