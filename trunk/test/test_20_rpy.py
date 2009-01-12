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
except exceptions.ImportError:
    print "simuRPy can not be imported. Either rpy or r is not installed properly."
    hasRPy = False

class TestRPy(unittest.TestCase):

    def testVarPlotterByRep(self):
        'Testing byRep parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=MigrByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', numRep=3, byRep=True,
                 varDim=3, win=10, title='subPop size')
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterByVal(self):
        'Testing byVal paramter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=MigrByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]',
                 numRep=3, byVal=True,
                 varDim=3, win=10, title='subPop size')
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterSeparate(self):
        'Testing separate parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=MigrByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', numRep=3,
                 byVal=True,
                 varDim=3, win=10,
                 title='subPop size',
                 separate=1)
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterSaveAs(self):
        'Testing saveAs parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=MigrByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', numRep=3, byRep=True,
                 varDim=3, win=10, update=5, title='subPop size',
                 saveAs='demo')
             ],
             gen = 31
        )
        for f in range(5,31,5):
            os.remove('demo%d.eps' % f)
        sleep(1)
        r.dev_off()

    def testVarPlotterYlim(self):
        'Testing ylim parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=MigrByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', numRep=3, byRep=True, ylim=[0,100],
                 separate=True, varDim=3, win=10, update=5, title='subPop size')
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterNoHistory(self):
        'Testing no history parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=MigrByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]',
                 history=False,
                 numRep=3, byRep=True,
                 win=10, update=5, title='subPop size')
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterNoHisSeparate(self):
        'Testing separate parameter in the no history case'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=MigrByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', history=False,
                 numRep=3, win=10, update=5, separate=1,
                 title='subPop size')
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterImage(self):
        'Testing the image parameter of varPlotter'
        if not hasRPy:
            return True
        # test the operator
        nr = 1
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=nr)
        migr = migrator([[0,.2,.1],[.25,0,.1],[.1,.2,0]],
                                        mode=MigrByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.setGen(0)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', history=True, numRep=nr, byRep=True,
                 win=10, update=5, plotType='image', varDim=3, title='subPop size')
             ],
             gen = 10
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterImageNoHistory(self):
        'Testing no history version of image parameter'
        if not hasRPy:
            return True
        # testing histroy = false and plotType='image'
        numRep=4
        popSize=500
        #turnOnDebug(DBG_MUTATOR)
        #turnOnDebug(DBG_SIMULATOR)
        simu = simulator(population(size=popSize, loci=[10], infoFields = ['migrate_to']),
            randomMating(), rep=numRep)
        simu.evolve(
            preOps = [ initByFreq([.2,.3,.5])],
            ops = [
                recombinator(rate=0.01),
                # count number of affected
                stat(LD=[ [x, y] for x in range(0,10) for y in range(0,10)]),
                # plot affected status
                varPlotter(expr='LD', plotType='image', numRep=4, byRep=True, update=5,
                    title='pairwise LD', history=False),
                # pause(rep=REP_LAST, step=5),
            ],
            gen = 20,
        )
        sleep(1)
        r.dev_off()

if __name__ == '__main__':
    unittest.main()
