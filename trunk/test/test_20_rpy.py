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
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', byRep=True, lty=[1, 2, 3],
                type='l', win=10, main='subPop size', xlab='gen', ylim=[0, 100], 
                ylab=['subPopSize (rep %d)' % x for x in range(3)])
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
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', ylab='sp', type='l',
                 byVal=True, win=10, main='subPop size')
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotter(self):
        'Testing byVal paramter of varPlotter'
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
                 win=10, main='subPop size', lty=range(1, 4))
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterByRepVal(self):
        'Testing byVal paramter of varPlotter'
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
             varPlotter('[x**2 for x in subPopSize]', ylab='sp', type='l', byVal=True, byRep=True,
                 win=10, main=['subPop size %d' % x for x in range(9)],
                 valArgs=['main'], lty=range(1, 4))
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
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('[x**2 for x in subPopSize]', byRep=True,
                 win=10, update=5, main='subPop size', ylab='sp',
                 saveAs='demo')
             ],
             gen = 31
        )
        for f in range(5,31,5):
            self.assertEqual(os.path.isfile('demo%d.eps' % f), True)
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
            mode=ByProbability)
        stator = stat(popSize=1, stage=PreMating)
        simu.evolve(
            preOps = [initSex()],
            ops = [
             migr,
             stator,
             varPlotter('subPopSize', byRep=True, ylim=[0,100],
                 ylab='sp', win=10, update=5, main='subPop size')
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

if __name__ == '__main__':
    unittest.main()
