#!/usr/bin/env python
#
# Testing demographic models
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 149 $
# $LastChangedDate: 2006-02-03 15:51:04 -0600 (Fri, 03 Feb 2006) $
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
from time import sleep

from simuPOP.demography import *

class TestDemography(unittest.TestCase):
    def testInstantChangeModel(self):
        InstantChangeModel(T=100, N0=200)._assertSize(
            {
                0: 200,
                99: 200,
            })
        InstantChangeModel(T=100, N0=[200, 300])._assertSize(
            {
                0: [200, 300],
                99: [200, 300]
            })
        InstantChangeModel(T=100, N0=200, G=10, NG=30)._assertSize(
            {
                9: 200,
                10: 30,
                99: 30
            })
        InstantChangeModel(T=100, N0=200, G=10, NG=30)._assertSize(
            {
                0: 200,
                9: 200,
                10: 30,
                99: 30
            }, initSize=500)
        InstantChangeModel(T=100, N0=200, G=10, NG=[30, 40])._assertSize(
            {
                0: 200,
                9: 200,
                10: (30, 40),
                99: (30, 40),
            }, initSize=500)
        InstantChangeModel(T=100, N0=[200, 300], G=10, NG=[30, 40])._assertSize(
            {
                0: (200, 300),
                9: (200, 300),
                10: (30, 40),
                99: (30, 40),
            }, initSize=500)
        InstantChangeModel(T=100, N0=[200, 300], G=10, NG=[[30, 40]])._assertSize(
            {
                0: (200, 300),
                9: (200, 300),
                10: (30, 40),
                99: (30, 40),
            }, initSize=500)
        self.assertRaises(ValueError, InstantChangeModel, T=100, N0=[200, 300],
            G=[10], NG=[30, 40])
        InstantChangeModel(T=100, N0=[200, 300], G=(10, 20), NG=[30, 40])._assertSize(
            {
                0: (200, 300),
                9: (200, 300),
                10: 30,
                20: 40,
                99: 40,
            }, initSize=500)
        InstantChangeModel(T=100, N0=[200, 300], G=(10, 20), NG=[30, [40, 50]])._assertSize(
            {
                0: (200, 300),
                9: (200, 300),
                10: 30,
                20: (40, 50),
                99: (40, 50)
            })
        self.assertRaises(ValueError, InstantChangeModel, T=100, N0=[200, 300], G=(10, 20),
            NG=[30, [40, 50], 40], )
        # testing the use of proportions
        InstantChangeModel(T=100, N0=[200, 300], G=(10, 20), NG=[0.5, [0.4, 0.5]])._assertSize(
            {
                0: (200, 300),
                9: (200, 300),
                10: 250,
                20: (100, 125),
                99: (100, 125)
            })
        InstantChangeModel(T=100, G=(10, 20), NG=[0.5, [0.4, 0.5]])._assertSize(
            {
                0: (200, 300),
                9: (200, 300),
                10: 250,
                20: (100, 125),
                99: (100, 125)
            }, initSize=[200, 300])
        InstantChangeModel(T=100, G=(10, 20), NG=[0.5, [0.4, None]])._assertSize(
            {
                0: (200, 300),
                9: (200, 300),
                10: 250,
                20: (100, 250),
                99: (100, 250)
            }, initSize=[200, 300])
        InstantChangeModel(T=100, G=(10, 20), NG=[0.5, [None, 50]])._assertSize(
            {
                0: (200, 300),
                9: (200, 300),
                10: 250,
                20: (250, 50),
                99: (250, 50)
            }, initSize=[200, 300])

    def testLinearGrowthModel(self):
        self.assertRaises(LinearGrowthModel)
        self.assertRaises(LinearGrowthModel, T=100, N0=200)
        #
        LinearGrowthModel(T=100, N0=200, NT=400)._assertSize(
            {
                0: 202,
                1: 204,
                99: 400,
            }, initSize=500)
        LinearGrowthModel(T=100, N0=200, r=0.01)._assertSize(
            {
                0: 202,
                2: 206,
                99: 400,
            }, initSize=500)
        LinearGrowthModel(N0=200, r=0.01, NT=400)._assertSize(
            {
                0: 202,
                2: 206,
                99: 400,
            }, initSize=500)
        #
        LinearGrowthModel(T=100, N0=[200, 200], NT=[400, 800])._assertSize(
            {
                0: [202, 206],
                1: [204, 212],
                99: [400, 800],
            }, initSize=500)
        LinearGrowthModel(T=100, N0=[200, 200], r=0.01)._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                99: [400, 400],
            }, initSize=500)
        # unknown NT
        LinearGrowthModel(T=100, N0=[200, 200], r=[0.01, 0.03])._assertSize(
            {
                0: [202, 206],
                2: [206, 218],
                99: [400, 800],
            }, initSize=500)
        # unknown T
        LinearGrowthModel(N0=[200, 200], r=0.01, NT=[400, 800])._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                99: [400, 400],
                100: [400, 402],
                298: [400, 798],
                300: [400, 800],
            }, initSize=500)
        # N0 size ...
        LinearGrowthModel(N0=[1., 1.], r=0.01, NT=[400, 800])._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                99: [400, 400],
                100: [400, 402],
                298: [400, 798],
                300: [400, 800],
            }, initSize=200)
        LinearGrowthModel(N0=[(None, 'name'), 1.], r=0.01, NT=[400, 800])._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                99: [400, 400],
                100: [400, 402],
                298: [400, 798],
                300: [400, 800],
            }, initSize=200)
        self.assertRaises(LinearGrowthModel, N0=[1., 1.], r=0.01, NT=[2, 3.])
        self.assertRaises(LinearGrowthModel, N0=[1., 1.], r=0.01, NT=[2000])
        self.assertRaises(LinearGrowthModel, N0=[1., 1.], r=[0.01], NT=[2000, 3000])

    def testExponentialGrowthModel(self):
        self.assertRaises(ExponentialGrowthModel)
        self.assertRaises(ExponentialGrowthModel, T=100, N0=200)
        #
        ExponentialGrowthModel(T=100, N0=200, NT=400)._assertSize(
            {
                0: 201,
                1: 202,
                99: 400,
            }, initSize=500)
        ExponentialGrowthModel(T=100, N0=200, r=0.01)._assertSize(
            {
                0: 202,
                2: 206,
                99: 540,
            }, initSize=500)
        ExponentialGrowthModel(N0=200, r=0.01, NT=400)._assertSize(
            {
                0: 202,
                2: 206,
                99: 540,
            }, initSize=500)
        #
        ExponentialGrowthModel(T=100, N0=[200, 200], NT=[400, 800])._assertSize(
            {
                0: [201, 202],
                1: [202, 205],
                99: [400, 800],
            }, initSize=500)
        ExponentialGrowthModel(T=100, N0=[200, 200], r=0.01)._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                99: [540, 540],
            }, initSize=500)
        # unknown NT
        ExponentialGrowthModel(T=100, N0=[200, 200], r=[0.01, 0.03])._assertSize(
            {
                0: [202, 206],
                2: [206, 218],
                99: [543, 4017],
            }, initSize=500)
        # unknown T
        ExponentialGrowthModel(N0=[200, 200], r=0.01, NT=[400, 800])._assertSize(
            {
                2: [206, 206],
                99: [400, 543],
                100: [400, 549],
            }, initSize=500)
        # N0 size ...
        ExponentialGrowthModel(N0=[1., 1.], r=0.01, NT=[400, 800])._assertSize(
            {
                2: [206, 206],
                99: [400, 543],
                100: [400, 549],
            }, initSize=200)
        ExponentialGrowthModel(N0=[(None, 'name'), 1.], r=0.01, NT=[400, 800])._assertSize(
            {
                2: [206, 206],
                99: [400, 543],
                100: [400, 549],
            }, initSize=200)
        self.assertRaises(ExponentialGrowthModel, N0=[1., 1.], r=0.01, NT=[2, 3.])
        self.assertRaises(ExponentialGrowthModel, N0=[1., 1.], r=0.01, NT=[2000])
        self.assertRaises(ExponentialGrowthModel, N0=[1., 1.], r=[0.01], NT=[2000, 3000])

    def testMultiStageModel(self):
        MultiStageModel([
            InstantChangeModel(T=1000, N0=1000, G=[500, 600], NG=[100, 1000]),  
            ExponentialGrowthModel(T=100, NT=10000)
        ])._assertSize(
            {
                100: 1000,
                500: 100,
                501: 100,
                600: 1000,
                999: 1000,
                1098: 9772,
                1099: 10000,
            })
        MultiStageModel([
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            ExponentialGrowthModel(T=100, N0=[0.4, 0.6], r=0.001),
            ExponentialGrowthModel(r=0.01, NT=[2000, 4000])
        ])._assertSize(
            {
                0: 1010,
                99: 2000,
                100: [800, 1201],
                102: [802, 1203],
                199: [884, 1326],
                250: [1472, 2208],
                300: [2000, 3640],
                310: [2000, 4000],
            })


    def testTerminator(self):
        model = MultiStageModel([
            InstantChangeModel(N0=1000, 
                ops=[
                    Stat(alleleFreq=ALL_AVAIL, numOfSegSites=ALL_AVAIL),
                    # terminate if the average allele frequency of segregating sites
                    # are more than 0.1 
                    TerminateIf('sum([x[1] for x in alleleFreq.values() if '
                        'x[1] != 0])/(1 if numOfSegSites==0 else numOfSegSites) > 0.1')
                ]
            ),
            ExponentialGrowthModel(N0=[0.5, 0.5], r=0.01, NT=[2000, 5000])
            ]
        )
        pop = Population(size=model.init_size, loci=100)
        pop.evolve(
            initOps=InitSex(),
            preOps=SNPMutator(u=0.01, v=0.01),
            matingScheme=RandomMating(subPopSize=model),
        )

    
if __name__ == '__main__':
    unittest.main()
