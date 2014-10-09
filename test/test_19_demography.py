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
        # test the removeEmptySubPops parameter
        InstantChangeModel(T=0, N0=[100, 0])._assertSize(
            {
                0: (100, 0),
            }, initSize=[200, 300])
        InstantChangeModel(T=0, N0=[0, 100], removeEmptySubPops=True)._assertSize(
            {
                0: (100,),
            }, initSize=[200, 300])
        InstantChangeModel(T=0, removeEmptySubPops=True)._assertSize(
            {
                0: (200, 300, 400, 500),
            }, initSize=[200, 300, 400, 0, 0, 500])

    def testAdmixtureModel(self):
        AdmixtureModel(T=10, model=('HI', 1, 2, 0.5))._assertSize(
            {
                0: (200, 300, 400, 600),
                4: (200, 300, 400, 600),
                5: (200, 300, 400, 600),
            }, initSize=[200, 300, 400])
        AdmixtureModel(T=10, model=('HI', 0, 2, 0.5))._assertSize(
            {
                0: (200, 300, 400, 400),
                4: (200, 300, 400, 400),
                5: (200, 300, 400, 400),
            }, initSize=[200, 300, 400])
        AdmixtureModel(T=10, model=('CGF', 0, 2, 0.8))._assertSize(
            {
                0: (200, 300, 400),
                4: (200, 300, 400),
                5: (200, 300, 400),
            }, initSize=[200, 300, 400])
        AdmixtureModel(T=10, model=('CGF', 1, 2, 0.8))._assertSize(
            {
                0: (200, 300, 400),
                4: (200, 300, 400),
                5: (200, 300, 400),
            }, initSize=[200, 300, 400])

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
        # carrying capacity
        LinearGrowthModel(T=100, N0=[200, 200], r=0.01, NT=[300, 350])._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                50: [300, 302],
                99: [300, 350],
            }, initSize=500)

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
        # carrying capacity
        ExponentialGrowthModel(T=100, N0=[200, 200], r=0.01, NT=[300, 350])._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                50: [300, 333],
                99: [300, 350],
            }, initSize=500)

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
        MultiStageModel([
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            ExponentialGrowthModel(T=100, N0=[0.4, 0.6], r=0.001),
            ExponentialGrowthModel(r=0.01, NT=[2000, 4000]),
            InstantChangeModel(N0=[1000, 0], T=100)
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
                350: [1000, 0]
            })
        MultiStageModel([
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            ExponentialGrowthModel(T=100, N0=[0.4, 0.6], r=0.001),
            ExponentialGrowthModel(r=0.01, NT=[2000, 4000]),
            AdmixtureModel(model=['HI', 0, 1, 0.3], T=1),
            InstantChangeModel(N0=[0, 0, None], T=100)
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
                311: [2000, 4000, 5714],
                350: [0, 0, 5714]
            })
        MultiStageModel([
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            ExponentialGrowthModel(T=100, N0=[0.4, 0.6], r=0.001),
            ExponentialGrowthModel(r=0.01, NT=[2000, 4000]),
            AdmixtureModel(model=['HI', 0, 1, 0.3], T=1),
            InstantChangeModel(N0=[0, 0, None], removeEmptySubPops=True, T=1)
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
                311: [2000, 4000, 5714],
                312: [5714]
            })
        MultiStageModel([
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            InstantChangeModel(N0=[1000, 1000], T=0),
            ExponentialGrowthModel(T=100, N0=[0.4, 0.6], r=0.001),
        ])._assertSize(
            {
                0: 1010,
                99: 2000,
                100: [400, 600],
            })
        MultiStageModel([
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            InstantChangeModel(N0=[1000, 1000], T=0),
            ExponentialGrowthModel(T=100, N0=[0.4, 0.6], r=0.001),
            InstantChangeModel(N0=[1000, 1000], T=1),
        ])._assertSize(
            {
                0: 1010,
                99: 2000,
                100: [400, 600],
                199: [442, 663],
                200: [1000, 1000],
            })
        self.assertRaises(ValueError, MultiStageModel, [
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            InstantChangeModel(N0=[1000, 1000], T=0),
            ExponentialGrowthModel(T=100, N0=[0.4, 0.6], r=0.001),
            InstantChangeModel(N0=[1000, 1000], T=0),
        ])

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

    
    def testRevertSimulation(self):
        #
        # There is a big problem with applying a demographic model
        # NOT in sequential manner. That is to say, when we apply
        # a demographic model to a population for the first time, it
        # will take its generation number as the initial generation
        # number. 
        pop = Population(100, loci=1)
        m = LinearGrowthModel(N0=[0.5, 0.5], r=[0.1, 0.2], T=100)
        for gen in [0, 10, 5, 7, 6, 8, 20]:
            pop.dvars().gen = gen
            self.assertEqual(m(pop), [50  + 5 * (gen+1), 50 + 10*(gen+1)])
        #
        m = InstantChangeModel(N0=100, T=100, G=[10, 40], NG=[30, 50])
        for gen in [0, 5, 60, 10, 45, 15, 70]:
            # the object has to be initialized with gen = 0
            pop.dvars().gen = gen
            self.assertEqual(list(m(pop)), 
                [{0:100, 5: 100, 10: 30, 15: 30, 45: 50, 60:50, 70:50}[gen]])
        #
        # for multiple stage stuff
        m = MultiStageModel([
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            InstantChangeModel(N0=[1000, 1000], T=1),
            LinearGrowthModel(T=100, N0=[400, 600], r=0.01),
        ])
        for gen in [0, 10, 40, 20, 50, 110, 120, 10, 40, 50]:
            # the object has to be initialized with gen = 0
            pop.dvars().gen = gen
            self.assertEqual(m(pop), 
                {0:[1000 + 10], 50: [1000 + 10*51],
                10:[1000 + 10*11],
                20:[1000 + 10*21],
                40:[1000 + 10*41],
                110:[400 + 4*10, 600+6*10],
                120:[400 + 4*20, 600+6*20],
                }[gen])
        #
        # for multiple stage stuff
        m = MultiStageModel([
            LinearGrowthModel(T=100, N0=1000, r=0.01),  
            InstantChangeModel(N0=[1000, 1000], T=1),
            LinearGrowthModel(T=100, N0=[0.4, 0.6], r=0.01),
        ])
        for gen in [0, 10, 40, 20, 50, 100, 101, 110, 120, 10, 40, 50]:
            # the object has to be initialized with gen = 0
            pop.dvars().gen = gen
            self.assertEqual(m(pop), 
                {0:[1000 + 10], 50: [1000 + 10*51],
                10:[1000 + 10*11],
                20:[1000 + 10*21],
                40:[1000 + 10*41],
                99:[1000 + 10*100],
                100:[1000, 1000],
                101:[400 + 4*1, 600+6*1],
                110:[400 + 4*10, 600+6*10],
                120:[400 + 4*20, 600+6*20],
                }[gen])
        #

    def testEventBasedModel(self):
        'Test event based implementation of demographic models'
        EventBasedModel(T=100, N0=200)._assertSize(
            {
                0: 200,
                99: 200,
            })

    def testExponentialGrowthEvent(self):
        EventBasedModel(T=100, N0=200,
            events=ExponentialGrowthEvent(rates=0.01)
            )._assertSize(
            {
                0: 202,
                2: 206,
                99: 466,
            }, initSize=500)
        EventBasedModel(T=100, N0=[200, 200], 
            events=ExponentialGrowthEvent(rates=0.01)
            )._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                99: [466, 466],
            }, initSize=500)
        # unknown NT
        EventBasedModel(T=100, N0=[200, 200], 
            events=ExponentialGrowthEvent(rates=[0.01, 0.03])
            )._assertSize(
            {
                0: [202, 206],
                2: [206, 218],
                99: [466, 3574],
            }, initSize=500)
      
    def testLinearGrowthEvent(self):
        EventBasedModel(T=100, N0=200,
            events=LinearGrowthEvent(rates=0.01)
            )._assertSize(
            {
                0: 202,
                2: 206,
                99: 400,
            }, initSize=500)
        EventBasedModel(T=100, N0=[200, 200], 
            events=LinearGrowthEvent(rates=0.01)
            )._assertSize(
            {
                0: [202, 202],
                2: [206, 206],
                99: [400, 400],
            }, initSize=500)
        # unknown NT
        EventBasedModel(T=100, N0=[200, 200], 
            events=LinearGrowthEvent(rates=[0.01, 0.03])
            )._assertSize(
            {
                0: [202, 206],
                2: [206, 218],
                99: [400, 800],
            }, initSize=500)

    def testHIAdmixtureEvent(self):
        EventBasedModel(T=10, 
            events=AdmixToNewPopEvent(subPops=[1, 2], proportions=[0.5, 0.5], at=2)
        )._assertSize(
            {
                0: (200, 300, 400),
                2: (200, 300, 400, 600),
                4: (200, 300, 400, 600),
                5: (200, 300, 400, 600),
            }, initSize=[200, 300, 400])
        EventBasedModel(T=10, 
            events=AdmixToNewPopEvent(subPops=[0, 2], proportions=[0.5, 0.5], at=2)
        )._assertSize(
            {
                0: (200, 300, 400),
                2: (200, 300, 400, 400),
                4: (200, 300, 400, 400),
                5: (200, 300, 400, 400),
            }, initSize=[200, 300, 400])
        EventBasedModel(T=10, 
            events=AdmixToNewPopEvent(subPops=[0, 2], proportions=[0.5, 0.1], at=2)
        )._assertSize(
            {
                0: (200, 300, 400),
                2: (200, 300, 400, 400),
                4: (200, 300, 400, 400),
                5: (200, 300, 400, 400),
            }, initSize=[200, 300, 400])
        EventBasedModel(T=10, 
            events=AdmixToNewPopEvent(subPops=[0, 1, 2], proportions=[0.2, 0.3, 0.5], at=2)
        )._assertSize(
            {
                0: (200, 300, 400),
                # take 160, 240, 400
                2: (200, 300, 400, 800),
                4: (200, 300, 400, 800),
                5: (200, 300, 400, 800),
            }, initSize=[200, 300, 400])

    def testCGFAdmixtureEvent(self):
        EventBasedModel(T=10, model=('CGF', 0, 2, 0.8))._assertSize(
            {
                0: (200, 300, 400),
                4: (200, 300, 400),
                5: (200, 300, 400),
            }, initSize=[200, 300, 400])
        EventBasedModel(T=10, model=('CGF', 1, 2, 0.8))._assertSize(
            {
                0: (200, 300, 400),
                4: (200, 300, 400),
                5: (200, 300, 400),
            }, initSize=[200, 300, 400])



    def testRevertAndDemo(self):
        'Test the use of demographic model with RevertIf operator'
        pop = Population(100, loci=1)
        initpop = pop.clone()
        pop.evolve(initOps=InitSex(), matingScheme=RandomMating(), gen=10)
        pop.save('burnin.pop')
        pop1 = initpop.clone()
        gens = pop1.evolve(
            preOps=RevertIf(True, 'burnin.pop', at=0),
            matingScheme=RandomMating(),
            gen=20
        )
        self.assertEqual(gens, 10)
        # in this case, the evolution starts at generation 10
        # and evolve for another 20 generations.
        pop1 = initpop.clone()
        self.assertEqual(pop1.evolve(
            preOps=RevertIf(True, 'burnin.pop', at=0),
            matingScheme=RandomMating(subPopSize=LinearGrowthModel(T=20, r=0.1))
        ), 20)
        #
        pop1 = initpop.clone()
        self.assertEqual(pop1.evolve(
            initOps=InitSex(),
            postOps=RevertIf(True, 'burnin.pop', at=0),
            matingScheme=RandomMating(subPopSize=LinearGrowthModel(T=20, r=0.1))
        ), 10)


    def testStockModels(self):
        'Test stock demographic models'
        OutOfAfricaModel(20000).plot()
        OutOfAfricaModel(20000, scale=5).plot()
        OutOfAfricaModel(20000, outcome='EU', scale=5).plot()
        OutOfAfricaModel(20000, outcome=['EU', 'AS'], scale=5).plot()
        SettlementOfNewWorldModel(20000).plot()
        SettlementOfNewWorldModel(20000, scale=5).plot()
        SettlementOfNewWorldModel(20000, outcome='MXL', scale=5).plot()
        SettlementOfNewWorldModel(20000, outcome='AF', scale=5).plot()
        SettlementOfNewWorldModel(20000, outcome=['EU', 'AS'], scale=5).plot()
        CosiModel(20000).plot()
        CosiModel(20000, scale=5).plot()


if __name__ == '__main__':
    unittest.main()
