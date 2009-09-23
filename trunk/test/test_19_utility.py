#!/usr/bin/env python
#
# Purpose:
#    testing of interfaces of random number selector
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#


import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
from simuUtil import *

import unittest, os, sys, exceptions
import random

# for memory leak testing.

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]


def memory(since=0.0):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:') - since

def resident(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since


def stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since


class TestUtility(unittest.TestCase):

    def interactiveTestPauseAtGen(self):
        'Testing resume to simulation'
        simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
            randomMating(), reps=5)
        print "\n\nUSER INTERACTION: Please press q\n\n"
        self.assertRaises( exceptions.SystemError, simu.evolve,
            ops=[ pause(at=[10]),
                        # should quite, can not reach generation 12
                        terminateIf("True", at=[12] ) ] )

    def interactiveTestExitToShell(self):
        'Testing exit to a shell'
        simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
            randomMating(), reps=5)
        print "\n\nUSER INTERACTION: Please press s and then Ctrl-D"
        print "Please check the existence of variable pop\n\n"
        simu.evolve(
            ops=[ pause(at=[10]) ], end=12)
        print "\n\nUSER INTERACTION: Please press s and then Ctrl-D"
        print "Please check the existence of variable tmpPop\n\n"
        simu.evolve(
            ops=[ pause(at=[20], popName='tmpPop') ], end=25)


    def testSetRNG(self):
        'Testing all RNG types'
        for rg in AvailableRNGs():
            SetRNG(rg)

    def testDefaultRNG(self):
        'Testing default RNG'
        rg = GetRNG()
        self.assertEqual(rg.name(), 'mt19937')

    def testBinomial(self):
        'Testing binomial distribution'
        rg = GetRNG()
        for n in range(1,10):
            rg.randBinomial(10, .7)

    def testUniform01(self):
        'Testing uniform distribution generator'
        rg = GetRNG()
        for n in range(1,10):
            rg.randUniform01()

    def testRandBit(self):
        'Testing random bit function'
        rg = GetRNG()
        sum = 0
        for n in range(10000):
            sum += rg.randBit()
        self.assertEqual(sum < 5100, True)
        self.assertEqual(sum > 4900, True)


    def testBernulliTrials(self):
        'Testing bernullitrials'
        import math
        rg = GetRNG()
        p = [0.00001, 0.001, 0.5, 0.99]
        N = 1000000
        bt = BernulliTrials(rg, p, N)
        bt.doTrial()
        for i in range(len(p)):
            prop = bt.trialSuccRate(i)
            # binomial, mean p, variance = p(1-p)/n
            std = math.sqrt(p[i]*(1.-p[i])/N)
            assert prop > p[i] - 3*std and prop < p[i] + 3*std
        # another test, for each trail
        for pi in p:
            bt = BernulliTrials(rg, [pi]*N, 10)
            bt.doTrial()
            for i in range(10):
                bt.trial();
                prop = bt.probSuccRate()
                # binomial, mean p, variance = p(1-p)/n
                std = math.sqrt(pi*(1.-pi)/N)
                assert prop > pi - 3*std and prop < pi + 3*std, \
                    'Obtain %f, expected %f (std: %f)' % (prop, pi, std)

        # test find_first and find_next
        bt = BernulliTrials(rg, p, N)
        bt.doTrial()
        for i in range(len(p)):
            pos = bt.trialFirstSucc(i)
            assert pos < N
            if pos != bt.npos:
                assert bt.trialSucc(i, pos)
                for j in range(pos):
                    assert not bt.trialSucc(i, j)
                while True:
                    last_pos = pos
                    pos = bt.trialNextSucc(i, pos)
                    if pos == bt.npos:
                        break;
                    assert pos < N
                    assert bt.trialSucc(i, pos)
                    for j in range(last_pos+1, pos):
                        assert not bt.trialSucc(i, j)
                for j in range(pos+1, N):
                    assert not bt.trialSucc(i, j)

    def testSeed(self):
        'Testing RNG::seed() and RNG::setSeed()'
        import random
        seed = []
        name = AvailableRNGs()[0]
        for i in range(100):
            SetRNG(name)
            sd = GetRNG().seed()
            assert sd not in seed
            seed.append(sd)
        # test set seed
        sd = random.randint(100, 10000)
        SetRNG(name, sd)
        self.assertEqual(GetRNG().seed(), sd)
        # test if sequences are the same once the seed is set
        sd = random.randint(100, 10000)
        SetRNG(name, sd)
        seq = [GetRNG().randInt(10000) for x in range(100)]
        SetRNG(name, sd)
        seq1 = [GetRNG().randInt(10000) for x in range(100)]
        self.assertEqual(seq, seq1)
        # randBit need to be treated separately because it uses
        # global variable of RNG().
        sd = random.randint(100, 10000)
        SetRNG(name, sd)
        seq = [GetRNG().randBit() for x in range(100)]
        SetRNG(name, sd)
        seq1 = [GetRNG().randBit() for x in range(100)]
        self.assertEqual(seq, seq1)

    def testWeightedSampler(self):
        'Testing weighted sampler'
        sampler = weightedSampler(GetRNG(), [1, 2, 3, 4])
        num = []
        for i in range(100000):
            num.append(sampler.get())
        #print [num.count(i) for i in range(4)]


    def TestLargePedigree(self):
        'Testing getting large pedigree, for simuUtils.ascertainPedigree'
        import simuUtil
        pop = population(100, ancestralDepth=2, infoFields=['father_idx', 'mother_idx'])
        simu = simulator(pop, randomMating(numOffspring=0.3, mode=MATE_GeometricDistribution))
        simu.evolve(ops=[parentsTagger()], end=5)
        pop = simu.population(0)
        #
        def comb(geno):
            return sum(geno)+1
        SampleLargePedigree(pop, numPedigree=10, minPedSize=5, minAffected=0,
            maxOffspring=5, output='ped', loci=[], combine=comb)
        simuUtil.VC_merlin('ped')


    def testMemoryLeakLoadPopulation(self):
        'Testing if loadPopulation leaks memory'
        # run this and see if memory usage goes up continuously
        pop = population(100, loci=[400])
        InitByFreq(pop, [0.2, 0.8])
        Stat(pop, alleleFreq=range(pop.totNumLoci()))
        pop.save('test.bin')
        for i in range(4):
            pop = LoadPopulation('test.bin')
            #pop = population(100, loci=[1000]*10)
            Stat(pop, alleleFreq=range(pop.totNumLoci()))
            #print pop.dvars().alleleFreq[100][0]
            if i < 2:
                # let m0 stablize
                m0 = memory()
                #print 'M0', m0
            else:
                m1 = memory(m0)
                #dprint 'M1', m1
                self.assertEqual(m1, 0.0)
            del pop
        os.remove('test.bin')

    def testTrajectorySimulatorSimple(self):
        'Testing trajectory simulation with constant population size and no selection'
        import logging
        logging.basicConfig(level=logging.DEBUG)
        logger = logging.getLogger()
        trajSimulator = trajectorySimulator(N=1000)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1)
        self.assertNotEqual(traj, None)
        # trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0.])
        self.assertAlmostEqual(traj.freq(3000, 0)[0], 0.1)
        # invalid inputs
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=3000, freqEnd=[0.1, 0.2])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=0, freqEnd=0.1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=400,
            maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500 and traj._beginGen() <= 2600)
        # one mutation
        self.assertEqual(len(traj.mutators()), 1)
        # the function form
        traj = BackwardTrajectory(N=1000, genEnd=3000, freqEnd=0.1) 
        #
        # simuForward
        # invalid input
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=0.09)
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.09])
        # test trajectory
        traj = trajSimulator.simuForward(genBegin=100, genEnd=200,
            freqBegin=0.1, freqEnd=[0.09, 0.11])
        self.assertAlmostEqual(traj.freq(100, 0)[0], 0.1)
        self.assertTrue(traj.freq(200, 0)[0] >= 0.09 and traj.freq(200, 0)[0] <= 0.11)
        # no mutation
        self.assertEqual(len(traj.mutators()), 0)
        # the function form
        traj = ForwardTrajectory(N=1000, genBegin=100, genEnd=200, freqBegin=0.1,
            freqEnd=[0.09, 0.11]) 
        #
        # simple with subpopulation
        trajSimulator = trajectorySimulator(N=[1000, 2000])
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1)
        self.assertNotEqual(traj, None)
        # trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0.])
        self.assertAlmostEqual(traj.freq(3000, 0)[0], 0.1)
        # invalid inputs
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=3000, freqEnd=[0.1, 0.2, 0.3])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=0, freqEnd=0.1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=[0.1, 0.2], minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=[0.1, 0.2], maxMutAge=1000)
        self.assertTrue(traj._beginGen() >= 2000)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=[0.1, 0.2], minMutAge=400,
            maxMutAge=1000)
        self.assertTrue(traj._beginGen() >= 2000 and traj._beginGen() <= 2600)
        # one mutation in each subpopulation
        self.assertEqual(len(traj.mutators()), 2)
        # the function form
        traj = BackwardTrajectory(N=[1000, 2000], genEnd=3000, freqEnd=0.1) 
        #
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=0.09)
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.09])
        # test trajectory
        traj = trajSimulator.simuForward(genBegin=100, genEnd=200,
            freqBegin=0.1, freqEnd=[0.09, 0.11])
        self.assertAlmostEqual(traj.freq(100, 0)[0], 0.1)
        fq = (traj.freq(200, 0)[0] + traj.freq(200, 1)[0]*2)/3.
        self.assertTrue(fq >= 0.09 and fq <= 0.11)
        # no mutation
        self.assertEqual(len(traj.mutators()), 0)
        # the function form
        traj = ForwardTrajectory(N=[1000, 2000], genBegin=100, genEnd=200, freqBegin=0.1,
            freqEnd=[0.09, 0.11]) 


    def testTrajectorySimulatorSimpleMultiLoci(self):
        'Testing trajectory simulation of multiple loci with constant population size and no selection'
        trajSimulator = trajectorySimulator(N=1000, nLoci=2)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1)
        self.assertNotEqual(traj, None)
        # trajectory ends at 0.
        for gen in range(3000):
            self.assertEqual(len(traj.freq(gen, 0)), 2)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0., 0.])
        # most likely the starting generation will be different (there is a zero...)
        self.assertEqual(traj.freq(traj._beginGen() + 1, 0).count(0.), 1)
        # invalid inputs
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=3000, freqEnd=[0.1, 0.2, 0.3])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=0, freqEnd=0.1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=400,
            maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500 and traj._beginGen() <= 2600)
        # one mutation
        self.assertEqual(len(traj.mutators()), 2)
        # the function form
        traj = BackwardTrajectory(N=1000, nLoci=2, genEnd=3000, freqEnd=0.1) 
        #
        # simuForward
        # invalid input
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=0.09)
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.09])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.08, 0.09])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[[0.18, 0.09]*2])
        # test trajectory
        traj = trajSimulator.simuForward(genBegin=100, genEnd=200,
            freqBegin=[0.1, 0.13], freqEnd=[[0.09, 0.11], [0.11, 0.12]])
        self.assertAlmostEqual(traj.freq(100, 0)[0], 0.1)
        self.assertAlmostEqual(traj.freq(100, 0)[1], 0.13)
        self.assertTrue(traj.freq(200, 0)[0] >= 0.09 and traj.freq(200, 0)[0] <= 0.11)
        self.assertTrue(traj.freq(200, 0)[1] >= 0.11 and traj.freq(200, 0)[0] <= 0.12)
        # no mutation
        self.assertEqual(len(traj.mutators()), 0)
        # the function form
        traj = ForwardTrajectory(N=1000, nLoci=2, genBegin=100, genEnd=200, freqBegin=0.1,
            freqEnd=[[0.09, 0.11]]*2) 
        #
        # cases with multiple subpopulation
        #
        trajSimulator = trajectorySimulator(N=[1000, 2000], nLoci=2)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1)
        self.assertNotEqual(traj, None)
        # trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0., 0.])
        # most likely the starting generation will be different (there is three zeros...)
        self.assertTrue(traj.freq(traj._beginGen() + 1, 0).count(0.) > 1)
        # invalid inputs
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=3000, freqEnd=[0.1, 0.2, 0.3])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=0, freqEnd=0.1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, maxMutAge=1500)
        self.assertTrue(traj._beginGen() >= 1500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=400,
            maxMutAge=1500)
        self.assertTrue(traj._beginGen() >= 1500 and traj._beginGen() <= 2600 )
        # one mutation
        self.assertEqual(len(traj.mutators()), 4)
        # the function form
        traj = BackwardTrajectory(N=[1000, 3000], nLoci=2, genEnd=3000, freqEnd=[0.1, 0.2]) 
        #
        # simuForward
        # invalid input
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=0.09)
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.09])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.08, 0.09])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[[0.18, 0.09]*2])
        # test trajectory
        traj = trajSimulator.simuForward(genBegin=100, genEnd=200,
            freqBegin=[0.1, 0.13], freqEnd=[[0.09, 0.13], [0.10, 0.14]])
        self.assertAlmostEqual((traj.freq(100, 0)[0] + traj.freq(100, 1)[0]*3)/4., 0.1)
        self.assertAlmostEqual((traj.freq(100, 0)[1] + traj.freq(100, 1)[1]*3)/4., 0.13)
        f1 = (traj.freq(200, 0)[0] + traj.freq(200, 1)[0]*3)/4.
        f2 = (traj.freq(200, 0)[1] + traj.freq(200, 1)[1]*3)/4.
        self.assertTrue(f1 >= 0.09 and f1 <= 0.13)
        self.assertTrue(f2 >= 0.10 and f2 <= 0.14)
        # no mutation
        self.assertEqual(len(traj.mutators()), 0)
        # the function form
        traj = ForwardTrajectory(N=[1000, 3000], nLoci=2, genBegin=100, genEnd=200, freqBegin=0.1,
            freqEnd=[[0.09, 0.11]]*2) 
        #


    def testTrajectorySimulatorSel(self):
        'Testing trajectory simulation with constant population size and natural selection'
        # invalid input
        self.assertRaises(exceptions.ValueError, trajectorySimulator, N=1000, fitness=[1, 1.01])
        trajSimulator = trajectorySimulator(N=1000, fitness=[1, 1.01, 1.03])
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1)
        self.assertNotEqual(traj, None)
        # trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0.])
        self.assertAlmostEqual(traj.freq(3000, 0)[0], 0.1)
        # invalid inputs
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=3000, freqEnd=[0.1, 0.2])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=0, freqEnd=0.1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=400,
            maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500 and traj._beginGen() <= 2600)
        # one mutation
        self.assertEqual(len(traj.mutators()), 1)
        # the function form
        traj = BackwardTrajectory(N=1000, fitness=[1, 1.01, 1.03], genEnd=3000, freqEnd=0.1) 
        #
        # simuForward
        # invalid input
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=0.09)
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.09])
        # test trajectory
        traj = trajSimulator.simuForward(genBegin=100, genEnd=200,
            freqBegin=0.1, freqEnd=[0.14, 0.15])
        self.assertAlmostEqual(traj.freq(100, 0)[0], 0.1)
        self.assertTrue(traj.freq(200, 0)[0] >= 0.14 and traj.freq(200, 0)[0] <= 0.15)
        # no mutation
        self.assertEqual(len(traj.mutators()), 0)
        # the function form
        traj = ForwardTrajectory(N=1000, fitness=[1, 1.01, 1.03], genBegin=100,
            genEnd=200, freqBegin=0.1, freqEnd=[0.14, 0.15]) 
        #
        # FIXME: multi-subpopulation case
        # 

    def testTrajectorySimulatorSelMultiLoci(self):
        'Testing trajectory simulation of multiple loci with constant population size and natural selection'
        trajSimulator = trajectorySimulator(N=1000, fitness=[1, 1.01, 1.03], nLoci=2)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1)
        self.assertNotEqual(traj, None)
        # trajectory ends at 0.
        for gen in range(3000):
            self.assertEqual(len(traj.freq(gen, 0)), 2)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0., 0.])
        # most likely the starting generation will be different (there is a zero...)
        self.assertEqual(traj.freq(traj._beginGen() + 1, 0).count(0.), 1)
        # invalid inputs
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=3000, freqEnd=[0.1, 0.2, 0.3])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuBackward,
            genEnd=0, freqEnd=0.1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=200)
        self.assertTrue(traj._beginGen() <= 2800)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1, minMutAge=400,
            maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500 and traj._beginGen() <= 2600)
        # one mutation
        self.assertEqual(len(traj.mutators()), 2)
        # the function form
        traj = BackwardTrajectory(N=1000, nLoci=2, genEnd=3000, freqEnd=0.1) 
        #
        # simuForward
        # invalid input
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=0.09)
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.09])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[0.08, 0.09])
        self.assertRaises(exceptions.ValueError, trajSimulator.simuForward,
            genBegin=100, genEnd=200, freqBegin=0.1, freqEnd=[[0.18, 0.09]*2])
        # test trajectory
        traj = trajSimulator.simuForward(genBegin=100, genEnd=200,
            freqBegin=[0.1, 0.13], freqEnd=[[0.09, 0.11], [0.11, 0.12]])
        self.assertAlmostEqual(traj.freq(100, 0)[0], 0.1)
        self.assertAlmostEqual(traj.freq(100, 0)[1], 0.13)
        self.assertTrue(traj.freq(200, 0)[0] >= 0.09 and traj.freq(200, 0)[0] <= 0.11)
        self.assertTrue(traj.freq(200, 0)[1] >= 0.11 and traj.freq(200, 0)[0] <= 0.12)
        # no mutation
        self.assertEqual(len(traj.mutators()), 0)
        # the function form
        traj = ForwardTrajectory(N=1000, nLoci=2, genBegin=100, genEnd=200, freqBegin=0.1,
            freqEnd=[[0.09, 0.11]]*2) 
        #
        # Fix: multi subpopulation.
        #



    def testSimuCase2(self):
        'Testing trajectory simulation with interaction'
        # 2: dependent case when freq = 0 in backward simulation
        # and [0, 1] in forward simulation. nLoci = 2, constant nSubPop = 3:
        def Nt(gen):
            return [2000] * 3
        def fitness(gen, sp):
            return [1] * 9
        trajSimulator = trajectorySimulator(N=Nt, fitness=fitness, nLoci = 2)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0] * 2)
        traj1 = trajSimulator.simuForward(genBegin=0, freqBegin = [0, 1], freqEnd = [[0, 1]]*2, genEnd = 200)
        self.assertEqual(traj1.freq(random.randint(0, 200), 1), [0, 1])
        
        
    def testSimuCase4(self):
        'Testing trajectory simulation without selection'
        # 4: In backward simulation, given a normal set of parameters, if current freq of
        # allele a is very small, the first element of traj[] should be 1/(ploidy * N) and
        # the second element of traj[] should be larger than 1/(ploidy * N):
        # nLoci = 1 and no subpopulation
        fitness = [1,1,1]
        #### ??? N = 101
        trajSimulator = trajectorySimulator(N=1000, fitness=fitness, nLoci = 1)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.01)
        #####
        #print len(traj.freq)
        #print [traj.freq[min(traj.freq.keys()) + i] for i in range(3)]
        #####
        #self.assertAlmostEqual(traj.freq[min(traj.freq.keys())+1][0], 1. / (2 * 1000))
        #if traj.freq[min(traj.freq.keys())+2][0] < 1. / (2 * 1000):
        #   raise ValueError('fail in test number 4 of testSimuBackward.')
        
    def testSimuCase5(self):
        'Testing backward trajectory simulator'
        # 5: test given variable number of subpopulations, if allele frequencies would be
        # recorded in the correct form. nLoci = 3.
        # In a backward simulation, 3 subPops merge when gen = genEnd - 5 in backward sense.
        # In another Forward simulation, a single population splits into 3 subPops when 
        # gen = genEnd - 5 in forward sense. 
        def Nt(gen):
            if gen > 295:
                return [1000] * 3
            else:   # merge from all subpops into one population
                return 3000
        trajSimulator = trajectorySimulator(N=Nt, fitness=[1,1,1], nLoci = 3)
        traj = trajSimulator.simuBackward(genEnd=300, freqEnd=0.01)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0] * 3)
        self.assertEqual(len(traj.freq(296, 0)), 3)
        self.assertEqual(len(traj.freq(295, 0)), 3)
        traj1 = trajSimulator.simuForward(genBegin=0, genEnd=300, freqBegin=0.5, freqEnd=[[0,1]]*3)
        self.assertEqual(len(traj1.freq(295, 0)), 3)
        self.assertEqual(len(traj1.freq(296, 0)), 3)
        
        # nLoci = 2.
        # In a backward simulation, number of populations splits from one to 4 when
        # gen = genEnd - 5 in backward sense.
        # In another forward simulation, 4 subpopulations merge to one when
        # gen = genEnd - 5 in forward sense.
        def Nt(gen):
            if gen > 295:
                return 4000
            else:   # merge from all subpops into one population
                return [1000] * 4
        trajSimulator = trajectorySimulator(N=Nt, fitness=[1,1,1], nLoci = 3)
        traj = trajSimulator.simuBackward(genEnd=300, freqEnd=0.01)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0] * 3)
        self.assertEqual(len(traj.freq(296, 0)), 3)
        self.assertEqual(len(traj.freq(295, 0)), 3)
        traj1 = trajSimulator.simuForward(genBegin=0, genEnd=300, freqBegin=0.5, freqEnd=[[0,1]]*3)
        self.assertEqual(len(traj1.freq(295, 0)), 3)
        self.assertEqual(len(traj1.freq(296, 0)), 3)

    def testSimuCase6(self):
        'Testing backward trajectory simulator'
        # 6: given a normal set of parameters, plot a backward trajectory and
        # another forward trajectory for the case of single locus without
        # subpopulations.
        trajSimulator = trajectorySimulator(N=3000, fitness=[1,1,1], nLoci = 1)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=0.1)
        #traj.plot()
        traj1 = trajSimulator.simuForward(genBegin=0, genEnd=200, freqBegin=0.5, freqEnd=[[0,1]])
        #traj1.plot()
        
    def testSimuCase7(self):
        'Testing backward trajectory simulator'
        # 7: given a normal set of parameters, considering changable subpopulation
        # sizes with mulitiple loci, plot a backward trajectory and another
        # forward trajectory. nSubPops = 3, nLoci = 2.
        def Nt(gen):
            if gen > 2500:
                return [1000] * 3
            else:   # merge from all subpops into one population in backward sense
                return 3000
        trajSimulator = trajectorySimulator(N=Nt, fitness=[1,1,1], nLoci = 2)
        traj = trajSimulator.simuBackward(genEnd=3000, freqEnd=[0.05, 0.1])
        #traj.plot()
        def Nt1(gen):
            if gen < 100:
                return 3000
            else:   # split from one population into 3 subpops in forward sense
                return [1000] * 3
        trajSimulator1 = trajectorySimulator(N=Nt1, fitness=[1,1,1], nLoci = 5)
        traj1 = trajSimulator1.simuForward(genBegin=0, genEnd = 200,
            freqBegin = [0.5, 0.6, 0.7, 0.8, 0.9],
            freqEnd = [[0,1]]*5)
        #traj1.plot()
            

if __name__ == '__main__':
    unittest.main()
