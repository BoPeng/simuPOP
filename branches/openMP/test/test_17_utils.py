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
from simuPOP.utils import *
from simuPOP.gsl import *

import unittest, os, sys
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
        simu = Simulator( Population(size=10, ploidy=2, loci=[2, 3]),
            RandomMating(), reps=5)
        print("\n\nUSER INTERACTION: Please press q\n\n")
        self.assertRaises( SystemError, simu.evolve,
            postOps=[ Pause(at=[10]),
                        # should quite, can not reach generation 12
                        TerminateIf("True", at=[12] ) ] )

    def interactiveTestExitToShell(self):
        'Testing exit to a shell'
        simu = Simulator( Population(size=10, ploidy=2, loci=[2, 3]),
            RandomMating(), reps=5)
        print("\n\nUSER INTERACTION: Please press s and then Ctrl-D")
        print("Please check the existence of variable pop\n\n")
        simu.evolve(
            postOps=[ Pause(at=[10]) ], end=12)
        print("\n\nUSER INTERACTION: Please press s and then Ctrl-D")
        print("Please check the existence of variable tmpPop\n\n")
        simu.evolve(
            postOps=[ Pause(at=[20], popName='tmpPop') ], end=25)


    def testsetRNG(self):
        'Testing all RNG types'
        for rg in moduleInfo()['availableRNGs']:
            getRNG().set(rg)

    def testDefaultRNG(self):
        'Testing default RNG'
        rg = getRNG()
        self.assertEqual(rg.name(), 'mt19937')

    def testBinomial(self):
        'Testing binomial distribution'
        rg = getRNG()
        for n in range(1,10):
            rg.randBinomial(10, .7)

    def testUniform(self):
        'Testing uniform distribution generator'
        rg = getRNG()
        for n in range(1,10):
            rg.randUniform()

    def testRandBit(self):
        'Testing random bit function'
        rg = getRNG()
        sum = 0
        for n in range(10000):
            sum += rg.randBit()
        self.assertTrue(sum < 5100)
        self.assertTrue(sum > 4900)


    def testBernullitrials(self):
        'Testing bernullitrials'
        import math
        rg = getRNG()
        p = [0.00001, 0.001, 0.5, 0.99]
        N = 1000000
        bt = Bernullitrials(rg, p, N)
        bt.doTrial()
        for i in range(len(p)):
            prop = bt.trialSuccRate(i)
            # binomial, mean p, variance = p(1-p)/n
            std = math.sqrt(p[i]*(1.-p[i])/N)
            self.assertTrue(prop > p[i] - 3*std and prop < p[i] + 3*std)
        # another test, for each trail
        for pi in p:
            bt = Bernullitrials(rg, [pi]*N, 10)
            bt.doTrial()
            for i in range(10):
                bt.trial();
                prop = bt.probSuccRate()
                # binomial, mean p, variance = p(1-p)/n
                std = math.sqrt(pi*(1.-pi)/N)
                self.assertTrue(prop > pi - 3*std and prop < pi + 3*std, "We are testing if the proportion falls within three standard deviation of SD. This test may fail from time to time.")

        # test find_first and find_next
        bt = Bernullitrials(rg, p, N)
        bt.doTrial()
        for i in range(len(p)):
            pos = bt.trialFirstSucc(i)
            self.assertTrue(pos < N)
            if pos != bt.npos:
                self.assertTrue(bt.trialSucc(i, pos))
                for j in range(pos):
                    self.assertFalse(bt.trialSucc(i, j))
                while True:
                    last_pos = pos
                    pos = bt.trialNextSucc(i, pos)
                    if pos == bt.npos:
                        break;
                    self.assertTrue(pos < N)
                    self.assertTrue(bt.trialSucc(i, pos))
                    for j in range(last_pos+1, pos):
                        self.assertFalse(bt.trialSucc(i, j))
                for j in range(pos+1, N):
                    self.assertFalse(bt.trialSucc(i, j))

    def testSeed(self):
        'Testing RNG::seed() and RNG::setSeed()'
        import random
        seed = []
        name = moduleInfo()['availableRNGs'][0]
        for i in range(100):
            getRNG().set(name=name)
            sd = getRNG().seed()
            self.assertFalse(sd in seed)
            seed.append(sd)
        # test set seed
        sd = random.randint(100, 10000)
        getRNG().set(name, sd)
        self.assertEqual(getRNG().seed(), sd)
        # test if sequences are the same once the seed is set
        sd = random.randint(100, 10000)
        getRNG().set(name, sd)
        seq = [getRNG().randInt(10000) for x in range(100)]
        getRNG().set(name, sd)
        seq1 = [getRNG().randInt(10000) for x in range(100)]
        self.assertEqual(seq, seq1)
        # randBit need to be treated separately because it uses
        # global variable of RNG().
        sd = random.randint(100, 10000)
        getRNG().set(name, sd)
        seq = [getRNG().randBit() for x in range(100)]
        getRNG().set(name, sd)
        seq1 = [getRNG().randBit() for x in range(100)]
        self.assertEqual(seq, seq1)

    def testWeightedSampler(self):
        'Testing weighted sampler'
        sampler = WeightedSampler([1, 2, 3, 4])
        num = []
        for i in range(100000):
            num.append(sampler.draw())
        for i in range(4):
            self.assertAlmostEqual(num.count(i) / 100000., 0.1 * (i+1), places=1)
        num = sampler.drawSamples(100000)
        for i in range(4):
            self.assertAlmostEqual(num.count(i) / 100000., 0.1 * (i+1), places=1)
        # the proportion version
        sampler = WeightedSampler([0.1, 0.2, 0.3, 0.4], 100000)
        num = sampler.drawSamples(100000)
        for i in range(4):
            # the count must be exact
            self.assertEqual(num.count(i), 10000 * (i+1))

    def TestLargePedigree(self):
        'Testing getting large Pedigree, for simuUtils.ascertainPedigree'
        import simuUtil
        pop = Population(100, ancestralDepth=2, infoFields=['father_idx', 'mother_idx'])
        simu = Simulator(pop, RandomMating(numOffspring=0.3, mode=MATE_GeometricDistribution))
        simu.evolve(duringOps=[ParentsTagger()], end=5)
        pop = simu.population(0)
        #
        def comb(geno):
            return sum(geno)+1
        SampleLargePedigree(pop, numPedigree=10, minPedSize=5, minAffected=0,
            maxOffspring=5, output='ped', loci=[], combine=comb)
        simuUtil.VC_merlin('ped')


    def testMemoryLeakloadPopulation(self):
        'Testing if loadPopulation leaks memory'
        # run this and see if memory usage goes up continuously
        pop = Population(100, loci=[400])
        initGenotype(pop, freq=[0.2, 0.8])
        stat(pop, alleleFreq=list(range(pop.totNumLoci())))
        pop.save('test.bin')
        for i in range(4):
            pop = loadPopulation('test.bin')
            #pop = Population(100, loci=[1000]*10)
            stat(pop, alleleFreq=list(range(pop.totNumLoci())))
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

    def testTrajectorySimulatorInputs(self):
        'Testing invalid inputs for Trajectory simulation functions'
        # both backward and forward
        self.assertRaises(ValueError, simulateForwardTrajectory,
            N=1000, fitness=[1, 1.01], beginGen=100, endGen=1000,
            beginFreq=0.1, endFreq=[0.1, 0.2])
        self.assertRaises(ValueError, simulateForwardTrajectory,
            N=1000, beginGen=2000, endGen=1000,
            beginFreq=0.1, endFreq=[0.2, 0.3])
        self.assertRaises(ValueError, simulateForwardTrajectory,
            N=1000, fitness=[1]*5, beginGen=0, endGen=1000,
            beginFreq=0.1, endFreq=[0.2, 0.3])
        #
        # forward Trajectory simulation
        self.assertRaises(ValueError, simulateForwardTrajectory,
            N=1000, beginGen=100, endGen=200, beginFreq=0.1, endFreq=0.09)
        self.assertRaises(ValueError, simulateForwardTrajectory,
            N=1000, beginGen=100, endGen=200, beginFreq=0.1, endFreq=[0.09])
        self.assertRaises(ValueError, simulateForwardTrajectory,
            N=1000, nLoci=2, beginGen=100, endGen=200, beginFreq=0.1,
            endFreq=[0.08, 0.09])
        self.assertRaises(ValueError, simulateForwardTrajectory,
            N=1000, nLoci=2, beginGen=100, endGen=200, beginFreq=0.1,
            endFreq=[[0.18, 0.09]*2])
        #
        # backward Trajectory simulation
        self.assertRaises(ValueError, simulateBackwardTrajectory,
            N=1000, endGen=3000, endFreq=[0.1, 0.2])
        self.assertRaises(ValueError, simulateBackwardTrajectory,
            N=1000, endGen=0, endFreq=0.1)
        self.assertRaises(ValueError, simulateBackwardTrajectory,
            N=1000, endGen=3000, endFreq=[0.1, 0.2, 0.3])
        self.assertRaises(ValueError, simulateBackwardTrajectory,
            N=1000, endGen=0, endFreq=0.1)

    def testBackTrajectorySimple(self):
        'Testing backward Trajectory simulation with constant population size and no selection'
        trajSimulator = TrajectorySimulator(N=1000)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1)
        # Trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0.])
        # Trajectory starts at 0.1
        self.assertAlmostEqual(traj.freq(3000, 0)[0], 0.1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=400,
            maxMutAge=1000)
        self.assertTrue(traj._beginGen() >= 2000 and traj._beginGen() <= 2600)
        # one mutation
        self.assertEqual(len(traj.mutators(loci=[0])), 1)

    def testBackTrajectorySubpop(self):
        'Testing backward Trajectory simulation with subpopulations and no selection'
        # simple with subpopulation
        trajSimulator = TrajectorySimulator(N=[1000, 2000])
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1)
        self.assertNotEqual(traj, None)
        # Trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0.])
        self.assertEqual(traj.freq(traj._beginGen(), 1), [0.])
        self.assertAlmostEqual(traj.freq(3000, 0)[0], 0.1)
        self.assertAlmostEqual(traj.freq(3000, 1)[0], 0.1)
        # different destination frequency
        trajSimulator = TrajectorySimulator(N=[1000, 2000])
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=[0.1, 0.2])
        # Trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0.])
        self.assertEqual(traj.freq(traj._beginGen(), 1), [0.])
        self.assertAlmostEqual(traj.freq(3000, 0)[0], 0.1)
        self.assertAlmostEqual(traj.freq(3000, 1)[0], 0.2)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=[0.1, 0.2], minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=[0.1, 0.2], maxMutAge=1000)
        self.assertTrue(traj._beginGen() >= 2000)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=[0.1, 0.2], minMutAge=400,
            maxMutAge=1000)
        self.assertTrue(traj._beginGen() >= 2000 and traj._beginGen() <= 2600)
        # one mutation in each subpopulation
        self.assertEqual(len(traj.mutators(loci=[0])), 2) 
    
    def testsimulateForwardTrajectorySimple(self):
        'Testing forward Trajectory simulation with constant population size and no selection'
        trajSimulator = TrajectorySimulator(N=1000)
        traj = trajSimulator.simuForward(beginGen=100, endGen=200,
            beginFreq=0.1, endFreq=[0.09, 0.11])
        self.assertTrue(traj.freq(200, 0)[0] >= 0.09 and traj.freq(200, 0)[0] <= 0.11)
        # no mutation
        self.assertEqual(len(traj.mutators(loci=[0])), 0)
        #

    def testsimulateForwardTrajectorySubPop(self):
        'Testing forward Trajectory simulation with subpopulations and no selection'
        # simple with subpopulation
        trajSimulator = TrajectorySimulator(N=[1000, 2000])
        # simuForward, test Trajectory
        traj = trajSimulator.simuForward(beginGen=100, endGen=200,
            beginFreq=0.1, endFreq=[0.09, 0.11])
        fq = (traj.freq(200, 0)[0] + traj.freq(200, 1)[0]*2)/3.
        self.assertTrue(fq >= 0.09 and fq <= 0.11)
        # no mutation
        self.assertEqual(len(traj.mutators(loci=[0])), 0)

    def testForwardSelBySubPop(self):
        'Testing forward Trajectory simulation with subpop specific selection'
        def fitness(gen, sp):
            if sp == 0:
                return [1, 1.001, 1.002]
            else:
                return [1, 0.999, 0.998]
        #
        traj = simulateForwardTrajectory(N=[2000, 4000], fitness=fitness, beginGen=0,
            endGen=150, beginFreq=[0.2, 0.2], endFreq=[[0.25, 0.35], [0.1, 0.15]])
        self.assertTrue(traj.freq(150, 0)[0] >= 0.25 and traj.freq(150, 0)[0] <= 0.35)
        self.assertTrue(traj.freq(150, 1)[0] >= 0.10 and traj.freq(150, 1)[0] <= 0.15)


    def testBackwardSimpleMultiLoci(self):
        'Testing backward Trajectory simulation of multiple loci with constant population size and no selection'
        trajSimulator = TrajectorySimulator(N=1000, nLoci=2)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1)
        self.assertNotEqual(traj, None)
        # Trajectory ends at 0.
        for gen in range(3000):
            self.assertEqual(len(traj.freq(gen, 0)), 2)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0., 0.])
        # most likely the starting generation will be different (there is a zero...)
        self.assertEqual(traj.freq(traj._beginGen() + 1, 0).count(0.), 1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=400,
            maxMutAge=1000)
        self.assertTrue(traj._beginGen() >= 2000 and traj._beginGen() <= 2600)
        # one mutation
        self.assertEqual(len(traj.mutators(loci=[0, 1])), 2)

    def testForwardSimpleMultiLoci(self):
        'Testing forward Trajectory simulation of multiple loci with constant population size and no selection'
        # test Trajectory
        trajSimulator = TrajectorySimulator(N=1000, nLoci=2)
        traj = trajSimulator.simuForward(beginGen=100, endGen=200,
            beginFreq=[0.1, 0.13], endFreq=[[0.09, 0.11], [0.11, 0.12]])
        self.assertTrue(traj.freq(200, 0)[0] >= 0.09 and traj.freq(200, 0)[0] <= 0.11)
        self.assertTrue(traj.freq(200, 0)[1] >= 0.11 and traj.freq(200, 0)[0] <= 0.12)
        # no mutation
        self.assertEqual(len(traj.mutators(loci=[0, 1])), 0)
        # the function form
        traj = simulateForwardTrajectory(N=1000, nLoci=2, beginGen=100, endGen=200, beginFreq=0.1,
            endFreq=[[0.09, 0.11]]*2) 

    def testBackwardSimpleMultiLociSubPop(self):
        'Testing backward Trajectory simulation of multiple loci with multiple populations and no selection'
        trajSimulator = TrajectorySimulator(N=[1000, 3000], nLoci=2)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1)
        self.assertNotEqual(traj, None)
        # Trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0., 0.])
        # most likely the starting generation will be different (there is three zeros...)
        self.assertTrue(traj.freq(traj._beginGen() + 1, 0).count(0.) > 1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, maxMutAge=1500)
        self.assertTrue(traj._beginGen() >= 1500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=400,
            maxMutAge=1500)
        self.assertTrue(traj._beginGen() >= 1500 and traj._beginGen() <= 2600 )
        # one mutation
        self.assertEqual(len(traj.mutators(loci=[0,1])), 4)

    def testForwardSimpleMultiLociSubPop(self):
        'Testing forward Trajectory simulation of multiple loci with multiple populations and no selection'
        # test Trajectory
        trajSimulator = TrajectorySimulator(N=[1000, 3000], nLoci=2)
        traj = trajSimulator.simuForward(beginGen=100, endGen=200,
            beginFreq=[0.1, 0.13], endFreq=[[0.08, 0.13], [0.10, 0.15]])
        f1 = (traj.freq(200, 0)[0] + traj.freq(200, 1)[0]*3)/4.
        f2 = (traj.freq(200, 0)[1] + traj.freq(200, 1)[1]*3)/4.
        self.assertTrue(f1 >= 0.08 and f1 <= 0.13)
        self.assertTrue(f2 >= 0.10 and f2 <= 0.15)
        # no mutation
        self.assertEqual(len(traj.mutators(loci=[0,1])), 0)
        # the function form
        traj = simulateForwardTrajectory(N=[1000, 3000], nLoci=2, beginGen=100, endGen=200, beginFreq=0.1,
            endFreq=[[0.09, 0.11]]*2) 
        #


    def testTrajectorySimulatorSel(self):
        'Testing Trajectory simulation with constant population size and natural selection'
        # invalid input
        trajSimulator = TrajectorySimulator(N=1000, fitness=[1, 1.01, 1.03])
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1)
        self.assertNotEqual(traj, None)
        # Trajectory ends at 0.
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0.])
        self.assertAlmostEqual(traj.freq(3000, 0)[0], 0.1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=500)
        self.assertTrue(traj._beginGen() <= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=400,
            maxMutAge=1000)
        self.assertTrue(traj._beginGen() >= 2000 and traj._beginGen() <= 2600)
        # one mutation
        self.assertEqual(len(traj.mutators(loci=[0])), 1)
        # the function form
        traj = simulateBackwardTrajectory(N=1000, fitness=[1, 1.01, 1.03], endGen=3000, endFreq=0.1) 
        #
        # simuForward
        # test Trajectory
        traj = trajSimulator.simuForward(beginGen=100, endGen=200,
            beginFreq=0.1, endFreq=[0.14, 0.15])
        self.assertTrue(traj.freq(200, 0)[0] >= 0.14 and traj.freq(200, 0)[0] <= 0.15)
        # no mutation
        self.assertEqual(len(traj.mutators(loci=[0])), 0)
        # the function form
        traj = simulateForwardTrajectory(N=1000, fitness=[1, 1.01, 1.03], beginGen=100,
            endGen=200, beginFreq=0.1, endFreq=[0.14, 0.15]) 
        #
        # FIXME: multi-subpopulation case
        # 

    def testTrajectorySimulatorSelMultiLoci(self):
        'Testing Trajectory simulation of multiple loci with constant population size and natural selection'
        trajSimulator = TrajectorySimulator(N=1000, fitness=[1, 1.01, 1.03], nLoci=2)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1)
        self.assertNotEqual(traj, None)
        # Trajectory ends at 0.
        for gen in range(3000):
            self.assertEqual(len(traj.freq(gen, 0)), 2)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0., 0.])
        # most likely the starting generation will be different (there is a zero...)
        self.assertEqual(traj.freq(traj._beginGen() + 1, 0).count(0.), 1)
        # checking minMutAge and maxMutAge
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=200)
        self.assertTrue(traj._beginGen() <= 2800)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, maxMutAge=500)
        self.assertTrue(traj._beginGen() >= 2500)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1, minMutAge=400,
            maxMutAge=1000)
        self.assertTrue(traj._beginGen() >= 2000 and traj._beginGen() <= 2600)
        # one mutation
        self.assertEqual(len(traj.mutators(loci=[0,1])), 2)
        # the function form
        traj = simulateBackwardTrajectory(N=1000, nLoci=2, endGen=3000, endFreq=0.1) 
        #
        # simuForward
        traj = trajSimulator.simuForward(beginGen=100, endGen=200,
            beginFreq=[0.1, 0.13], endFreq=[[0.12, 0.14], [0.15, 0.19]])
        self.assertTrue(traj.freq(200, 0)[0] >= 0.12 and traj.freq(200, 0)[0] <= 0.14)
        self.assertTrue(traj.freq(200, 0)[1] >= 0.15 and traj.freq(200, 0)[0] <= 0.19)
        # no mutation
        self.assertEqual(len(traj.mutators(loci=[0,1])), 0)
        # the function form
        traj = simulateForwardTrajectory(N=1000, nLoci=2, beginGen=100, endGen=200, beginFreq=0.1,
            endFreq=[[0.09, 0.11]]*2) 
        #
        # Fix: multi subpopulation.
        #



    def testSimuCase2(self):
        'Testing Trajectory simulation with interaction'
        # 2: dependent case when freq = 0 in backward simulation
        # and [0, 1] in forward simulation. nLoci = 2, constant nSubPop = 3:
        def Nt(gen):
            return [2000] * 3
        def fitness(gen, sp):
            return [1] * 9
        trajSimulator = TrajectorySimulator(N=Nt, fitness=fitness, nLoci = 2)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0] * 2)
        traj1 = trajSimulator.simuForward(beginGen=0, beginFreq = [0, 1], endFreq = [[0, 1]]*2, endGen = 200)
        self.assertEqual(traj1.freq(random.randint(0, 200), 1), [0, 1])
        
        
    def testSimuCase4(self):
        'Testing Trajectory simulation without selection'
        # 4: In backward simulation, given a normal set of parameters, if current freq of
        # allele a is very small, the first element of traj[] should be 1/(ploidy * N) and
        # the second element of traj[] should be larger than 1/(ploidy * N):
        # nLoci = 1 and no subpopulation
        fitness = [1,1,1]
        #### ??? N = 101
        trajSimulator = TrajectorySimulator(N=1000, fitness=fitness, nLoci = 1)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.01)
        #####
        #print len(traj.freq)
        #print [traj.freq[min(traj.freq.keys()) + i] for i in range(3)]
        #####
        #self.assertAlmostEqual(traj.freq[min(traj.freq.keys())+1][0], 1. / (2 * 1000))
        #if traj.freq[min(traj.freq.keys())+2][0] < 1. / (2 * 1000):
        #   raise ValueError('fail in test number 4 of testSimuBackward.')
        
    def testSimuCase5(self):
        'Testing backward Trajectory Simulator'
        # 5: test given variable number of subpopulations, if allele frequencies would be
        # recorded in the correct form. nLoci = 3.
        # In a backward simulation, 3 subPops merge when gen = endGen - 5 in backward sense.
        # In another Forward simulation, a single population splits into 3 subPops when 
        # gen = endGen - 5 in forward sense. 
        def Nt(gen):
            if gen > 295:
                return [1000] * 3
            else:   # merge from all subpops into one population
                return 3000
        trajSimulator = TrajectorySimulator(N=Nt, fitness=[1,1,1], nLoci = 3)
        traj = trajSimulator.simuBackward(endGen=300, endFreq=0.01)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0] * 3)
        self.assertEqual(len(traj.freq(296, 0)), 3)
        self.assertEqual(len(traj.freq(295, 0)), 3)
        traj1 = trajSimulator.simuForward(beginGen=0, endGen=300, beginFreq=0.5, endFreq=[[0,1]]*3)
        self.assertEqual(len(traj1.freq(295, 0)), 3)
        self.assertEqual(len(traj1.freq(296, 0)), 3)
        # nLoci = 2.
        # In a backward simulation, number of populations splits from one to 4 when
        # gen = endGen - 5 in backward sense.
        # In another forward simulation, 4 subpopulations merge to one when
        # gen = endGen - 5 in forward sense.
        def Nt(gen):
            if gen > 295:
                return 4000
            else:   # merge from all subpops into one population
                return [1000] * 4
        trajSimulator = TrajectorySimulator(N=Nt, fitness=[1,1,1], nLoci = 3)
        traj = trajSimulator.simuBackward(endGen=300, endFreq=0.01)
        self.assertEqual(traj.freq(traj._beginGen(), 0), [0] * 3)
        self.assertEqual(len(traj.freq(296, 0)), 3)
        self.assertEqual(len(traj.freq(295, 0)), 3)
        traj1 = trajSimulator.simuForward(beginGen=0, endGen=300, beginFreq=0.5, endFreq=[[0,1]]*3)
        self.assertEqual(len(traj1.freq(295, 0)), 3)
        self.assertEqual(len(traj1.freq(296, 0)), 3)

    def testSimuCase6(self):
        'Testing backward Trajectory Simulator'
        # 6: given a normal set of parameters, plot a backward Trajectory and
        # another forward Trajectory for the case of single locus without
        # subpopulations.
        trajSimulator = TrajectorySimulator(N=3000, fitness=[1,1,1], nLoci = 1)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=0.1)
        #traj.plot()
        traj1 = trajSimulator.simuForward(beginGen=0, endGen=200, beginFreq=0.5, endFreq=[[0,1]])
        #traj1.plot()
        
    def testSimuCase7(self):
        'Testing backward Trajectory Simulator'
        # 7: given a normal set of parameters, considering changable subpopulation
        # sizes with mulitiple loci, plot a backward Trajectory and another
        # forward Trajectory. nSubPops = 3, nLoci = 2.
        def Nt(gen):
            if gen > 2500:
                return [1000] * 3
            else:   # merge from all subpops into one population in backward sense
                return 3000
        trajSimulator = TrajectorySimulator(N=Nt, fitness=[1,1,1], nLoci = 2)
        traj = trajSimulator.simuBackward(endGen=3000, endFreq=[0.05, 0.1])
        #traj.plot()
        def Nt1(gen):
            if gen < 100:
                return 3000
            else:   # split from one population into 3 subpops in forward sense
                return [1000] * 3
        trajSimulator1 = TrajectorySimulator(N=Nt1, fitness=[1,1,1], nLoci = 5)
        traj1 = trajSimulator1.simuForward(beginGen=0, endGen = 200,
            beginFreq = [0.5, 0.6, 0.7, 0.8, 0.9],
            endFreq = [[0,1]]*5)
        #traj1.plot()
    
    def testGSL(self):
        'Testing GSL functions'
        gsl_cdf_gaussian_P(0.5, 1)
        gsl_cdf_gaussian_Q(0.5, 1)
        gsl_cdf_gaussian_Pinv(0.5, 1)
        gsl_cdf_gaussian_Qinv(0.5, 1)
        gsl_cdf_ugaussian_P(0.5)
        gsl_cdf_ugaussian_Q(0.5)
        gsl_cdf_ugaussian_Pinv(0.5)
        gsl_cdf_ugaussian_Qinv(0.5)
        #
        gsl_cdf_exponential_P(0.5, 0.1)
        gsl_cdf_exponential_Q(0.5, 0.1)
        gsl_cdf_exponential_Pinv(0.5, 0.1)
        gsl_cdf_exponential_Qinv(0.5, 0.1)
        #
        gsl_cdf_chisq_P(0.5, 0.1)
        gsl_cdf_chisq_Q(0.5, 0.1)
        gsl_cdf_chisq_Pinv(0.5, 0.1)
        gsl_cdf_chisq_Qinv(0.5, 0.1)
        #
        gsl_cdf_gamma_P(0.5, 0.1, 0.2)
        gsl_cdf_gamma_Q(0.5, 0.1, 0.2)
        gsl_cdf_gamma_Pinv(0.5, 0.1, 0.2)
        gsl_cdf_gamma_Qinv(0.5, 0.1, 0.2)
        #
        gsl_cdf_beta_P(0.5, 0.1, 0.2)
        gsl_cdf_beta_Q(0.5, 0.1, 0.2)
        gsl_cdf_beta_Pinv(0.5, 0.1, 0.2)
        gsl_cdf_beta_Qinv(0.5, 0.1, 0.2)
        #
        gsl_cdf_binomial_P(5, 0.1, 100)
        gsl_cdf_binomial_Q(5, 0.1, 100)
        gsl_ran_binomial_pdf(5, 0.1, 100)
        #
        gsl_cdf_poisson_P(2, 3.4)
        gsl_cdf_poisson_Q(2, 3.4)
        gsl_ran_poisson_pdf(2, 3.4)

if __name__ == '__main__':
    unittest.main()
