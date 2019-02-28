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

import random
import gzip
import unittest, os, sys
from simuOpt import setOptions
# this line also tests the use of parameter version in setOptions
setOptions(quiet=True, gui='batch', version='1.0.1')
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
from simuPOP.sampling import *
from simuPOP.utils import *
from simuPOP.gsl import *

# for memory leak testing.

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}


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

class TestUtils(unittest.TestCase):

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

    def testDefDictType(self):
        'Testing type defdict'
        d = defdict()
        self.assertEqual(d[1], 0)
        self.assertEqual(d[10], 0)
        self.assertEqual(d[20], 0)
        d[4] = 5
        self.assertEqual(d[4], 5)

    def testsetRNG(self):
        'Testing all RNG types'
        old_rng = getRNG().name()
        for rg in moduleInfo()['availableRNGs']:
            getRNG().set(rg)
        setRNG(name=old_rng)

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
                # 
                # We need a more formatl test for distribution
                #
                #self.assertTrue(prop > pi - 3*std and prop < pi + 3*std, "We are testing if the proportion falls within three standard deviation of SD. This test may fail from time to time.")

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

    def testBernullitrials_T(self):
        'Testing bernullitrials_T'
        import math
        rg = getRNG()
        p = [0.001, 0.001, 0.5, 0.99]
        N = 1000000
        bt = Bernullitrials_T(rg, p, N)
        bt.doTrial()
        for i in range(len(p)):
            prop = bt.trialSuccRate(i)
            # binomial, mean p, variance = p(1-p)/n
            std = math.sqrt(p[i]*(1.-p[i])/N)
            self.assertTrue(prop > p[i] - 3*std and prop < p[i] + 3*std)
        # another test, for each trail
        for pi in p:
            bt = Bernullitrials_T(rg, [pi]*N, 10)
            bt.doTrial()
            for i in range(10):
                bt.trial();
                prop = bt.probSuccRate()
                # binomial, mean p, variance = p(1-p)/n
                std = math.sqrt(pi*(1.-pi)/N)
                # 
                # FIXME: we need a more rigirous test for distribution
                #
                #self.assertTrue(prop > pi - 3*std and prop < pi + 3*std, "Observed proprotin of bernullitrials {} not within {} of designed proportion {}".format(prop, pi, 2*std))

        # test find_first and find_next
        nP = 10000
        bt = Bernullitrials_T(rg, [0.001]*10000, 10)
        for i in range(100):
            bt.doTrial()
            #all_pos = []
            pos = bt.probFirstSucc()
            if pos == bt.npos:
                continue
            #all_pos.append(pos)
            self.assertTrue(pos < nP)
            if pos != bt.npos:
                self.assertTrue(bt.trialSucc(pos))
                for j in range(pos):
                    self.assertFalse(bt.trialSucc(j))
                while True:
                    last_pos = pos
                    pos = bt.probNextSucc(pos)
                    if pos == bt.npos:
                        break;
                    #all_pos.append(pos)
                    self.assertTrue(pos < nP)
                    self.assertTrue(bt.trialSucc(pos))
                    for j in range(last_pos+1, pos):
                        self.assertFalse(bt.trialSucc(j))
                for j in range(pos+1, N):
                    self.assertFalse(bt.trialSucc(j))
            #print(all_pos)



    def testSeed(self):
        'Testing RNG::seed() and RNG::setSeed()'
        import random
        seed = []
        name = moduleInfo()['availableRNGs'][0]
        old_rng = getRNG().name()
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
        setRNG(name=old_rng)

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

    def testWeightedSamplerWithZero(self):
        'Testing weighted sampler with Zero'
        sampler = WeightedSampler([0, 1, 2, 0, 0, 3, 4, 0])
        num = []
        for i in range(100000):
            num.append(sampler.draw())
        self.assertAlmostEqual(num.count(1) / 100000., 0.1 * 1, places=1)
        self.assertAlmostEqual(num.count(2) / 100000., 0.1 * 2, places=1)
        self.assertAlmostEqual(num.count(5) / 100000., 0.1 * 3, places=1)
        self.assertAlmostEqual(num.count(6) / 100000., 0.1 * 4, places=1)
        self.assertEqual(num.count(0), 0)
        self.assertEqual(num.count(3), 0)
        self.assertEqual(num.count(4), 0)
        self.assertEqual(num.count(7), 0)
        num = sampler.drawSamples(100000)
        self.assertAlmostEqual(num.count(1) / 100000., 0.1 * 1, places=1)
        self.assertAlmostEqual(num.count(2) / 100000., 0.1 * 2, places=1)
        self.assertAlmostEqual(num.count(5) / 100000., 0.1 * 3, places=1)
        self.assertAlmostEqual(num.count(6) / 100000., 0.1 * 4, places=1)
        self.assertEqual(num.count(0), 0)
        self.assertEqual(num.count(3), 0)
        self.assertEqual(num.count(4), 0)
        self.assertEqual(num.count(7), 0)
        # the proportion version
        sampler = WeightedSampler([0, 0.1, 0.2, 0, 0, 0.3, 0.4, 0], 100000)
        num = sampler.drawSamples(100000)
        self.assertEqual(num.count(1), 10000)
        self.assertEqual(num.count(2), 20000)
        self.assertEqual(num.count(5), 30000)
        self.assertEqual(num.count(6), 40000)
        self.assertEqual(num.count(0), 0)
        self.assertEqual(num.count(3), 0)
        self.assertEqual(num.count(4), 0)
        self.assertEqual(num.count(7), 0)
        #
        sampler = WeightedSampler([0, 1, 1, 0, 0, 1, 1, 0])
        num = []
        for i in range(100000):
            num.append(sampler.draw())
        self.assertAlmostEqual(num.count(1) / 100000., 0.25 , places=1)
        self.assertAlmostEqual(num.count(2) / 100000., 0.25 , places=1)
        self.assertAlmostEqual(num.count(5) / 100000., 0.25 , places=1)
        self.assertAlmostEqual(num.count(6) / 100000., 0.25 , places=1)
        self.assertEqual(num.count(0), 0)
        self.assertEqual(num.count(3), 0)
        self.assertEqual(num.count(4), 0)
        self.assertEqual(num.count(7), 0)
        num = sampler.drawSamples(100000)

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
        #self.assertTrue(traj._beginGen() >= 2000 and traj._beginGen() <= 2600)
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

    def lineOfFile(self, filename, line):
        with open(filename) as input:
            for i in range(line):
                res = input.readline()
        return res

    def testExportFStat(self):
        '''Testing export population in fstat format'''
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'], infoFields=['loc', 'pheno']) 
        initGenotype(pop, haplotypes=[0,1])
        initSex(pop, sex=[MALE, FEMALE])
        #
        # basic export
        export(pop, format='fstat', output='fstat.txt')
        self.assertEqual(self.lineOfFile('fstat.txt', 1), '2 6 2 1\n')
        self.assertEqual(self.lineOfFile('fstat.txt', 2), 'a\n')
        self.assertEqual(self.lineOfFile('fstat.txt', 8), '1 11 22 11 22 11 22\n')
        # 
        # test parameter lociNames
        self.assertRaises(ValueError, export, pop, format='fstat', output='fstat.txt',
            lociNames=['aa', 'bb'])
        export(pop, format='fstat', output='fstat.txt', lociNames=['aa', 'bb', 'c', 'd', 'e', 'f'])
        self.assertEqual(self.lineOfFile('fstat.txt', 2), 'aa\n')
        #
        # test parameter adjust
        export(pop, format='fstat', output='fstat.txt', adjust=0)
        self.assertEqual(self.lineOfFile('fstat.txt', 8), '1 00 11 00 11 00 11\n')
        #
        # test the case with no loci name and no specified name
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            infoFields=['loc', 'pheno']) 
        initGenotype(pop, haplotypes=[0,1])
        initSex(pop, sex=[MALE, FEMALE])
        export(pop, format='fstat', output='fstat.txt')
        self.assertEqual(self.lineOfFile('fstat.txt', 8), '1 11 22 11 22 11 22\n')

        # cleanup
        os.remove('fstat.txt')

    def testImportFStat(self):
        '''Testing import population from a file in FStat format'''
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'], infoFields=['loc', 'pheno']) 
        initGenotype(pop, freq=[0.2, 0.3, 0.4, 0.1])
        genotypes = list(pop.genotype())
        initSex(pop, sex=[MALE, FEMALE])
        export(pop, format='fstat', output='fstat.txt')
        #
        pop = importPopulation(format='fstat', filename='fstat.txt')
        self.assertEqual(pop.numSubPop(), 2)
        self.assertEqual(pop.numChrom(), 1)
        self.assertEqual(pop.totNumLoci(), 6)
        self.assertEqual(pop.ploidy(), 2)
        self.assertEqual(pop.lociNames(), ('a', 'b', 'c', 'd', 'e', 'f'))
        if moduleInfo()['alleleType'] == 'binary':
            self.assertEqual(list(pop.genotype()), [x+1 >= 1 for x in genotypes])
        else:
            self.assertEqual(list(pop.genotype()), [x+1 for x in genotypes])
        del pop
        # test for parameter adjust
        pop = importPopulation(format='fstat', filename='fstat.txt', adjust=-1)
        self.assertEqual(pop.numSubPop(), 2)
        self.assertEqual(pop.numChrom(), 1)
        self.assertEqual(pop.totNumLoci(), 6)
        self.assertEqual(pop.ploidy(), 2)
        self.assertEqual(pop.lociNames(), ('a', 'b', 'c', 'd', 'e', 'f'))
        self.assertEqual(list(pop.genotype()), [x for x in genotypes])
        # cleanup
        os.remove('fstat.txt')

    def testExportGenePop(self):
        '''Testing export population in genepop format'''
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'], infoFields=['loc', 'pheno']) 
        initGenotype(pop, haplotypes=[0,1])
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        #
        # test parameter title
        export(pop, format='genepop', output='genepop.txt', title='something')
        self.assertEqual(self.lineOfFile('genepop.txt', 1), 'something\n')
        # pending \n does not matter
        export(pop, format='genepop', output='genepop.txt', title='something\n')
        self.assertEqual(self.lineOfFile('genepop.txt', 1), 'something\n')
        # first line?
        export(pop, format='genepop', output='genepop.txt')
        self.assertEqual(self.lineOfFile('genepop.txt', 1)[:10], 'Outputted ')
        # loci names
        self.assertEqual(self.lineOfFile('genepop.txt', 2), 'a, b, c, d, e, f\n')
        # pop?
        self.assertEqual(self.lineOfFile('genepop.txt', 3), 'POP\n')
        # genotype?
        self.assertEqual(self.lineOfFile('genepop.txt', 4), 'SubPop0-1, 0101 0202 0101 0202 0101 0202\n')
        # large allele?
        if moduleInfo()['alleleType'] != 'binary':
            pop.genotype()[0] = 100
            export(pop, format='genepop', output='genepop.txt')
            self.assertEqual(self.lineOfFile('genepop.txt', 4), 'SubPop0-1, 101001 002002 001001 002002 001001 002002\n')
        # cleanup
        os.remove('genepop.txt')
        
    def testImportGenePop(self):
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'], infoFields=['loc', 'pheno']) 
        initGenotype(pop, freq=[0.2, 0.3, 0.4, 0.1])
        genotypes = list(pop.genotype())
        initSex(pop, sex=[MALE, FEMALE])
        export(pop, format='genepop', output='genepop.txt')
        #
        pop = importPopulation(format='genepop', filename='genepop.txt')
        self.assertEqual(pop.numSubPop(), 2)
        self.assertEqual(pop.numChrom(), 1)
        self.assertEqual(pop.totNumLoci(), 6)
        self.assertEqual(pop.ploidy(), 2)
        self.assertEqual(pop.lociNames(), ('a', 'b', 'c', 'd', 'e', 'f'))
        if moduleInfo()['alleleType'] == 'binary':
            self.assertEqual(list(pop.genotype()), [x+1 >= 1 for x in genotypes])
        else:
            self.assertEqual(list(pop.genotype()), [x+1 for x in genotypes])
        del pop
        # test for parameter adjust
        pop = importPopulation(format='genepop', filename='genepop.txt', adjust=-1)
        self.assertEqual(pop.numSubPop(), 2)
        self.assertEqual(pop.numChrom(), 1)
        self.assertEqual(pop.totNumLoci(), 6)
        self.assertEqual(pop.ploidy(), 2)
        self.assertEqual(pop.lociNames(), ('a', 'b', 'c', 'd', 'e', 'f'))
        self.assertEqual(list(pop.genotype()), [x for x in genotypes])
        # cleanup
        os.remove('genepop.txt')


    def testExportStructure(self):
        '''Testing export population in structure format'''
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'], infoFields=['loc', 'pheno']) 
        initGenotype(pop, haplotypes=[0,1])
        initInfo(pop, [20], infoFields='loc')
        initInfo(pop, [30], infoFields='pheno')
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        export(pop, format='structure', output='stru.txt')
        self.assertEqual(self.lineOfFile('stru.txt', 1), 'a\tb\tc\td\te\tf\n')
        # test for parameter markerNames
        self.assertRaises(ValueError, export, pop, format='structure', output='stru.txt', markerNames='cc')
        self.assertRaises(ValueError, export, pop, format='structure', output='stru.txt', markerNames=['cc'])
        export(pop, format='structure', output='stru.txt', markerNames=['f', 'e', 'd', 'c', 'b', 'a'])
        self.assertEqual(self.lineOfFile('stru.txt', 1), 'f\te\td\tc\tb\ta\n')
        # test for parameter recessive alleles
        self.assertRaises(ValueError, export, pop, format='structure', output='stru.txt',
            recessiveAlleles=2)
        export(pop, format='structure', output='stru.txt', recessiveAlleles=0)
        self.assertEqual(self.lineOfFile('stru.txt', 2), '0\n')
        # test for parameter inter marker distance
        export(pop, format='structure', output='stru.txt')
        self.assertEqual(self.lineOfFile('stru.txt', 2), '-1\t1.0\t-1\t1.0\t1.0\t1.0\n')
        # test for parameter phaseInformation
        self.assertRaises(ValueError, export, pop, format='structure', output='stru.txt',
            phaseInformation=2)
        export(pop, format='structure', output='stru.txt', phaseInformation=0)
        self.assertEqual(self.lineOfFile('stru.txt', 3), '0\n')
        # test for parameter label
        export(pop, format='structure', output='stru.txt', label=False)
        self.assertEqual(self.lineOfFile('stru.txt', 3), '1\t0\t1\t0\t1\t0\t1\n')
        # test for parameter popData
        export(pop, format='structure', output='stru.txt', label=False, popData=False)
        self.assertEqual(self.lineOfFile('stru.txt', 3), '0\t1\t0\t1\t0\t1\n')
        # test for parameter popFlag
        export(pop, format='structure', output='stru.txt', popFlag=1)
        self.assertEqual(self.lineOfFile('stru.txt', 3), '1\t1\t1\t0\t1\t0\t1\t0\t1\n')
        # test for parameter locData
        self.assertRaises(ValueError, export, pop, format='structure', output='stru.txt',
            locData='a')
        export(pop, format='structure', output='stru.txt', locData='loc')
        self.assertEqual(self.lineOfFile('stru.txt', 3), '1\t1\t20\t0\t1\t0\t1\t0\t1\n')
        # test for parameter phenotype
        self.assertRaises(ValueError, export, pop, format='structure', output='stru.txt',
            phenotype='a')
        export(pop, format='structure', output='stru.txt', locData='loc', phenotype='pheno')
        self.assertEqual(self.lineOfFile('stru.txt', 3), '1\t1\t20\t30\t0\t1\t0\t1\t0\t1\n')
        # test for parameter subPops
        self.assertEqual(len(open('stru.txt').read().split('\n')), 21)
        export(pop, format='structure', output='stru.txt', locData='loc',
            subPops=0)
        self.assertEqual(len(open('stru.txt').read().split('\n')), 11)
        # test for virtual subPops
        export(pop, format='structure', output='stru.txt', locData='loc',
            subPops=[(0, 0)])
        self.assertEqual(len(open('stru.txt').read().split('\n')), 7)
        # cleanup
        os.remove('stru.txt')

    def testExportCSV(self):
        '''Testing export population in CSV format'''
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'], infoFields=['loc', 'pheno']) 
        initGenotype(pop, haplotypes=[0,1])
        initInfo(pop, [20], infoFields='loc')
        initInfo(pop, [30], infoFields='pheno')
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        # test basic
        export(pop, format='csv', output='pop.csv')
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex,aff,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), 'M,U,0,0,1,1,0,0,1,1,0,0,1,1\n')
        # test parameter header
        export(pop, format='csv', output='pop.csv', header='header')
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'header\n')
        export(pop, format='csv', output='pop.csv', header=['sex', 'aff', 'header'])
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex,aff,header\n')
        export(pop, format='csv', output='pop.csv', header=False)
        self.assertEqual(self.lineOfFile('pop.csv', 1),'M,U,0,0,1,1,0,0,1,1,0,0,1,1\n')
        # test parameter delimiter
        export(pop, format='csv', output='pop.csv', delimiter='\t')
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex\taff\ta_1\ta_2\tb_1\tb_2\tc_1\tc_2\td_1\td_2\te_1\te_2\tf_1\tf_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), 'M\tU\t0\t0\t1\t1\t0\t0\t1\t1\t0\t0\t1\t1\n')
        # test parameter sexFormatter
        export(pop, format='csv', output='pop.csv', sexFormatter=None)
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'aff,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), 'U,0,0,1,1,0,0,1,1,0,0,1,1\n')
        export(pop, format='csv', output='pop.csv', sexFormatter={MALE: '1', FEMALE: '0'})
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex,aff,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), '1,U,0,0,1,1,0,0,1,1,0,0,1,1\n')
        # test parameter affectionFormatter
        export(pop, format='csv', output='pop.csv', affectionFormatter=None)
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), 'M,0,0,1,1,0,0,1,1,0,0,1,1\n')
        export(pop, format='csv', output='pop.csv', affectionFormatter={True: '1', False: '0'})
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex,aff,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), 'M,0,0,0,1,1,0,0,1,1,0,0,1,1\n')
        # test parameter infoFormatter and infoFields
        export(pop, format='csv', output='pop.csv', infoFields='loc')
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'loc,sex,aff,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), '20.0,M,U,0,0,1,1,0,0,1,1,0,0,1,1\n')
        export(pop, format='csv', output='pop.csv', infoFields='loc', infoFormatter='%.0f')
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'loc,sex,aff,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), '20,M,U,0,0,1,1,0,0,1,1,0,0,1,1\n')
        export(pop, format='csv', output='pop.csv', infoFields=('loc', 'pheno'), infoFormatter='%.0f,%.0f')
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'loc,pheno,sex,aff,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), '20,30,M,U,0,0,1,1,0,0,1,1,0,0,1,1\n')
        # test parameter genoFormatter
        export(pop, format='csv', output='pop.csv', genoFormatter={(0,0):0, (0,1):1, (1,0):1, (1,1):2})
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex,aff,a,b,c,d,e,f\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), 'M,U,0,2,0,2,0,2\n')
        def gen(geno):
            return sum(geno) + 2
        export(pop, format='csv', output='pop.csv', genoFormatter=gen)
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex,aff,a,b,c,d,e,f\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), 'M,U,2,4,2,4,2,4\n')
        # cleanup
        os.remove('pop.csv')
        export(pop, format='csv', output='pop.csv', subPopFormatter='mine')
        self.assertEqual(self.lineOfFile('pop.csv', 1), 'sex,aff,a_1,a_2,b_1,b_2,c_1,c_2,d_1,d_2,e_1,e_2,f_1,f_2,mine\n')
        self.assertEqual(self.lineOfFile('pop.csv', 2), 'M,U,0,0,1,1,0,0,1,1,0,0,1,1,0\n')
        self.assertEqual(self.lineOfFile('pop.csv', 6), 'M,U,0,0,1,1,0,0,1,1,0,0,1,1,1\n')
        # cleanup
        #os.remove('pop.csv')

    def testExportPED(self):
        '''Testing export population in PED format'''
        IdTagger().reset(1)
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'])
        initGenotype(pop, haplotypes=[0,1])
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        # test basic
        export(pop, format='ped', output='pop.ped')
        self.assertEqual(self.lineOfFile('pop.ped', 1), '1 0 0 0 1 1 1 1 2 2 1 1 2 2 1 1 2 2\n')
        self.assertEqual(self.lineOfFile('pop.ped', 2), '2 0 0 0 2 1 1 1 2 2 1 1 2 2 1 1 2 2\n')
        # test parameter subPops
        export(pop, format='ped', output='pop.ped', subPops=[(0,0)])
        self.assertEqual(self.lineOfFile('pop.ped', 1), '1 0 0 0 1 1 1 1 2 2 1 1 2 2 1 1 2 2\n')
        self.assertEqual(self.lineOfFile('pop.ped', 2), '2 0 0 0 1 1 1 1 2 2 1 1 2 2 1 1 2 2\n')
        # test ind_id 
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'], infoFields=['ind_id', 'pheno'])
        initGenotype(pop, haplotypes=[0,1])
        initSex(pop, sex=[MALE, FEMALE])
        tagID(pop)
        pop.setVirtualSplitter(SexSplitter())
        # test basic
        export(pop, format='ped', output='pop.ped')
        self.assertEqual(self.lineOfFile('pop.ped', 1), '1 1 0 0 1 1 1 1 2 2 1 1 2 2 1 1 2 2\n')
        self.assertEqual(self.lineOfFile('pop.ped', 2), '2 2 0 0 2 1 1 1 2 2 1 1 2 2 1 1 2 2\n')
        # test parameter subPops
        export(pop, format='ped', output='pop.ped', subPops=[(0,0)])
        self.assertEqual(self.lineOfFile('pop.ped', 1), '1 1 0 0 1 1 1 1 2 2 1 1 2 2 1 1 2 2\n')
        self.assertEqual(self.lineOfFile('pop.ped', 2), '2 3 0 0 1 1 1 1 2 2 1 1 2 2 1 1 2 2\n')
        # test pheno field
        pop.setIndInfo(range(pop.popSize()), field='pheno')
        export(pop, format='ped', output='pop.ped', subPops=[(0,0)], phenoField='pheno')
        self.assertEqual(self.lineOfFile('pop.ped', 1), '1 1 0 0 1 0.0 1 1 2 2 1 1 2 2 1 1 2 2\n')
        self.assertEqual(self.lineOfFile('pop.ped', 2), '2 3 0 0 1 2.0 1 1 2 2 1 1 2 2 1 1 2 2\n')
        # create a large population
        pop = Population(size=[5000, 20000], ploidy=2, loci=[5,10],
                ancGen=2, infoFields=['fitness', 'father_id', 'mother_id', 'ind_id'])
        pop.evolve(
            initOps=[
                 InitSex(),
                 InitGenotype(freq=[.2, .8], loci=[0]),
                 InitGenotype(freq=[.2]*5, loci=list(range(1, pop.totNumLoci()))),
                 IdTagger(),
            ],
            matingScheme=RandomMating(
                numOffspring=(UNIFORM_DISTRIBUTION, 2, 5),
                ops=[MendelianGenoTransmitter(),
                IdTagger(),
                PedigreeTagger(),
                ]),
            postOps = [
                Stat( alleleFreq=[0,1], genoFreq=[0,1]),
                MapPenetrance(loci=0,
                    penetrance={(0,0):0.1, (0,1):.7, (1,1):1}),
            ],
            gen = 10
        )
        # export pedigree
        s = drawThreeGenFamilySample(pop, families=10, pedSize=(3, 20),
            numOffspring=(1,5), numOfAffected=(0, 5))
        # write in PED format?
        export(s, format='ped', output='pop.ped')
        # get number of families
        # cleanup
        os.remove('pop.ped')

    def testExportMAP(self):
        'Testing export loci information in MAP format'''
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2, 
            lociNames=['a', 'b', 'c', 'd', 'e', 'f'])
        export(pop, format='map', output='pop.map')
        self.assertEqual(self.lineOfFile('pop.map', 1), '1 a 1\n')
        self.assertEqual(self.lineOfFile('pop.map', 3), '2 c 1\n')
        self.assertEqual(self.lineOfFile('pop.map', 6), '2 f 4\n')
        # test no loci name
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2)
        export(pop, format='map', output='pop.map')
        self.assertEqual(self.lineOfFile('pop.map', 1), '1 . 1\n')
        self.assertEqual(self.lineOfFile('pop.map', 3), '2 . 1\n')
        self.assertEqual(self.lineOfFile('pop.map', 6), '2 . 4\n')
        # test float position
        pop = Population(size=[4, 5], loci=[2, 4], ploidy=2,
            lociPos=[0.5, 1.4, 0.001, 1, 2, 4])
        export(pop, format='map', output='pop.map')
        self.assertEqual(self.lineOfFile('pop.map', 1), '1 . 0.5\n')
        self.assertEqual(self.lineOfFile('pop.map', 3), '2 . 0.001\n')
        self.assertEqual(self.lineOfFile('pop.map', 6), '2 . 4\n')
        # test multiplier
        export(pop, format='map', output='pop.map', posMultiplier=1000)
        self.assertEqual(self.lineOfFile('pop.map', 1), '1 . 500\n')
        self.assertEqual(self.lineOfFile('pop.map', 3), '2 . 1\n')
        self.assertEqual(self.lineOfFile('pop.map', 6), '2 . 4000\n')
        # cleanup
        os.remove('pop.map')

    def testExportPhylip(self):
        'Testing export genotype in phylip format'''
        pop = Population(size=[4,5], loci=[20, 90], ploidy=2)
        self.assertRaises(ValueError, export, pop, format='phylip', output='pop.phy')
        initGenotype(pop, freq=[0.25]*4)
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG')
        self.assertEqual(self.lineOfFile('pop.phy', 1), '18 110\n')
        self.assertEqual(len(self.lineOfFile('pop.phy', 2).rstrip()), 100)
        self.assertEqual(len(self.lineOfFile('pop.phy', 3).rstrip()), 20)
        self.assertEqual(len(self.lineOfFile('pop.phy', 4).rstrip()), 100)
        # no seqname, diploid
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'S2_1      ')
        self.assertRaises(ValueError, export, pop, format='phylip', output='pop.phy', seqNames=['1', '2'])
        # individual names
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', seqNames=['N%d' % x for x in range(9)])
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'N1_1      ')
        # seq names
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', seqNames=['N%d' % x for x in range(18)])
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'N2        ')
        # haploid
        pop = Population(size=[4,5], loci=[20, 90], ploidy=1)
        self.assertRaises(ValueError, export, pop, format='phylip', output='pop.phy')
        initGenotype(pop, freq=[0.25]*4)
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG')
        self.assertEqual(self.lineOfFile('pop.phy', 1), '9 110\n')
        self.assertEqual(len(self.lineOfFile('pop.phy', 2).rstrip()), 100)
        self.assertEqual(len(self.lineOfFile('pop.phy', 3).rstrip()), 20)
        self.assertEqual(len(self.lineOfFile('pop.phy', 4).rstrip()), 100)
        # no seqname, diploid
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'S3        ')
        # individual names
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', seqNames=['N%d' % x for x in range(9)])
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'N2        ')
        # 
        # another style (interleaved)
        #
        pop = Population(size=[4,5], loci=[20, 90], ploidy=2)
        self.assertRaises(ValueError, export, pop, format='phylip', output='pop.phy', alleleNames='ACTG', style='asd')
        initGenotype(pop, freq=[0.25]*4)
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', style='interleaved')
        self.assertEqual(self.lineOfFile('pop.phy', 1), '18 110\n')
        self.assertEqual(len(self.lineOfFile('pop.phy', 2).rstrip()), 100)
        self.assertEqual(len(self.lineOfFile('pop.phy', 3).rstrip()), 100)
        self.assertEqual(len(self.lineOfFile('pop.phy', 4).rstrip()), 100)
        self.assertEqual(self.lineOfFile('pop.phy', 20), '\n')
        self.assertEqual(len(self.lineOfFile('pop.phy', 21).rstrip()), 20)
        # no seqname, diploid
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'S3_1      ')
        self.assertRaises(ValueError, export, pop, format='phylip', output='pop.phy', seqNames=['1', '2'])
        # individual names
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', seqNames=['N%d' % x for x in range(9)], style='interleaved')
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'N2_1      ')
        # seq names
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', seqNames=['N%d' % x for x in range(18)], style='interleaved')
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'N4        ')
        # haploid
        pop = Population(size=[4,5], loci=[20, 90], ploidy=1)
        initGenotype(pop, freq=[0.25]*4)
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', style='interleaved')
        self.assertEqual(self.lineOfFile('pop.phy', 1), '9 110\n')
        self.assertEqual(len(self.lineOfFile('pop.phy', 2).rstrip()), 100)
        self.assertEqual(len(self.lineOfFile('pop.phy', 3).rstrip()), 100)
        self.assertEqual(len(self.lineOfFile('pop.phy', 4).rstrip()), 100)
        self.assertEqual(len(self.lineOfFile('pop.phy', 11).rstrip()), 0)
        self.assertEqual(len(self.lineOfFile('pop.phy', 12).rstrip()), 20)
        # no seqname, diploid
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'S5        ')
        # individual names
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', seqNames=['N%d' % x for x in range(9)], style='interleaved')
        self.assertEqual(self.lineOfFile('pop.phy', 6)[:10], 'N4        ')
        # cleanup
        os.remove('pop.phy')

    def testImportPhylip(self):
        'Testing export genotype in phylip format'''
        pop = Population(size=[4,5], loci=[20, 90], ploidy=2)
        initGenotype(pop, freq=[0.25]*4)
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG')
        pop1 = importPopulation(format='phylip', filename='pop.phy', alleleNames='ACTG')
        self.assertEqual(pop.genotype(), pop1.genotype())
        self.assertEqual(pop1.ploidy(), 1)
        #
        pop1 = importPopulation(format='phylip', filename='pop.phy', alleleNames='ACTG', ploidy=2)
        self.assertEqual(pop.genotype(), pop1.genotype())
        self.assertEqual(pop1.ploidy(), 2)
        #
        # another style (interleaved)
        #
        pop = Population(size=[4,5], loci=[20, 90], ploidy=2)
        initGenotype(pop, freq=[0.25]*4)
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', style='interleaved')
        pop1 = importPopulation(format='phylip', filename='pop.phy', alleleNames='ACTG')
        self.assertEqual(pop.genotype(), pop1.genotype())
        self.assertEqual(pop1.ploidy(), 1)
        #
        pop1 = importPopulation(format='phylip', filename='pop.phy', alleleNames='ACTG', ploidy=2)
        self.assertEqual(pop.genotype(), pop1.genotype())
        self.assertEqual(pop1.ploidy(), 2)
        #
        # starting from a haploid population
        pop = Population(size=[4,5], loci=[20, 90], ploidy=1)
        initGenotype(pop, freq=[0.25]*4)
        export(pop, format='phylip', output='pop.phy', alleleNames='ACTG', style='interleaved')
        pop1 = importPopulation(format='phylip', filename='pop.phy', alleleNames='ACTG')
        self.assertEqual(pop.genotype(), pop1.genotype())
        self.assertEqual(pop1.ploidy(), 1)
        # import as diploid?
        self.assertRaises(ValueError, importPopulation, format='phylip', filename='pop.phy',
            alleleNames='ACTG', ploidy=2)
        # cleanup
        os.remove('pop.phy')

    def pyFunc(self, line):
        pass

    def testExportToFunc(self):
        'Testing export genotype to python func, to file handle, and to gzip'
        # starting from a haploid population
        pop = Population(size=[4,5], loci=[20, 90], ploidy=1)
        initGenotype(pop, freq=[0.25]*4)
        # 1: export to a python func
        export(pop, format='phylip', output=self.pyFunc, alleleNames='ACTG', style='interleaved')
        # 2: export to a file handle
        with open('pop.phy', 'w') as phy:
            export(pop, format='phylip', output=phy, alleleNames='ACTG', style='interleaved')
        pop1 = importPopulation(format='phylip', filename='pop.phy', alleleNames='ACTG')
        self.assertEqual(pop.genotype(), pop1.genotype())
        self.assertEqual(pop1.ploidy(), 1)
        # 3: export to a gzip file handle
        with gzip.open('pop.phy.gz', 'wb') as phy:
            export(pop, format='phylip', output=WithMode(phy, 'b'), alleleNames='ACTG', style='interleaved')
        # unzip first?
        f_in = gzip.open('pop.phy.gz', 'rb')
        f_out = open('pop1.phy', 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        #
        pop1 = importPopulation(format='phylip', filename='pop1.phy', alleleNames='ACTG')
        self.assertEqual(pop.genotype(), pop1.genotype())
        self.assertEqual(pop1.ploidy(), 1)
        #
        # cleanup
        os.remove('pop.phy')
        os.remove('pop1.phy')
        os.remove('pop.phy.gz')

    def testExportDyncOutput(self):
        'Testing the Exporter operator with output starting with !'
        pop = Population(100)
        pop.evolve(
            initOps=[InitSex()],
            preOps=
                Exporter(format='CSV', output='!"output%d.csv" % gen'),
            matingScheme=RandomMating(),
            gen=5
        )
        for i in range(5):
            self.assertTrue(os.path.isfile('output{}.csv'.format(i)))
            os.remove('output{}.csv'.format(i))

if __name__ == '__main__':
    unittest.main()
