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
            randomMating(), rep=5)
        print "\n\nUSER INTERACTION: Please press q\n\n"
        self.assertRaises( exceptions.SystemError, simu.evolve,
            ops=[ pause(at=[10]),
                        # should quite, can not reach generation 12
                        terminateIf("True", at=[12] ) ] )

    def interactiveTestExitToShell(self):
        'Testing exit to a shell'
        simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
            randomMating(), rep=5)
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

    def testLogging(self):
        'Testing logging output to a file'
        setLogOutput("session.log")
        # output stuff
        setLogOutput()
        os.remove('session.log')


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
        pop = population(100, loci=[1000]*10)
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

    # Test Backward and Forward Trajectories
    def testSimuCase1(self):
        # 1: independent case when freq = 0 in backward simulation
        # and [0, 1, 0] in forward simulation. nLoci = 3, constant nSubPop = 3:
        def Nt(gen):
            return [1000] * 3
        def fitness(gen):
            return [1] * 3
        trajSimulator = trajectorySimulator(N=Nt, fitness=fitness, nLoci = 3)
        traj = trajSimulator.simuBackward(genEnd=3000, freq=0, minMutAge = 0, maxMutAge = 100000, ploidy = 2,
                     restartIfFail = False, maxAttempts=1000, logger=None)
        self.assertEqual(traj.traj[min(traj.traj.keys())], [0] * 9)
        traj1 = trajSimulator.simuForward(freq = [0, 1, 0], destFreq = [[0,1]]*3, genEnd = 100)
        self.assertEqual(traj1.freq(random.randint(0, 100)), [0,0,0,1,1,1,0,0,0])        

    def testSimuCase2(self):
        # 2: dependent case when freq = 0 in backward simulation
        # and [0, 1] in forward simulation. nLoci = 2, constant nSubPop = 3:
        def Nt(gen):
            return [2000] * 3
        def fitness(gen):
            return [1] * 9
        trajSimulator = trajectorySimulator(N=Nt, fitness=fitness, nLoci = 2)
        traj = trajSimulator.simuBackward(genEnd=3000, freq=0)
        self.assertEqual(traj.traj[min(traj.traj.keys())], [0] * 6)
        traj1 = trajSimulator.simuForward(freq = [0, 1], destFreq = [[0, 1]]*2, genEnd = 200)
        self.assertEqual(traj1.freq(random.randint(0, 200)), [0,0,0,1,1,1])
        
    def testSimuCase3(self):
        # 3: nLoci = 2, fitness = [1, 0.0001, 0.0001, 1, 0.0002, 0.0002], freq = [0.5, 0.3],
        # no subpopulation. If fitness for disease alleles are close to 0, backward simulation 
        # should reach maxAttempts and freq should hit 0 before genEnd in forward simulation.
        fitness = [1, 0.0001, 0.0001, 1, 0.0002, 0.0002]
        trajSimulator = trajectorySimulator(N=1000, fitness=fitness, nLoci = 2)
        traj = trajSimulator.simuBackward(genEnd=3000, freq=[0.5, 0.3], maxAttempts=100)
        self.assertEqual(len(traj.traj), 0)
        self.assertEqual(traj.freq(0), [0] * 2)
        traj1 = trajSimulator.simuForward(genEnd = 200, freq = [0.2, 0.4],
                                          destFreq = [[0, 1]] * 2)
        self.assertEqual(traj1.freq(200), [0, 0])
        
    def testSimuCase4(self):
        # 4: In backward simulation, given a normal set of parameters, if current freq of
        # allele a is very small, the first element of traj[] should be 1/(ploidy * N) and
        # the second element of traj[] should be larger than 1/(ploidy * N):
        # nLoci = 1 and no subpopulation
        fitness = [1,1,1]
        #### ??? N = 101
        trajSimulator = trajectorySimulator(N=1000, fitness=fitness, nLoci = 1)
        traj = trajSimulator.simuBackward(genEnd=3000, freq=0.01)
        #####
        #print len(traj.traj)
        #print [traj.traj[min(traj.traj.keys()) + i] for i in range(3)]
        #####
        self.assertEqual(traj.traj[min(traj.traj.keys())+1][0], 1. / (2 * 1000))
        if traj.traj[min(traj.traj.keys())+2][0] < 1. / (2 * 1000):
            raise ValueError('fail in test number 4 of testSimuBackward.')
        
    def testSimuCase5(self):
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
        traj = trajSimulator.simuBackward(genEnd=300, freq=0.01)
        self.assertEqual(traj.traj[min(traj.traj.keys())], [0] * 3)
        self.assertEqual(len(traj.freq(296)), 9)
        self.assertEqual(len(traj.freq(295)), 3)
        traj1 = trajSimulator.simuForward(genEnd=300, freq=0.5, destFreq=[[0,1]]*3)
        self.assertEqual(len(traj1.freq(295)), 3)
        self.assertEqual(len(traj1.freq(296)), 9)
        
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
        traj = trajSimulator.simuBackward(genEnd=300, freq=0.01)
        self.assertEqual(traj.traj[min(traj.traj.keys())], [0] * 12)
        self.assertEqual(len(traj.freq(296)), 3)
        self.assertEqual(len(traj.freq(295)), 12)
        traj1 = trajSimulator.simuForward(genEnd=300, freq=0.5, destFreq=[[0,1]]*3)
        self.assertEqual(len(traj1.freq(295)), 12)
        self.assertEqual(len(traj1.freq(296)), 3)

    def testSimuCase6(self):
        # 6: given a normal set of parameters, plot a backward trajectory and
        # another forward trajectory for the case of single locus without
        # subpopulations.
        trajSimulator = trajectorySimulator(N=3000, fitness=[1,1,1], nLoci = 1)
        traj = trajSimulator.simuBackward(genEnd=3000, freq=0.1)
        #traj.plot()
        traj1 = trajSimulator.simuForward(genEnd=200, freq=0.5, destFreq=[[0,1]])
        #traj1.plot()
        
    def testSimuCase7(self):
        # 7: given a normal set of parameters, considering changable subpopulation
        # sizes with mulitiple loci, plot a backward trajectory and another
        # forward trajectory. nSubPops = 3, nLoci = 2.
        def Nt(gen):
            if gen > 2500:
                return [1000] * 3
            else:   # merge from all subpops into one population in backward sense
                return 3000
        trajSimulator = trajectorySimulator(N=Nt, fitness=[1,1,1], nLoci = 2)
        traj = trajSimulator.simuBackward(genEnd=3000, freq=[0.05, 0.1])
        #traj.plot()
        def Nt1(gen):
            if gen < 100:
                return 3000
            else:   # split from one population into 3 subpops in forward sense
                return [1000] * 3
        trajSimulator1 = trajectorySimulator(N=Nt1, fitness=[1,1,1], nLoci = 5)
        traj1 = trajSimulator1.simuForward(genEnd = 200, freq = [0.5, 0.6, 0.7, 0.8, 0.9],
                                           destFreq = [[0,1]]*5)
        #traj1.plot()
            

if __name__ == '__main__':
    unittest.main()
    

if __name__ == '__main__':
    unittest.main()
