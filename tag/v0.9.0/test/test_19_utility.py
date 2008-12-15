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
import unittest, os, sys, exceptions

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
        for rg in ListAllRNG():
            SetRNG(rg)

    def testDefaultRNG(self):
        'Testing default RNG'
        rg = rng()
        self.assertEqual(rg.name(), 'mt19937')

    def testBinomial(self):
        'Testing binomial distribution'
        rg = rng()
        for n in range(1,10):
            rg.randBinomial(10, .7)

    def testUniform01(self):
        'Testing uniform distribution generator'
        rg = rng()
        for n in range(1,10):
            rg.randUniform01()

    def testRandBit(self):
        'Testing random bit function'
        rg = rng()
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
        rg = rng()
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
        'repeated set rng() without seed, and see if the seed repeat'
        import random
        seed = []
        name = ListAllRNG()[0]
        for i in range(100):
            SetRNG(name)
            sd = rng().seed()
            assert sd not in seed
            seed.append(sd)
        # test set seed
        sd = random.randint(100, 10000)
        SetRNG(name, sd)
        self.assertEqual(rng().seed(), sd)


    def testWeightedSampler(self):
        'Testing weighted sampler'
        sampler = Weightedsampler(rng(), [1, 2, 3, 4])
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


    def TestLeakLoadPopulation(self):
        'Testing if loadPopulation leaks'
        # run this and see if memory usage goes up continuously
        pop = population(100, loci=[1000]*10)
        Stat(pop, alleleFreq=range(pop.totNumLoci()))
        pop.savePopulation('test.bin')
        for i in range(10000):
            print 'Loading %i' % i
            pop = LoadPopulation('test.bin')
            print pop.dvars().alleleFreq[100][0]


if __name__ == '__main__':
    unittest.main()
