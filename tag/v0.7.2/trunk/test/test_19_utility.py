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

class TestRNG(unittest.TestCase):

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
            prop = bt.succRate(i)
            # binomial, mean p, variance = p(1-p)/n
            std = math.sqrt(p[i]*(1.-p[i])/N)
            assert prop > p[i] - 3*std and prop < p[i] + 3*std
        # another test, for each trail
        for pi in p:
            bt = BernulliTrials(rg, [pi]*N, 10)
            bt.doTrial()
            for i in range(10):
                bt.trial();
                prop = bt.trialRate()
                # binomial, mean p, variance = p(1-p)/n
                std = math.sqrt(pi*(1.-pi)/N)
                assert prop > pi - 3*std and prop < pi + 3*std


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
            
    
if __name__ == '__main__':
    unittest.main()