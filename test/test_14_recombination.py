#!/usr/bin/env python
#
#  This is a unittest file for operator recombinator
#
#  Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
# 

import simuOpt
simuOpt.setOptions(quiet=False)

from simuPOP import *
import unittest, os, sys, exceptions

class TestRecombinator(unittest.TestCase):
    
    def testRecRate(self):
        'Testing to see if we actually recombine at this rate '
        a1, a2 = 0, 1
            
        pop = population(10000, loci=[2,3,2])
        InitByValue(pop, value=[a1]*7+[a2]*7)
        simu = simulator(pop, randomMating())
        simu.step( [ 
            stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]), 
            recombinator(rate = 0.1) ] )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        #print simu.dvars(0).haploFreq
        assert abs(simu.dvars(0).haploFreq['0-1']['%s-%s'%(a1,a2)] - 0.05) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['%s-%s'%(a1,a2)] - 0.05) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['%s-%s'%(a1,a2)] - 0.05) < 0.01
        assert abs(simu.dvars(0).haploFreq['4-5']['%s-%s'%(a1,a2)] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['%s-%s'%(a1,a2)] - 0.05) < 0.01
        # compare to the next test
        pop = population(10000, loci=[3,4])
        InitByValue(pop, value=[a1]*7+[a2]*7)
        simu = simulator(pop, randomMating())
        simu.step( [ 
            stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]), 
            recombinator(rate = 0.4, loci=[1,3]) ] )
        # 0.5
        assert abs(simu.dvars(0).haploFreq['0-1']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['0-1']['1-1'] - 0.5) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['1-2']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-1'] - 0.30) < 0.01
        # 0.25
        assert abs(simu.dvars(0).haploFreq['2-3']['0-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['0-1'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-1'] - 0.25) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['3-4']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-1'] - 0.30) < 0.01
        # no recombination.
        assert abs(simu.dvars(0).haploFreq['4-5']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['0-0'] - 0.5) < 0.01
        #
        # alrorithm 0?
        # 
        pop = population(10000, loci=[3,10])
        InitByValue(pop, value=[a1]*13+[a2]*13)
        simu = simulator(pop, randomMating())
        simu.step( [ 
            stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]), 
            recombinator(rate = 0.4, loci=[1,3,8]) ] )
        # 0.5
        assert abs(simu.dvars(0).haploFreq['0-1']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['0-1']['1-1'] - 0.5) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['1-2']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-1'] - 0.30) < 0.01
        # 0.25
        assert abs(simu.dvars(0).haploFreq['2-3']['0-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['0-1'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-1'] - 0.25) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['3-4']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-1'] - 0.30) < 0.01
        # no recombination.
        assert abs(simu.dvars(0).haploFreq['4-5']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['0-0'] - 0.5) < 0.01

    def testConversionRate(self):
        'Testing to see if we actually convert at this rate '
        a1, a2 = 0, 1
        pop = population(10000, loci=[3,4])
        InitByValue(pop, value=[a1]*7+[a2]*7)
        simu = simulator(pop, randomMating())
        rec = recombinator(rate = 0.4, convProb=1, loci=[1,3], convParam=1)
        simu.step( [
            stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]), 
            rec ] )
        #print simu.dvars(0).haploFreq
        #print rec.convCounts()
        # 0.5
        assert abs(simu.dvars(0).haploFreq['0-1']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['0-1']['1-1'] - 0.5) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['1-2']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-1'] - 0.30) < 0.01
        # 0.25
        assert abs(simu.dvars(0).haploFreq['2-3']['0-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['0-1'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-1'] - 0.25) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['3-4']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-1'] - 0.30) < 0.01
        # r/2 recombination induced recombination
        assert abs(simu.dvars(0).haploFreq['4-5']['0-0'] - 0.3) < 0.01
        assert abs(simu.dvars(0).haploFreq['4-5']['0-1'] - 0.2) < 0.01
        assert abs(simu.dvars(0).haploFreq['4-5']['1-1'] - 0.3) < 0.01
        assert abs(simu.dvars(0).haploFreq['4-5']['1-0'] - 0.2) < 0.01
        # copied
        assert abs(simu.dvars(0).haploFreq['5-6']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['1-1'] - 0.5) < 0.01
        #
        # length 2
        pop = population(10000, loci=[3,4])
        InitByValue(pop, value=[a1]*7+[a2]*7)
        simu = simulator(pop, randomMating())
        simu.step( [ 
            stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]), 
            recombinator(rate = 0.4, convProb=1, convParam=2, loci=[1,3]) ] )
        #
        assert abs(simu.dvars(0).haploFreq['0-1']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['0-1']['1-1'] - 0.5) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['1-2']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-1'] - 0.30) < 0.01
        # 0.25
        assert abs(simu.dvars(0).haploFreq['2-3']['0-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['0-1'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-1'] - 0.25) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['3-4']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-1'] - 0.30) < 0.01
        # copied
        assert abs(simu.dvars(0).haploFreq['4-5']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['4-5']['1-1'] - 0.5) < 0.01
        # r/2 recombination induced recombination
        assert abs(simu.dvars(0).haploFreq['5-6']['0-0'] - 0.3) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['0-1'] - 0.2) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['1-1'] - 0.3) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['1-0'] - 0.2) < 0.01
        #
        # algorithm 0??
        pop = population(10000, loci=[3,10])
        InitByValue(pop, value=[a1]*13+[a2]*13)
        simu = simulator(pop, randomMating())
        simu.step( [ 
            stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]), 
            recombinator(rate = 0.4, convProb=1, convParam=2, loci=[1,3,8]) ] )
        #
        assert abs(simu.dvars(0).haploFreq['0-1']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['0-1']['1-1'] - 0.5) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['1-2']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['1-2']['1-1'] - 0.30) < 0.01
        # 0.25
        assert abs(simu.dvars(0).haploFreq['2-3']['0-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['0-1'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-0'] - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['2-3']['1-1'] - 0.25) < 0.01
        # r/2
        assert abs(simu.dvars(0).haploFreq['3-4']['0-1'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-0'] - 0.20) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['0-0'] - 0.30) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-4']['1-1'] - 0.30) < 0.01
        # copied
        assert abs(simu.dvars(0).haploFreq['4-5']['0-0'] - 0.5) < 0.01
        assert abs(simu.dvars(0).haploFreq['4-5']['1-1'] - 0.5) < 0.01
        # r/2 recombination induced recombination
        assert abs(simu.dvars(0).haploFreq['5-6']['0-0'] - 0.3) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['0-1'] - 0.2) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['1-1'] - 0.3) < 0.01
        assert abs(simu.dvars(0).haploFreq['5-6']['1-0'] - 0.2) < 0.01

                
    def testAtLociRecRates(self):
        'Testing loci parameter'
        if AlleleType() == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        pop = population(10000, loci=[2,3,2])
        InitByValue(pop, value=[a1]*7+[a2]*7)
        simu = simulator(pop, randomMating())
        simu.step( [ 
            stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]), 
            recombinator(rate = [0.1,0.15,0.3], afterLoci=[0,2,5] ) ] )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        #print simu.dvars(0).haploFreq
        assert (simu.dvars(0).haploFreq['0-1']['%s-%s'%(a1,a2)] - 0.05) < 0.01
        assert (simu.dvars(0).haploFreq['2-3']['%s-%s'%(a1,a2)] - 0.075)< 0.01
        try:    # do not have this haplotype
            simu.dvars(0).haploFreq['3-4']['%s-%s'%(a1,a2)]
        except exceptions.KeyError:
            pass
        assert (simu.dvars(0).haploFreq['4-5']['%s-%s'%(a1,a2)] - 0.25) < 0.01
        assert (simu.dvars(0).haploFreq['5-6']['%s-%s'%(a1,a2)] - 0.15) < 0.01
                    
    def testRecIntensity(self):
        'Testing recombination intensity'
        if AlleleType() == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        pop = population(10000, loci=[2,3,2], lociPos=[0,1,0,2,4,0,4] )
        InitByValue(pop, value=[a1]*7+[a2]*7)
        simu = simulator(pop, randomMating())
        simu.step( [ 
            stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]), 
            recombinator(intensity = 0.1) ] )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        #print simu.dvars(0).haploFreq
        assert (simu.dvars(0).haploFreq['0-1']['%s-%s'%(a1,a2)] - 0.05) < 0.01
        assert (simu.dvars(0).haploFreq['2-3']['%s-%s'%(a1,a2)] - 0.1)    < 0.01
        assert (simu.dvars(0).haploFreq['3-4']['%s-%s'%(a1,a2)] - 0.1)    < 0.01
        assert (simu.dvars(0).haploFreq['4-5']['%s-%s'%(a1,a2)] - 0.25) < 0.01
        assert (simu.dvars(0).haploFreq['5-6']['%s-%s'%(a1,a2)] - 0.2)    < 0.01
    

    def testRecIntensityAfterLoci(self):
        'Testing RecIntensity after loci'
        if AlleleType() == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        pop = population(10000, loci=[2,3,2], lociPos=[0,1,0,2,4,0,4] )
        InitByValue(pop, value=[a1]*7+[a2]*7)
        simu = simulator(pop, randomMating())
        simu.step( [ 
            stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]), 
            recombinator(intensity = 0.1, afterLoci=[2,5]) ] )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        #print simu.dvars(0).haploFreq
        assert simu.dvars(0).haploFreq['0-1'].setdefault('%s-%s'%(a1,a2),0) == 0
        assert (simu.dvars(0).haploFreq['2-3']['%s-%s'%(a1,a2)] - 0.1)    < 0.01
        assert simu.dvars(0).haploFreq['3-4'].setdefault('%s-%s'%(a1,a2),0) == 0
        assert (simu.dvars(0).haploFreq['4-5']['%s-%s'%(a1,a2)] - 0.25) < 0.01
        assert (simu.dvars(0).haploFreq['5-6']['%s-%s'%(a1,a2)] - 0.2) < 0.01        


    def testNoNewAllele(self):
        'Testing that no new allele appear because of recombination'
        if AlleleType() == 'binary':
            geno = [x%2 for x in [1,2,3,4,5,6,7] ]
        else:
            geno = [1,2,3,4,5,6,7]
        pop = population(1000, loci=[2,3,2])
        InitByValue(pop, value=geno)
        simu = simulator(pop, randomMating())
        simu.evolve( [ recombinator(rate = 0.4) ], end=100)
        Stat(simu.population(0), alleleFreq=range(0,7))
        for i in range(7):
            self.assertEqual(simu.dvars(0).alleleFreq[i][geno[i]], 1.)


    def testCrossBetweenChrom(self):
        'Testing if chromsomes are crossed by default'
        if AlleleType() == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        pop = population(100000, loci=[2,3,2])
        InitByValue(pop, value=[a1]*7+[a2]*7)
        simu = simulator(pop, randomMating())
        #TurnOnDebug(DBG_RECOMBINATOR)
        simu.step( [ 
            stat( haploFreq = [[0,1], [2,3], [2,4], [5,6], [0,2], [0,6], [3,6]],
                stage=PrePostMating), 
            ## for debug purpose, output haploFreq
            #pyEval('haploFreq', stage=PrePostMating),
            recombinator(rate = 0) ] )
        self.assertEqual(simu.dvars(0).haploFreq['0-1'].setdefault('%s-%s'%(a1,a2), 0), 0.)
        self.assertEqual(simu.dvars(0).haploFreq['2-3'].setdefault('%s-%s'%(a1,a2), 0), 0.)
        self.assertEqual(simu.dvars(0).haploFreq['2-4'].setdefault('%s-%s'%(a1,a2), 0), 0.)
        self.assertEqual(simu.dvars(0).haploFreq['5-6'].setdefault('%s-%s'%(a1,a2), 0), 0.)
        assert abs(simu.dvars(0).haploFreq['0-2'].setdefault('%s-%s'%(a1,a2), 0) - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['0-6'].setdefault('%s-%s'%(a1,a2), 0) - 0.25) < 0.01
        assert abs(simu.dvars(0).haploFreq['3-6'].setdefault('%s-%s'%(a1,a2), 0) - 0.25) < 0.01


    def testRecProportion(self):
        'Testing table 4 of H&C 3nd edition P49 '
        N = 10000
        r = 0.1
        genoDad = [[0,0],[0,0],[0,0],[0,0],[0,1],[0,1],[0,1],[1,0],[1,0],[1,1]]
        genoMom = [[0,0],[0,1],[1,0],[1,1],[0,1],[1,0],[1,1],[1,0],[1,1],[1,1]]
        prop = [ [1, 0, 0, 0], [.5, .5, 0, 0], [.5, 0, 0.5, 0],
            [0.5-r/2, r/2, r/2, 0.5-r/2], [0, 1, 0, 0], [r/2, .5-r/2, .5-r/2, r/2], 
            [0, .5, 0, .5], [0, 0, 1, 0], [0, 0, .5, .5], [0, 0, 0, 1] ]
        for i in range(0, len(genoDad)):
            pop = population(size=N, loci=[2])
            InitByValue(pop, value=genoDad[i]+genoMom[i])
            simu = simulator(pop, randomMating())
            simu.step(
                    [ recombinator(rate=r),
                        stat(haploFreq=[0,1]) 
                    ])
            hf = simu.dvars(0).haploFreq['0-1']
            assert not ((hf.setdefault('0-0',0) - prop[i][0]) > 0.01 or \
                (hf.setdefault('0-1',0) - prop[i][1]) > 0.01 or \
                (hf.setdefault('1-0',0) - prop[i][2]) > 0.01 or \
                (hf.setdefault('1-1',0) - prop[i][3]) > 0.01), \
                "Recombination results in potentially wrong proportions." +    \
                str( genoDad[i]) + ' crossing ' + str(genoMom[i])


    def testLDDecay(self):
        'Testing formula Dn=(1-r)^n D0 '
        r = 0.1
        N = 10000
        pop = population(size=N, loci=[2])
        # genotype 11/22, with LD=1
        if AlleleType() == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        InitByValue(pop, value=[a1,a1,a2,a2])
        simu = simulator(pop, randomMating())
        simu.evolve( 
            ops    = [ recombinator(rate=r) ],
            postOps = [    stat(LD=[0,1]) ], 
            end=9)
        # check the change of LD, hopefully, the variation is not too high.
        assert abs(simu.dvars(0).LD[0][1] - 0.25*(1-r)**10) < 0.02, \
            "Decay of LD is not as expected: " + str(simu.dvars(0).LD[0][1]) + " vs expected " \
            + str( 0.25*(1-r)**10 )

    def testRecCount(self):
        'Testing count of recombination events '
        r = 0.01
        N = 1000
        G = 100
        rec = recombinator(rate=r)
        pop = population(size=N, loci=[10,10])
        simu = simulator(pop, randomMating())
        simu.evolve( [ rec ], end=G)
        # number of recombination event should be bionomial(ploidy*N*r, 0.1) with mean 10000
        # at end should be bionomial(ploidy*N*r, 0.5)
        assert abs( rec.recCount(0) - 2*N*r*G ) < 150, \
            'Number of recombination event is not as expected. %d %d' % (rec.recCount(0), 2*N*r*G)
        assert abs( rec.recCount(9) - 2*N*0.5*G ) < 2000, \
            'Number of recombination event is not as expected. %d %d' % (rec.recCount(9), 2*N*0.5*G)
        

    def testNoMaleRec(self):
        'Testing recombination of male chromosome'
        # create such an situation where female has 11111/22222, male has 1111/33333
        # we will test if 3333 is untouched under recombination
        r = 0.1
        N = 100
        if AlleleType() == 'binary':
            a1, a2, a3 = 0, 1, 1
        else:
            a1, a2, a3 = 1, 2, 3
        pop = population(size=N, loci=[2,5], sexChrom=True)
        # male     1 3
        # female 1 2
        InitByValue(pop, indRange=[0,N/2], sex=[Male]*(N/2), atPloidy=0, value=[a1]*7)
        InitByValue(pop, indRange=[0,N/2], sex=[Male]*(N/2), atPloidy=1, value=[a1]*2+[a3]*5)
        InitByValue(pop, indRange=[N/2,N], sex=[Female]*(N/2), value=[a1]*7+[a2]*7)
        # now let us recombine
        simu = simulator(pop, randomMating())
        simu.evolve( [ recombinator(rate=r) ], end=100)
        pop = simu.population(0)
        #
        for i in range( pop.popSize() ):
            ind = pop.individual(i)
            if ind.sex() == Male:
                # check the second chromosome
                # arrGenotype(ploidy, chrom), no dict parameter for efficiency purpose
                #print ind.arrGenotype()
                self.assertEqual(ind.arrGenotype(1, 1), [a3]*5)
            else: 
                # there is no allele 3 anywhere (only non-binary...)
                if AlleleType() != 'binary':
                    self.assertEqual(ind.arrGenotype().count(a3), 0)


    def testConversionCount(self):
        'Testing count of conversion events '
        r = 0.01
        N = 1000
        G = 100
        for mode, param in \
                [(CONVERT_GeometricDistribution, 0.3),
                 (CONVERT_TractLength, 2.5),
                 (CONVERT_NumMarkers, 2),
                 (CONVERT_ExponentialDistribution, 1)]:
            rec = recombinator(rate=r, convProb=0.2, 
                convMode=mode, convParam=param)
            pop = population(size=N, loci=[10,10])
            simu = simulator(pop, randomMating())
            simu.evolve( [ rec ], end=G-1)
            # at end should be bionomial(ploidy*N*r, 0.5)
            recCount = sum(rec.recCounts()) - max(rec.recCounts())
            convCount = sum(rec.convCounts().values())
            assert abs(recCount*0.2 - convCount) < 100
        


if __name__ == '__main__':
    unittest.main()   
