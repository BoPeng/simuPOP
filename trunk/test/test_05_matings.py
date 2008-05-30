#!/usr/bin/env python
#
# unittests for mating schemes
# 
# Author:
#   Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, random, math, sets

def setGen(pop, off, dad, mom):
    off.setAllele(pop.gen(), 0)
    return True
        
class TestMatingSchemes(unittest.TestCase):

    def testNoMating(self):
        'Testing noMating mating scheme'
        simu = simulator(population(10, loci=[1], ploidy=1),
            noMating())
        # during mating operator will be applied
        simu.step( ops=[ pyOperator(func=setGen, stage=DuringMating) ])
        self.assertEqual( simu.population(0).arrGenotype(True),
            [0]*10)
        simu.step( ops=[ pyOperator(func=setGen, stage=DuringMating) ])
        self.assertEqual( simu.population(0).arrGenotype(True),
            [1]*10)

    def testBinomialSelection(self):
        'Testing binomialSelection mating scheme (FIXME: imcomplete)'
        simu = simulator(population(10, loci=[1], ploidy=1),
            binomialSelection())
        
    def testSelection(self):
        'Testing selections (FIXME: imcomplete)'
        pass
     
    def testPopSizeChange(self):
        'Testing means to change population size (FIXME: imcomplete)'
        pass
        
    def getFamSize(self, mate, endGen=1, size=1000):
        simu = simulator(population(size, loci=[1]), mate)
        simu.evolve(ops=[], gen=endGen)
        return simu.population(0).dvars().famSizes
        
    def testNumOffspring(self):
        'Testing means to control number of offspring (FIXME: check distribution)'
        TurnOnDebug(DBG_MATING)
        self.assertEqual( 
            self.getFamSize(binomialSelection(numOffspring=2)),
            [2]*500)
        # numOffspringFunc
        def nos(gen):
            return gen%2+1
        self.assertEqual( 
            self.getFamSize(binomialSelection(numOffspringFunc=nos), endGen=2),
            [2]*500)
        self.assertEqual( 
            self.getFamSize(binomialSelection(numOffspringFunc=nos), endGen=3),
            [1]*1000)
        # what if each family have different number of offspring?
        def nos(gen):
            return random.randrange(1,4)
        #
        cnt = self.getFamSize(randomMating(numOffspringFunc=nos))
        self.assertEqual( sum(cnt), 1000)
        num = [ cnt.count(i) for i in range(1,4) ]
        # test for uniform?
        mean = sum(num)/3.
        for i in range(3):
            assert num[i] < mean + 50 and num[i] > mean - 50
        #
        # MATE_GeometricDistribution
        cnt = self.getFamSize( randomMating(numOffspring=.3, 
                mode=MATE_GeometricDistribution))
        #print cnt    
        # MATE_BinomialDistribution
        cnt = self.getFamSize( randomMating(numOffspring=.3, 
            maxNumOffspring=5, mode=MATE_BinomialDistribution))
        #print cnt
        # MATE_PoissonDistribution
        cnt = self.getFamSize( randomMating(numOffspring=.3, 
            mode=MATE_PoissonDistribution))
        #print cnt
        TurnOffDebug(DBG_MATING)
        
##     def testTrajectory(self):
##         'Testing trajectory prediction functions'
##         sel, Ne, freq, h, selection = 0.5, 100, 0.50, 2, 1
##         path = FreqTrajectorySelSim(sel, Ne, freq, h, selection)
##         # the second method, forward, with population expansion
##         low, high = 0.5, 0.55
##         mutage, grate, N0, sco = 840, 0.01, 1000000, 0.0
##         path = FreqTrajectoryForward(low, high, mutage, grate, N0, sco)
## 
##     def testTrajectoryStoch(self):
##         'Testing the trajectory obtained from backward binomial sampling'
##         # fitness
##         #     AA         Aa            aa
##         #        1         1+s1        1+s2
##         # constant population size
##         # s is default to neutral process
##         path = FreqTrajectoryStoch(freq=0.3, N=10000)
##         # advantageous allele, s2>s1>0 
##         path = FreqTrajectoryStoch(freq=0.3, N=10000,fitness=[1, 1, 1.01])
##         # overdominance, s1 > s2 > 0
##         path = FreqTrajectoryStoch(freq=0.3, N=10000,fitness=[1, 1.02, 1])
##         # with week purifying selection (additive)
##         path = FreqTrajectoryStoch(freq=0.3, N=10000,fitness=[1, 0.9999, 0.9998])
##         # population growth
##         def NtFunc(gen):
##             if gen < 5000:
##                 return [10000]
##             else:
##                 return [10000*math.exp(0.001*(gen-5000))]
##         # neutral
##         path = FreqTrajectoryStoch(curGen=10000, freq=0.3, NtFunc=NtFunc)
##         # advantageous allele, s2>s1>0 
##         path = FreqTrajectoryStoch(curGen=10000, freq=0.3, NtFunc=NtFunc,fitness=[1, 1, 1.01])
##         # overdominance, s1 > s2 > 0
##         path = FreqTrajectoryStoch(curGen=10000, freq=0.3, NtFunc=NtFunc,fitness=[1, 1.02, 1])
##         # with week purifying selection (additive)
##         path = FreqTrajectoryStoch(curGen=10000, freq=0.3, NtFunc=NtFunc,fitness=[1, 0.9999, 0.9998])
##         #
##         # changing selection pressure
##         def fitnessFunc(gen):
##             if gen < 9000:    # previously positive selection
##                 return [1, 1.01, 1.02]
##             else:                     # then under purifying selection
##                 return [1, 0.999, 0.998]
##         # neutral
##         path = FreqTrajectoryStoch(curGen=10000, freq=0.3, NtFunc=NtFunc, fitnessFunc=fitnessFunc)
##         # print path
## 
##     def testTrajectoryMultiStoch(self):
##         'Testing the trajectory obtained from backward binomial sampling'
##         #path = FreqTrajectoryMultiStoch(freq=[0.1], N=10000, 
##         # fitness=[1, 1,01, 1.02], maxMutAge=100000)
##         path = FreqTrajectoryMultiStoch(freq=[0.05, 0.1], N=10000, 
##            fitness=[1, 1.01, 1.02, 1, 1.002, 1.002],
##            maxMutAge=100000)
##         # using sFunc
##         def s(gen, freq):
##             if gen < 9000:
##                 return [1, 1.01, 1.02, 1, 1.002, 1.002]
##             else:
##                 return [1, 0.99, 0.98, 1, 0.999, 0.998]
##         path = FreqTrajectoryMultiStoch(curGen=10000, 
##             freq=[0.05, 0.1], N=10000, 
##             fitnessFunc=s, maxMutAge=10000)
##         # then , with frequency dependent?
##         #print path.numTraj(), path.maxLen(), path.traj(0), path.traj(1)
##    

    def testControlledMating(self):
        'Testing controlled mating'
        # planned trajectory
        freq = FreqTrajectoryStoch(freq=0.05, N=100)
        #print freq
        # staring from when?
        burnin = 100
        mutAge = len(freq)
        # trajectory function
        # 0 ...., 100, 101, .... 100+mutAge
        #                            x                 freq
        def freqRange(gen):
            if gen <= burnin:
                # whatever
                return [0,1]
            expected = freq[gen-1-burnin]
            return [expected, expected + 0.05]
        #
        # turn On debug
        simu = simulator( population(100, loci=[1], ploidy=2), 
            controlledMating( matingScheme=randomMating(), 
                locus=0, allele=1, freqFunc=freqRange ) 
            )
        #print "Simulator created"
        simu.evolve( 
            preOps=[
                initByValue([0])
                ],
            ops=[
                pointMutator(loci=[0], 
                    toAllele=1, 
                    inds = [0],
                    at = [burnin+1],
                    stage = PreMating),
                stat(alleleFreq=[0]),
                # pyEval(r'"%d %6.4f\n"%(gen, 1-alleleFreq[0][0])', begin=burnin)
            ], 
            gen=burnin+mutAge
        )
            

        
    def testControlledRandomMating(self):
        'Testing controlled random mating'
        # planned trajectory
        freq = FreqTrajectoryStoch(freq=0.05, N=100)
        # staring from when?
        burnin = 100
        mutAge = len(freq)
        # trajectory function
        # 0 ...., 100, 101, .... 100+mutAge
        #                            x                 freq
        def freqRange(gen):
            if gen <= burnin:
                return [0]
            else:
                return [freq[gen-1-burnin]]
        #
        # turn On debug
        simu = simulator( population(100, loci=[1], ploidy=2), 
            controlledRandomMating( locus=0, allele=1, freqFunc=freqRange ) 
            )
        #print "Simulator created"
        simu.evolve( 
            preOps=[
                initByValue([0])
                ],
            ops=[
                pointMutator(loci=[0], 
                    toAllele=1, 
                    inds = [0],
                    at = [burnin+1],
                    stage = PreMating),
                stat(alleleFreq=[0]),
                #pyEval(r'"%d %6.4f\n"%(gen, 1-alleleFreq[0][0])', begin=burnin)
            ], 
            gen=burnin+mutAge
        )
        
    def testControlledMultiRandomMating(self):
        'Testing the multi-locus version of controlled random mating'
        N = 5000
        # planned trajectory
        traj = FreqTrajectoryMultiStoch(freq=[0.05, 0.10], N=N, 
            maxMutAge=500, restartIfFail=True)        
        # staring from when?
        burnin = 100
        mutAge = max([len(x) for x in traj])
        # trajectory function
        # 0 ...., 100, 101, .... 100+mutAge
        #                            x                 traj
        endingGen = burnin + mutAge
        from simuUtil import trajFunc
        expectedFreq = trajFunc(endingGen, traj)
        #
        simu = simulator( population(N, loci=[1,1], ploidy=2), 
            controlledRandomMating( loci=[0,1], 
                alleles=[1]*2, freqFunc=expectedFreq ) 
            )
        #print "Simulator created"
        simu.evolve( 
            preOps=[
                initByValue([0]*2)
                ],
            ops=[
                pointMutator(loci=[0], 
                    toAllele=1, 
                    inds = [0],
                    at = [endingGen-len(traj[0])+1],
                    stage = PreMating),
                pointMutator(loci=[1], 
                    toAllele=1, 
                    inds = [1],
                    at = [endingGen-len(traj[1])+1],
                    stage = PreMating),
                stat(alleleFreq=[0,1]),
                #pyEval(r'"%d %6.4f %6.4f\n"%(gen, 1-alleleFreq[0][0], 1-alleleFreq[1][0])', begin=burnin)
            ], 
            gen=endingGen
        )


    def testForwardTrajectorySimulation(self):
        'Testing forward trajectory simulation'
        # fTurnOnDebug(DBG_DEVEL)
        def Nt(gen, oldSize):
            if gen == 10:
                return [10000, 20000]
            else:
                return [x+gen for x in oldSize]
        #
        def fitness(gen, frq):
            return [1, 1+gen/1000., 1+gen/500.]*3
        #
        def fitness1(gen, frq):
            return [1, 1+gen/10000., 1+gen/5000.]
        # one subpopulation
        freq = [[0.08, 0.15], [0.12, 0.24], [0.15, 0.27]]
        traj = ForwardFreqTrajectory(
            curGen = 10,
            endGen = 20,
            # three loci, one subpop
            curFreq = [0.1, 0.2, 0.3],
            freq = freq,
            N = [100],
            fitnessFunc=fitness)
        if len(traj) == 0:
            for i in range(3):
                assert traj[i] >= freq[i][0] and traj[i] <= freq[i][1] 
        # for t in traj:
        #    print ', '.join(['%.2f' % x for x in t])
        # one locus
        freq = [[0.08, 0.15]]
        traj = ForwardFreqTrajectory(
            curGen=10,
            endGen=20,
            # one locus, three subpop
            curFreq=[0.1, 0.2, 0.1],
            freq=freq,
            N = [1000]*3,
            fitnessFunc=fitness1)
        if len(traj) == 0:
            lastFrq = [x[-1] for x in traj]
            lastSize = [1000]*3
            avgFrq = sum([lastFrq[x]*lastSize[x] for x in range(3)]) / sum(lastSize)
            assert avgFrq >= freq[0][0] and avgFrq <= freq[0][1] 
        # for t in traj:
        #    print ', '.join(['%.2f' % x for x in t])
        # planned trajectory
        #TurnOnDebug(DBG_DEVEL)
        freq = [[0.05, 0.15], [0.05, 0.15], [0.05, 0.15]]
        traj = ForwardFreqTrajectory(
            curGen = 10,
            endGen = 20,
            # three loci, two subpop
            curFreq = [0.05, 0.10]*3,
            freq = freq,
            NtFunc = Nt,
            fitnessFunc = fitness)
        if len(traj) == 0:
            for loc in range(3):
                lastFrq = [traj[loc*2+sp][-1] for sp in range(2)]
                lastSize = Nt[20]
                avgFrq = sum([lastFrq[x]*lastSize[x] for x in range(2)]) / sum(lastSize)
                assert avgFrq >= freq[loc][0] and avgFrq <= freq[loc][1] 
        #for t in traj:
        #    print ', '.join(['%.2f' % x for x in t])



    def testSelfMating(self):
        'Test selfing mating scheme'
        pop = population(200, loci=[3,5])
        simu = simulator(pop, selfMating())
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        #
        simu = simulator(pop, selfMating(numOffspring=2))
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)

    
    def testPyMating(self):
        'Test pyMating mating scheme'
        ver = sys.version_info[:3]
        # only python >= 2.4 supports pymating
        if ver[0] <= 2 and ver[1] < 4:
            return
        pop = population(200, loci=[3,5])
        simu = simulator(pop, pyMating(
            randomParentChooser(),
            selfingOffspringGenerator()))
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        #
        self.assertRaises(exceptions.ValueError, pyMating,
            randomParentsChooser(),
            selfingOffspringGenerator())
        #
        simu = simulator(pop, pyMating(
            randomParentChooser(),
            cloneOffspringGenerator(numOffspring=3)))
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        #
        simu = simulator(pop, pyMating(
            randomParentsChooser(),
            mendelianOffspringGenerator(numOffspring=3)))
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        #
        def pc(pop, sp):
            while True:
                for ind in range(pop.subPopSize(sp)):
                    yield ind
        simu = simulator(pop, pyMating(
            pyParentsChooser(pc),
            cloneOffspringGenerator(numOffspring=3)))
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)


    def testCloneMating(self):
        'Testing clone mating scheme'
        TurnOnDebug(DBG_MATING)
        pop = population(size=[100, 200])
        InitByFreq(pop, [0.3, 0.7])
        simu = simulator(pop, cloneMating(numOffspring=2))
        simu.step(ops=[])
        pop1 = simu.population(0)
        self.assertEqual(pop.individual(0), pop1.individual(0))
        self.assertEqual(pop.individual(0), pop1.individual(1))
        self.assertEqual(pop.individual(1), pop1.individual(2))
        self.assertEqual(pop.individual(1), pop1.individual(3))
        self.assertEqual(pop.individual(2), pop1.individual(4))
        # ...
        self.assertEqual(simu.population(0).dvars().famSizes,
            [2]*150)
        pop = population(size=[100, 200])
        InitByFreq(pop, [.3, .7])
        simu = simulator(pop, cloneMating())
        simu.step()
        self.assertEqual(simu.population(0), pop)
        TurnOffDebug(DBG_MATING)


    def testWeightSystem(self):
        'Testing the weigting system used by heterogeneous mating scheme'
        TurnOnDebug(DBG_MATING)
        def getOffSize(N, weights):
            pop = population(N)
            ms = []
            for idx,w in enumerate(weights):
                ms.append(cloneMating(numOffspring=(idx+1)*10, weight=w))
            simu = simulator(pop, heteroMating(ms))
            simu.step()
            fs = simu.dvars(0).famSizes
            ret = [0]*len(weights)
            v0 = fs[0]
            idx = 0
            for v in fs:
                if v > v0:
                    v0 = v
                    idx += 1
                ret[idx] += v
            return ret
        self.assertEqual(getOffSize(1000, [0,0]), [500, 500])
        self.assertEqual(getOffSize(1000, [-1,0]), [1000, 0])
        self.assertEqual(getOffSize(1000, [-0.3, -0.7]), [300, 700])
        self.assertEqual(getOffSize(1000, [-0.5, 0, 0]), [500, 250, 250])
        self.assertEqual(getOffSize(1000, [-0.5, 2, 3]), [500, 200, 300])
        self.assertEqual(getOffSize(1000, [-0.5, 2, 3]), [500, 200, 300])
        self.assertEqual(getOffSize(1000, [2, -0.1, 3]), [360, 100, 540])
        self.assertEqual(getOffSize(1000, [-0.2, -0.1, 0]), [200, 100, 700])
        TurnOffDebug(DBG_MATING)
    
    def testHeteroMating(self):
        'Testing heterogeneous mating schemes'
        TurnOnDebug(DBG_MATING)
        pop = population(size=[100, 200])
        simu = simulator(pop,
            heteroMating(
                [randomMating(numOffspring=1, subPop=0),
                randomMating(numOffspring=2, subPop=1)])
        )
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        # ...
        self.assertEqual(simu.population(0).dvars().famSizes,
            [1]*100+[2]*100)
        #
        simu = simulator(pop,
            heteroMating(
                [selfMating(numOffspring=1, subPop=0),
                selfMating(numOffspring=4, subPop=1)])
        )
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        # ...
        self.assertEqual(simu.population(0).dvars().famSizes,
            [1]*100+[4]*50)
        #
        simu = simulator(pop,
            heteroMating(
                [selfMating(numOffspring=1, subPop=0),
                selfMating(numOffspring=2, subPop=0),
                selfMating(numOffspring=4, subPop=1)])
        )
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        # ...
        self.assertEqual(simu.population(0).dvars().famSizes,
            [1]*50+[2]*25+[4]*50)
        # test weight
        simu = simulator(pop,
            heteroMating(
                [selfMating(numOffspring=1, subPop=0, weight=4),
                selfMating(numOffspring=2, subPop=0, weight=1),
                selfMating(numOffspring=4, subPop=1)])
        )
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        # ...
        self.assertEqual(simu.population(0).dvars().famSizes,
            [1]*80+[2]*10+[4]*50)
        #
        # set proportional splitter
        pop.setVirtualSplitter(proportionSplitter([0.6, 0.4]), 0)
        simu = simulator(pop,
            heteroMating(
                [selfMating(numOffspring=1, subPop=0, virtualSubPop=0),
                selfMating(numOffspring=2, subPop=0, virtualSubPop=1),
                selfMating(numOffspring=4, subPop=1)])
        )
        simu.evolve(
            preOps=[initByFreq([0.3, 0.7])],
            ops=[],
            gen=10)
        # ...
        self.assertEqual(simu.population(0).dvars().famSizes,
            [1]*60+[2]*20+[4]*50)

    def testSequentialParentsChooser(self):
        'Testing sequential parent chooser'
        pop = population(size=[100, 200])
        InitByFreq(pop, [.3, .7])
        self.assertRaises(exceptions.ValueError, pyMating,
            sequentialParentChooser(),
            mendelianOffspringGenerator())
        simu = simulator(pop, pyMating(
            sequentialParentsChooser(),
            mendelianOffspringGenerator()))
        simu.step()

    def testPedigreeMating(self):
        'Testing pedigree mating'
        pop = population(size=[100, 100], loci=[2,5])
        InitByFreq(pop, [0.2, 0.8])
        simu = simulator(pop, randomMating())
        simu.evolve(
            ops = [parentsTagger(output='>>pedigree.dat',
                infoFields=[])],
            gen=10
        )
        ped = pedigree('pedigree.dat')
        #
        simu1 = simulator(pop, pedigreeMating(pedigree=ped))
        simu1.evolve(
            ops = [parentsTagger(output='>>pedigree_rep.dat',
                infoFields=[])],
            gen=10
        )
        self.assertEqual(open('pedigree.dat').read(),
            open('pedigree_rep.dat').read())
        #
        ped.markUnrelated()
        ped.removeUnrelated()
        ped.save('ped_shrink.dat')
        simu2 = simulator(pop, pedigreeMating(pedigree=ped))
        simu2.evolve(
            ops = [parentsTagger(output='>>ped_shrink_rep.dat',
                infoFields=[])],
            gen=10
        )
        self.assertEqual(open('ped_shrink.dat').read(),
            open('ped_shrink_rep.dat').read())
        for file in ['pedigree.dat', 'pedigree_rep.dat',
                'ped_shrink.dat', 'ped_shrink_rep.dat']:
            os.remove(file)
        #

    def testOneParentPedigreeMating(self):
        'Testing pedigree mating in the one parent cases'
        # testing haploid case
        pop = population(size=[100, 100], loci=[2,5])
        InitByFreq(pop, [0.2, 0.8])
        simu = simulator(pop, selfMating())
        simu.evolve(
            ops = [parentTagger(output='>>pedigree.dat',
                infoFields=[])],
            gen=10
        )
        ped = pedigree('pedigree.dat')
        #
        simu1 = simulator(pop, 
            pedigreeMating(generator=selfingOffspringGenerator(), 
                pedigree=ped))
        simu1.evolve(
            ops = [parentTagger(output='>>pedigree_rep.dat',
                infoFields=[])],
            gen=10
        )
        self.assertEqual(open('pedigree.dat').read(),
            open('pedigree_rep.dat').read())
        #
        ped.markUnrelated()
        ped.removeUnrelated()
        ped.save('ped_shrink.dat')
        simu2 = simulator(pop, 
            pedigreeMating(generator=selfingOffspringGenerator(), 
            pedigree=ped))
        simu2.evolve(
            ops = [parentTagger(output='>>ped_shrink_rep.dat',
                infoFields=[])],
            gen=10
        )
        self.assertEqual(open('ped_shrink.dat').read(),
            open('ped_shrink_rep.dat').read())
        for file in ['pedigree.dat', 'pedigree_rep.dat',
                'ped_shrink.dat', 'ped_shrink_rep.dat']:
            os.remove(file)
        

    def testMateRandomSex(self):
        '''Testing assigning offspring sex by random'''
        pop = population(10000)
        simu = simulator(pop, randomMating())
        simu.evolve(
            preOps = [initByFreq([0.5, 0.5])],
            ops = [
                stat(numOfMale=True),
                # number of male should be variable, but not too much
                terminateIf('numOfMale < 4500 or numOfMale > 6500'),
            ],
            end = 10
        )
        self.assertEqual(simu.gen(), 11)


    def testMateProbOfMale(self):
        '''Testing assigning offspring sex by probability'''
        pop = population(10000)
        simu = simulator(pop, randomMating(sexParam=0.3, sexMode=MATE_ProbOfMale))
        simu.evolve(
            preOps = [initByFreq([0.5, 0.5])],
            ops = [
                stat(numOfMale=True),
                # number of male should be variable, but not too much
                terminateIf('numOfMale < 2500 or numOfMale > 3500'),
            ],
            end = 10
        )
        self.assertEqual(simu.gen(), 11)


    def testMateNumOfMale(self):
        '''Testing assigning certain number of males among offspring'''
        pop = population(10000)
        simu = simulator(pop, randomMating(numOffspring=4, sexParam=1,
            sexMode=MATE_NumOfMale))
        simu.evolve(
            preOps = [initByFreq([0.5, 0.5])],
            ops = [
                stat(numOfMale=True),
                # number of male should be variable, but not too much
                terminateIf('numOfMale != 2500'),
            ],
            end = 10
        )
        self.assertEqual(simu.gen(), 11)


    def testMateNumOfFemale(self):
        '''Testing assigning certain number of females among offspring'''
        pop = population(10000)
        simu = simulator(pop, randomMating(numOffspring=4, sexParam=1,
            sexMode=MATE_NumOfFemale))
        simu.evolve(
            preOps = [initByFreq([0.5, 0.5])],
            ops = [
                stat(numOfMale=True),
                # number of male should be variable, but not too much
                terminateIf('numOfMale != 7500'),
            ],
            end = 10
        )
        self.assertEqual(simu.gen(), 11)


##   def testFreqTrajWithSubPop(self):
##     'Testing trajctory simulation with subpopulation structure'
##     from simuUtil import FreqTrajectoryMultiStochWithSubPop
##     initSize = 10000
##     endingSize = 200000
##     burninGen = 4000
##     splitGen = 6000
##     mixingGen = 9000
##     endingGen = 10000
##     numSubPop = 3
##     numLoci = 5
##     def popSizeFunc(gen, curSize=[]):
##       if gen < burninGen:
##         return [initSize]
##       rate =  (math.log(endingSize)-math.log(initSize))/(endingGen-burninGen)
##       if gen < splitGen:
##         return [int(initSize*math.exp((gen-burninGen)*rate))]
##       else:
##         return [int(initSize*math.exp((gen-burninGen)*rate)/numSubPop)]*numSubPop
##     (traj, gens, trajFunc) = FreqTrajectoryMultiStochWithSubPop(
##       curGen = endingGen,
##       # five dsl, five subpopulation
##       numLoci=numLoci,
##       freq=[0.5]*(numSubPop*5), 
##       NtFunc=popSizeFunc, 
##       fitness=[1, 1.0007, 1.0014]*5, 
##       minMutAge=endingGen-splitGen, 
##       maxMutAge=endingGen-burninGen, 
##       restartIfFail=True)
##     self.assertEqual( len(trajFunc(splitGen)), numLoci*numSubPop)
##     self.assertEqual( len(trajFunc(splitGen-1)), numLoci)
##     for i in range(len(gens)):
##       assert trajFunc(gens[i])[i] > 0
##       assert trajFunc(gens[i]-1)[i] == 0
  
    def testHaplodiploid(self):
        'Testing recombination in haplodiploid populations'
        pop = population(size=[20, 20], ploidy=Haplodiploid, loci=[3,5])
        simu = simulator(pop, haplodiploidMating())
        simu.evolve(
            preOps = [
                initByValue([0]*8 + [1]*8)
                ],
            ops = [recombinator(rate=0.3)],
            gen = 1)
        # all individuals get the second copy from the first copy of male parents
        # which are all zero
        for ind in simu.population(0).individuals():
            self.assertEqual(ind.arrGenotype(1), [0]*8)


    def testAlphaMating(self):
        'Testing alpha mating schemes'
        # use a random alpha parent
        pop = population(size=[200, 200], loci=[3,5])
        simu = simulator(pop, alphaMating(alphaSex=Male, alphaNum=1))
        simu.evolve(
            preOps = [
                initByFreq([0.2, 0.3, 0.5])
                ],
            ops = [],
            gen = 1)
        # there is only one Male...
        ch1 = []
        ch2 = []
        for ind in simu.population(0).individuals(0):
            genotype = ind.arrGenotype(1, 0)
            if genotype not in ch1:
                ch1.append(genotype)
            genotype = ind.arrGenotype(1, 1)
            if genotype not in ch2:
                ch2.append(genotype)
        assert len(ch1) <= 2
        assert len(ch2) <= 2
        # another subpopulation
        ch1 = []
        ch2 = []
        for ind in simu.population(0).individuals(1):
            genotype = ind.arrGenotype(1, 0)
            if genotype not in ch1:
                ch1.append(genotype)
            genotype = ind.arrGenotype(1, 1)
            if genotype not in ch2:
                ch2.append(genotype)
        assert len(ch1) <= 2
        assert len(ch2) <= 2
        #
        # use a information field
        pop = population(size=[200], loci=[3,5], infoFields=['alpha'])
        InitByFreq(pop, [0.2, 0.3, 0.5])
        pop.individual(0).setInfo(1, 'alpha')
        pop.individual(1).setInfo(1, 'alpha')
        pop.individual(0).setSex(Male)
        pop.individual(1).setSex(Male)
        simu = simulator(pop, alphaMating(alphaSex=Male, alphaField='alpha'))
        simu.evolve(
            preOps = [
                #infoEval('alpha'),
                ],
            ops = [],
            gen = 1)
        # there is only one Male...
        ch1 = []
        ch2 = []
        for ind in simu.population(0).individuals(0):
            genotype = ind.arrGenotype(1, 0)
            if genotype not in ch1:
                ch1.append(genotype)
            genotype = ind.arrGenotype(1, 1)
            if genotype not in ch2:
                ch2.append(genotype)
        assert len(ch1) <= 4
        assert len(ch2) <= 4


    def testMonoMating(self):
        'Testing monogemous mating scheme'
        pop = population(size=[2000], loci=[3,5], infoFields=['father_idx', 'mother_idx'])
        InitByFreq(pop, [0.2, 0.3, 0.5])
        # exactly 100 males and 100 females
        for i in range(1000):
            pop.individual(i).setSex(Male)
            pop.individual(1000+i).setSex(Female)
        simu = simulator(pop, monogamousMating(numOffspring=2, replenish=False))
        simu.evolve(
            preOps = [],
            ops = [parentsTagger()],
            gen = 1)
        # there is only one Male...
        self.assertEqual(len(sets.Set(simu.population(0).indInfo('father_idx'))), 1000)
        self.assertEqual(len(sets.Set(simu.population(0).indInfo('mother_idx'))), 1000)
        #
        # the next generation is unlikely to have exact number of male and female so replenish is needed
        self.assertRaises(exceptions.IndexError, simu.step)
        #


    def testPolyMating(self):
        'Testing polygemous mating scheme'
        pop = population(size=[200], loci=[3,5], infoFields=['father_idx', 'mother_idx'])
        InitByFreq(pop, [0.2, 0.3, 0.5])
        # exactly 100 males and 100 females
        for i in range(100):
            pop.individual(i).setSex(Male)
            pop.individual(100+i).setSex(Female)
        simu = simulator(pop, polygamousMating(polySex=Male, polyNum=3, numOffspring=2, replacement=True))
        simu.evolve(
            preOps = [],
            ops = [parentsTagger()],
            gen = 1)
        # there is only one Male...
        fi = simu.population(0).indInfo('father_idx')
        self.assertEqual(fi[0], fi[1])
        self.assertEqual(fi[0], fi[5])
        self.assertNotEqual(fi[0], fi[6])
        mi = simu.population(0).indInfo('mother_idx')
        self.assertEqual(mi[0], mi[1])
        self.assertNotEqual(mi[0], mi[2])
        #
        #

        
if __name__ == '__main__':
  unittest.main()
  sys.exit(0)


