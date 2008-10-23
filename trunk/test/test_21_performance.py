#!/usr/bin/env python
#
# Purpose:
#    testing of performance of various parts
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 475 $
# $LastChangedDate: 2006-10-04 01:00:38 -0500 (Wed, 04 Oct 2006) $
#


#
# These tests will be time consuming so they are not automatically
# run through run_tests.py (or run_tests.sh).
#
# The test results are recorded so we will be able to compare
# performance for different revisions.
#
# Many of the functions are defined in their respective
# .cpp files, and will be available ONLY when SIMUDEBUG is defined
#
import simuOpt
simuOpt.setOptions(quiet=False)

from simuPOP import *
import unittest, os, sys, exceptions, time

class TestPerformance(unittest.TestCase):

    def TestGetinfo(self):
        '''Comparing two ways to access info fields '''
        # There are two ways to get infomation from a population
        # From individual or from population as a whole
        # Let us see which one is faster
        for N in [10000, 100000, 1000000]:
            print "N=%d" % N
            pop = population(N, loci=[1], infoFields=['a', 'b'])
            c1 = time.clock()
            for rep in range(100):
                testGetinfoFromInd(pop)
            c2 = time.clock()
            print "From ind: %f " % (c2 - c1)
            c1 = time.clock()
            for rep in range(100):
                testGetinfoFromPop(pop, False)
            c2 = time.clock()
            print "From pop, no rearrange: %f " % (c2 - c1)
            c1 = time.clock()
            for rep in range(100):
                testGetinfoFromPop(pop, True)
            c2 = time.clock()
            print "From pop, with rearrange: %f " % (c2 - c1)
        # result: (optimized module)

        # Conclusion: (strange)
        #
        # N=10000
        # From ind: 0.030000
        # From pop, no rearrange: 0.030000
        # From pop, with rearrange: 0.040000
        # N=100000
        # From ind: 0.910000
        # From pop, no rearrange: 0.820000
        # From pop, with rearrange: 1.250000
        # N=1000000
        # From ind: 10.490000
        # From pop, no rearrange: 9.710000
        # From pop, with rearrange: 14.150000
        #
        # This discourage the use of gappedIteartor
        # and enougrage the use of individual iterator
        # I am not sure why, but mating schemes have been
        # rewitten according to this.
        #
        # version 0.8.2, gcc 4.2
        #
        # N=10000
        # From ind: 0.030000
        # From pop, no rearrange: 0.030000
        # From pop, with rearrange: 0.040000
        # N=100000
        # From ind: 0.850000
        # From pop, no rearrange: 0.790000
        # From pop, with rearrange: 1.200000
        # N=1000000
        # From ind: 10.410000
        # From pop, no rearrange: 9.560000
        # From pop, with rearrange: 13.860000


    def TestRandomMating(self):
        'Testing the performance of random mating '
        sel = maSelector(loci=[0], fitness=[1, 1-0.001/2, 1-0.001], wildtype=[0])
        r = 0.001
        migr = migrator(rate=[[1-r,r],[r,1-r]])
        p = 0.4
        for N in [10000, 100000, 1000000]:
            print "N=%d" % N
            pop = population(size=[N/2]*2, loci=[1], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, randomMating())
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [],
                end=100
            )
            c2 = time.clock()
            print "From ind (no sel): %f " % (c2 - c1)
            #
            # with sel
            print "N=%d" % N
            pop = population(N, loci=[1], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, randomMating())
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [sel],
                end = 100
            )
            c2 = time.clock()
            print "From ind (sel): %f " % (c2 - c1)
            # with migr and sel
            pop = population(N, loci=[1], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, randomMating())
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [migr, sel],
                end = 100
            )
            c2 = time.clock()
            print "From ind (with migration): %f " % (c2 - c1)

            # before optimization
            # Time: (N=10^4,5,6)
            #   5.19, 6.3, 7.95
            #   54.72, 67.75, 85.74
            #   581, ,,,
            #
            # After using ind iterator
            #
            #  4.09, 5.13, 6.85
            #  42.08, 57.24, 76.16
            #  479, 621, ...
            #
            # 10/30/2006, op module
            #
            # 0.28, 0.33, 0.50
            # 3.73, 6.95, 9.07
            # 81.5, ...
            #
            # 10/30/2006, baop module
            # 0.38, 0.50, 0.69
            # 5.19, 8.73, 11.43
            # 86.81, ...
            #
            # 10/30/2006, laop module (is setSex is implemented in bitset)
            # 0.27, 0.33, 0.5 (0.29, 0.33, 0.52)
            # 3.98, 7.43, 9.8 (4.5, 7.6, 10.4)
            # 84.79, ...
            #
            #
            # Tried to change maleIndex to IndIterator, but
            # the performance is slightly worse. (42,34, 64.35 etc)
            #
            # Change maSelector and get a little performance gain of
            #    4.26, 6.02, 7.65
            #    43.66, 64.6, 83.77
            #
            #
            # Version 0.8.2, rewrite offspringGenerators, GCC 4.1.2
            #
            # 11/8/2007
            #
            # op (standard)
            # 0.31, 0.37, 0.55   (0.51, 0.66, 0.99)
            # 4.31, 7.81, 9.93   (6.15,...)
            # 93.81, 123, 153
            #
            # baop
            # 0.41, 0.51, 0.71
            # 5.6, 9.61, 12.05
            # 98.84, ...
            #
            # laop
            # 0.32, 0.35, 0.54
            # 4.43, 8.42, 10.55
            # 97.44
            #
            # Using IndIterator with m_allInds (no obvious performance penalty)
            #
            # op
            # 0.31, 0.38, 0.58
            # 4.22, 7.94, 10.42
            # 92.87, 121.83, 152.49
            #
            # baop
            # 0.41, 0.53, 0.76
            # 5.64, 9.70, 12.66
            # 100.48,
            #
            # Using IndIterator.valid(), should be quicker,
            # but there is no performance gain.
            #
            # op
            # 0.31, 0.38, 0.61
            # 4.26, 8.00, 10.8
            # 92.62, 121.88, 155.84
            #
            # baop
            # 0.42, 0.51, 0.78
            # 5.59, 9.72, 12.91
            # 97.35,
            #
            # Store IndIteartor.rawIter in parentsChooser. No improvement.
            # 0.31, 0.39, 0.64
            # 4.22, 7.70, 10.88
            # 92.92
            #
            # Using The visible algorithm of IndIterator,
            # Storing IndIterator.rawIter() in parentsChooser
            #
            # op
            # 0.32, 0.38, 0.71
            # 4.31, 7.96, 12.12
            # 95.14,
            #
            # After massive mating scheme reconstruction, they should be slower,
            # but how slow?
            #
            # op
            # 0.34, 0.41, 0.66
            # 4,67, 8.48, 11.34
            # 97.75
            #
            # laop
            #
            # baop


    def TestInfoIterator(self):
        'Testing the performance of info iterator'
        sel = maSelector(loci=[0], fitness=[1, 1-0.001/2, 1-0.001], wildtype=[0])
        r = 0.001
        p = 0.4
        for N in [10000, 100000, 1000000]:
            # with sel
            pop = population(N, loci=[1], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, randomMating())
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [sel],
                end = 100
            )
            c2 = time.clock()
            print " %.2f " % (c2 - c1),
        print
        for N in [10000, 100000, 1000000]:
            # with sel
            pop = population(N, loci=[1], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, binomialSelection())
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [sel],
                end = 100
            )
            c2 = time.clock()
            print " %.2f " % (c2 - c1),
        print
        for N in [10000, 100000, 1000000]:
            # with sel
            pop = population(N, loci=[1], infoFields=['a', 'fitness'])
            c1 = time.clock()
            for i in range(200):
                MaPenetrance(pop, locus = 0, penetrance=[0.2, 0.4, 0.8])
                info = pop.indInfo('fitness')
            c2 = time.clock()
            print " %.2f " % (c2 - c1),
        print

        #
        # With the old infoIterator,
        # op
        # 0.37   7.82   114.52
		# 0.19   4.65   58.93
		# 0.18   2.90   30.34
		#
        # laop
        # 0.38   8.26   116.15
		# 0.21   5.13   70.09
		# 0.18   2.97   30.64
		#
        # baop
        # 0.53   9.56   125.22
		# 0.33   5.99   75.34
		# 0.47   5.72   58.07
		#
		#
        # With the new individual based IndInfoIterator.
		# (2007/11/12, r1316) The performance loss is noticable.
		#
        # op:
        # 0.37   7.76  112.79
        # 0.23   4.77  61.21
        # 0.21   3.20  32.99
        #
        # laop:
        # 0.38   8.35   117.35
        # 0.25   5.37   72.33
        # 0.22   3.27   33.54
        #
        # baop:
        # 0.55   9.85   127.83
        # 0.38   6.50   80.08
        # 0.51   6.12   62.07
		#
        # With the combined iterator, we have
        #
        # op:
        # 0.38   7.94   113.23
        # 0.23   4.70   59.12
        # 0.21   3.02   31.17
        #
        # laop:
        # 0.37   8.18   114.80
        # 0.23   5.12   69.29
        # 0.20   3.00   30.77
        #
        # baop:
        # 0.54   9.67   124.93
        # 0.35   6.19   76.71
        # 0.48   5.77   58.70
        #

    def TestAlleleIterator(self):
        'Testing the performance of the new combined allele iterator'
        for N in [10000, 100000, 1000000]:
            c1 = time.clock()
            for i in range(10):
                pop = population(N, loci=[20,40])
                InitByFreq(pop, [0.2, 0.8])
                InitByFreq(pop, [0.4, 0.6])
                InitByFreq(pop, [0.5, 0.5])
                Stat(pop, alleleFreq=range(0, pop.totNumLoci()))
                Stat(pop, heteroFreq=range(0, pop.totNumLoci()))
                Stat(pop, genoFreq=range(0, pop.totNumLoci()))
            c2 = time.clock()
            print "%.2f" % (c2-c1)
        print
        #
        # Original iterator
        #
        # op:   2.20, 23.29, 231.69
        # laop: 3.04  31.96, 318.13
        # baop: 3.05 29.46 294.36
        #
        #
        # Combined iterator
        #
        # op:   2.42, 24.81, 247.89
        # laop: 3.21  32.89  326.94
        # baop: 3.78  36.89  367.50
        #
        # Use ptr directly (without gappediterator interface)
        # very strangely, there is some slight performance loss. I apply the patch any
        # way to remove the gapped iterator interface.
        #
        # op: 2.46  25.18 250.65
        # baop: 3.85 37.41 372.60


    def TestLongGenome(self):
        'Testing the performance of recombination with long genome'
        sel = maSelector(loci=[0], fitness=[1, 1-0.001/2, 1-0.001], wildtype=[0])
        r = 0.001
        migr = migrator(rate=[[1-r,r],[r,1-r]])
        p = 0.4
        for N in [10000]:
            print "N=%d" % N
            pop = population(size=[N/2]*2, loci=[1000], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, randomMating())
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [],
                end=100
            )
            c2 = time.clock()
            print "Random mating: %f " % (c2 - c1)
            #
            # with recombination
            pop = population(size=[N/2]*2, loci=[1000], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, randomMating())
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [recombinator(rate=0.0001)],
                end=100
            )
            c2 = time.clock()
            print "Low recombination: %f " % (c2 - c1)
            # with high recombination
            pop = population(size=[N/2]*2, loci=[1000], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, randomMating())
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [recombinator(rate=0.5)],
                end=100
            )
            c2 = time.clock()
            print "High recombination: %f " % (c2 - c1)
    # binary (before optimization, copy bits one by one)
    # 52.93, 76, 103
    #
    # after optimization the binaries (use block copying)
    # 2.52, 28.58, 106.9 (special treatment is given for high recombination rate cases)
    #                    (otherwise 140.0 s)
    #
    # short
    # 4.38, 33.73, 56.03
    #
    # long allele
    # 5.53, 34.81, 57.75
    #
    #
    # using the special case for low recombination rate (use the search-valid algorithm),
    #
    # binary: 2.54,  18.37, 104.5
    # short:  4.44,  21.26, 56.61
    # long:   5.52,  22.19, 57.57
    #
    # Using vector<bool> for bernulli trials, and cache actual pointers,
    #
    # binary: 2.54,  16.05, 100.48
    # short:  4.52,  18.64, 50.28
    # long:   5.59,  20.04, 52.00
    #
    # 2007/11/9, version 0.8.2
    # binary: 2.28, 15.61, 98.95
    # short: 5.32, 18.62, 50.72
    # long: 6.32, 19.82, 52.03


    def TestBernulliTrials(self):
        'Testing the performance of bernulli trials'
        rg = rng()
        p = [0.00001, 0.001, 0.5, 0.99]
        N = 1000000
        for pi in p:
            c1 = time.clock()
            bt = BernulliTrials(rg, [pi]*100, N)
            for rep in range(400):
                bt.doTrial()
            c2 = time.clock()
            print "p = %f: %f " % (pi, c2 - c1)
        #
    #
    #
    # using vector<bool> as storage engine
    #
    # p = 0.000010: 2.870000
    # p = 0.001000: 13.650000
    # p = 0.100000: 1212.650000
    # p = 0.500000: 139.860000
    # p = 0.990000: 112.390000
    #
    # Using dynamic bitset
    #
    # p = 0.000010: 2.500000
    # p = 0.001000: 14.400000
    # p = 0.500000: 142.570000
    # p = 0.990000: 111.770000

    #





    def TestMatingAlgorithm(self):
        'Testing the performance of mating algorithm'
        sel = maSelector(loci=[0], fitness=[1, 1-0.001/2, 1-0.001], wildtype=[0])
        for N in [40000, 400000]:
            print "N=%d" % N
            pop = population(size=[N/2]*2, loci=[1], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, binomialSelection())
            p = 0.4
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [],
                end = 100
            )
            c2 = time.clock()
            print "From ind: %f " % (c2 - c1)
            print "Random mating"
            pop = population(size=[N/2]*2, loci=[1], infoFields=['a', 'fitness'])
            c1 = time.clock()
            simu = simulator(pop, randomMating())
            p = 0.4
            simu.evolve(
                preOps = [initByFreq([1-p]+[p/10.]*10)],
                ops = [],
                end = 100
            )
            c2 = time.clock()
            print "From ind: %f " % (c2 - c1)

        # The version 0.7.3:
        #      0.69 (bin), 1.94 (random mating), 14.63, 34.45
        #   binary:
        #      0.98, 2.41, 17.9, 36.77
        #
        # change indIterator to reference
        # use iterator comparison
        # separate generateOffspring function
        #
        # std:
        #   0.72, 1.96, 15.05, 34.57
        # binary:
        #   0.97, 2.42, 17.22, 39.64
        #
        # STRANGE: I can not reproduce this wonderful result
        # something went wrong with this run...
        #   0.17, 1.26, 1.67, 12.68
        #
        # long:
        #   0.78, 1.38, 17.89, 53.31
        # binary:
        #   1.00, 1.79, 17.56, 34.66
        #

    def TestSimuComplexDisease(self):
        ''' recording the runtime of simuComplexDisease.py '''
        cmd = '''time python /home/bpeng/simuPOP/scripts/simuComplexDisease.py --noDialog
            --numChrom='8' --optimized
            --numLoci='10' --markerType='SNP' --DSL='(15, 26, 37, 48)'
            --DSLLoc='(0.5, 0.5, 0.5, 0.5)' --initSize='1000' --endingSize='50000'
            --growthModel='exponential' --burninGen='300' --splitGen='500' --mixingGen='800'
            --endingGen='1000' --numSubPop='1' --migrModel='stepping stone' --migrRate='0.0001'
            --alleleDistInSubPop='even' --curAlleleFreq="[0.1]*4" --minMutAge='0' --maxMutAge='0'
            --fitness='(1, 1, 1)' --selMultiLocusModel='none' --mutaRate='0.0001'
            --recRate='[0.0005]' --savePop='[]' --simuName='simu' --saveFormat='txt'
            '''
        cmd = ' '.join(cmd.split())
        print cmd
        os.system(cmd)
        #
        # 0.7.3  (baop):       the same (134u, 2min 14s)
        #
        # 0.7.4: use no-stack algorithm: 2min 29s (why more?)
        #
        # Use stack (mixed method)
        # 0.7.4 (baop):       89u, 1min 17s
        #  direct gepped allele count does oot improve anything
        #
        # 0.7.5 (fast binary operation)
        # 44.5u
        #
        # Remove rebombinator: 57.5s
        #
        # Something wrong. only 41s is needed.
        #

    def TestRecombinationAlgorithm0(self):
        ''' Testing the performance of different recombination algorithms '''
        # form 0
        c1 = time.clock()
        simu = simulator(
            population(size=[1000]*10, loci=[200]),
            randomMating(numOffspring=1/3., mode=MATE_GeometricDistribution),
            rep=1)
        simu.evolve(
            preOps = [initByValue([50]*200)],
            ops = [
               #recombinator(rate=0.5),
               smmMutator(rate=0.0001),
               pyEval(r'"%s\n"%gen', step=100),
           ],
           end = 500
        )
        c2 = time.clock()
        print "One chromosome with no recombination: %f " % (c2 - c1)

    def TestRecombinationAlgorithm1(self):
        ''' Testing the performance of different recombination algorithms '''
        # form 1
        c1 = time.clock()
        simu = simulator(
            population(size=[1000]*10, loci=[200]),
            randomMating(numOffspring=1/3., mode=MATE_GeometricDistribution),
            rep=1)
        simu.evolve(
            preOps = [initByValue([50]*200)],
            ops = [
               recombinator(rate=0.5),
               smmMutator(rate=0.0001),
               pyEval(r'"%s\n"%gen', step=100),
           ],
           end = 500
        )
        c2 = time.clock()
        print "Recombinator with rate 0.5: %f " % (c2 - c1)
        #

    def TestRecombinationAlgorithm2(self):
        # form 2
        c1 = time.clock()
        simu = simulator(
            population(size=[1000]*10, loci=[1]*200),
            randomMating(numOffspring=1/3., mode=MATE_GeometricDistribution),
            rep=1)
        simu.evolve(
            preOps = [initByValue([50]*200)],
            ops = [
               smmMutator(rate=0.0001),
               pyEval(r'"%s\n"%gen', step=100),
           ],
           end = 500
        )
        c2 = time.clock()
        print "200 chromosomes: %f " % (c2 - c1)

    def TestRecombinationAlgorithm3(self):
        # form 1
        c1 = time.clock()
        simu = simulator(
            population(size=[1000]*10, loci=[200]*10),
            randomMating(numOffspring=1/3., mode=MATE_GeometricDistribution),
            rep=1)
        simu.evolve(
            preOps = [initByValue([50]*200)],
            ops = [
               recombinator(rate=0.0005),
               smmMutator(rate=0.0001),
               pyEval(r'"%s\n"%gen', step=10),
           ],
           end = 100
        )
        c2 = time.clock()
        print "20x10 with low rec rate: %f " % (c2 - c1)
        # case 0: no rec, 1 chrom:    6.74
        # case 1: 0.5 rec, 200 chrom: 49.33
        # case 2: 200 chrom, no rec:  27.28
        # case 3: 20 chrom, low rec:  57.06
        #
        # algorith2: (rec list first):
        #
        # case 1: 48.14s
        # case 3: 59s
        #
        # get trail result directly (no BoolResults)
        #   case 0: 5.52
        #   case 1L 47.96
        #   case 2: 25.83
        #   case 3: 59.00, 55.8
        #
        # conclusion: accept this.

    def TestSerialization(self):
        import stat
        pop = LoadPopulation('exp3_0_9.txt')
        #pop = population(size=100000, loci=[10]*8, ancestralDepth=2)
        #InitByFreq(pop, [.2, .8])
        #simu = simulator(pop, randomMating(), rep=1)
        #simu.evolve(ops=[], end=2)
        #pop = simu.population(0)
        print 'Start saving file'
        for format in ['txt', 'bin', 'xml']:
            for comp in [True, False]:
                if comp:
                    label = 'small'
                else:
                    label = 'large'
                c1 = time.clock()
                pop.savePopulation('exp_%s.%s' % (label, format), compress=comp)
                c2 = time.clock()
                print "%s, save, %s: %.1f, size: %.2fM " % (format, label, c2 - c1,
                    os.stat('exp_%s.%s' % (label, format) )[stat.ST_SIZE]/1024./1024.)
                c1 = time.clock()
                pop1 = LoadPopulation('exp_%s.%s' % (label, format))
                self.assertEqual(pop, pop1)
                c2 = time.clock()
                print "%s, load, %s: %.1f " % (format, label, c2 - c1)

        #
        #  xml, save, comp: 842.8, size: 166.50M
        #  xml, load, comp: 1843.7
        #  xml, save, nocomp: 554.7, size: 8432.21M
        #  xml, load, nocomp: 1801.3
        #
        #  txt, save, comp: 225.2, size: 112.29M
        #  txt, load, comp: 128.4
        #  txt, save, nocomp: 80.4, size: 1023.27M
        #  txt, load, nocomp: 129.2
        #  bin, save, comp: 155.5, size: 99.22M
        #  bin, load, comp: 30.5
        #  bin, save, nocomp: 23.1, size: 548.83M
        #  bin, load, nocomp: 24.8

        #
        # BINARY MODULES
        # Using the new 32bit/block save/load method
        #
        # (default)
        # txt, save, small: 72.4, size: 98.98M
        # txt, load, small: 33.9
        # txt, save, large: 24.8, size: 234.12M
        # txt, load, large: 31.8
        # bin, save, small: 50.4, size: 98.20M
        # bin, load, small: 16.6
        # bin, save, large: 8.5, size: 192.96M
        # bin, load, large: 13.5
        # xml, save, small: 191.0, size: 120.16M
        
    
    def TestArrGenotype(self):
        '''Comparing the time for genotype and arrGenotype '''
        for N in [100, 1000, 10000]:
            print "N=%d" % N
            pop = population(10000, loci=[N])
            c1 = time.clock()
            for ind in pop.individuals():
				ind.setGenotype([1]*(pop.totNumLoci()*2))
            c2 = time.clock()
            print "For genotype writing: %f " % (c2 - c1)
            c1 = time.clock()
            for ind in pop.individuals():
				sum(ind.genotype())
            c2 = time.clock()
            print "For genotype reading: %f " % (c2 - c1)
            c1 = time.clock()
            for ind in pop.individuals():
				geno = ind.arrGenotype()
				geno[:] = [1]*(pop.totNumLoci()*2)
            c2 = time.clock()
            print "For arrGenotype writing: %f " % (c2 - c1)
            c1 = time.clock()
            for ind in pop.individuals():
				geno = ind.arrGenotype()
				sum(geno)
            c2 = time.clock()
            print "For arrGenotype reading: %f " % (c2 - c1)
        
        #
        # Compare the performance of returning a tuple of
        # genotype and using arr based genotype access functions.
        #
        # RESULT:
        #   arrGenotype() and genotype() are comparable in reading
        #   arrGenotype() is faster than genotype() in writing under windows
        #           but slower under linux.
        #   Because genotype() needs to allocate memory for a tuple
        #       arrGenotype is preferred in reading.
        #   We add setGenotype() because this function is easier to use.
        #
        # Windows result
        # Short modual:
		# N=100
		# For genotype writing: 0.331818
		# For genotype reading: 0.160001
		# For arrGenotype writing: 0.273920
		# For arrGenotype reading: 0.148057
		# N=1000
		# For genotype writing: 2.415126
		# For genotype reading: 1.356521
		# For arrGenotype writing: 2.303960
		# For arrGenotype reading: 1.294205
		# N=10000
		# For genotype writing: 23.452248
		# For genotype reading: 13.992072
		# For arrGenotype writing: 22.342265
		# For arrGenotype reading: 12.868873
		# Long modual:
		# N=100
		# For genotype writing: 0.321894
		# For genotype reading: 0.160523
		# For arrGenotype writing: 0.269319
		# For arrGenotype reading: 0.148007
		# N=1000
		# For genotype writing: 2.324584
		# For genotype reading: 1.376434
		# For arrGenotype writing: 2.301938
		# For arrGenotype reading: 1.295467
		# N=10000
		# For genotype writing: 23.234146
		# For genotype reading: 13.494426
		# For arrGenotype writing: 22.797773
		# For arrGenotype reading: 12.826390
		# Binary modual:
		# N=100
		# For genotype writing: 0.435293
		# For genotype reading: 0.183319
		# For arrGenotype writing: 0.274531
		# For arrGenotype reading: 0.169256
		# N=1000
		# For genotype writing: 3.648646
		# For genotype reading: 1.568618
		# For arrGenotype writing: 2.374255
		# For arrGenotype reading: 1.476079
		# N=10000
		# For genotype writing: 35.354781
		# For genotype reading: 15.545612
		# For arrGenotype writing: 23.066307
		# For arrGenotype reading: 14.525967
		
		# Linux result
		# Short modual:
		# N=100
		# For genotype writing: 0.240000 
		# For genotype reading: 0.190000 
		# For arrGenotype writing: 0.310000 
		# For arrGenotype reading: 0.160000 
		# N=1000
		# For genotype writing: 1.230000 
		# For genotype reading: 1.270000 
		# For arrGenotype writing: 2.250000 
		# For arrGenotype reading: 1.010000 
		# N=10000
		# For genotype writing: 11.050000 
		# For genotype reading: 11.950000 
		# For arrGenotype writing: 22.360000 
		# For arrGenotype reading: 9.510000
		# Long modual
		# N=100
		# For genotype writing: 0.240000 
		# For genotype reading: 0.190000 
		# For arrGenotype writing: 0.310000 
		# For arrGenotype reading: 0.160000 
		# N=1000
		# For genotype writing: 1.210000 
		# For genotype reading: 1.240000 
		# For arrGenotype writing: 2.320000 
		# For arrGenotype reading: 1.040000 
		# N=10000
		# For genotype writing: 11.300000 
		# For genotype reading: 11.700000 
		# For arrGenotype writing: 22.870000 
		# For arrGenotype reading: 9.800000
		# Binary modual:
		# N=100
		# For genotype writing: 0.230000 
		# For genotype reading: 0.230000 
		# For arrGenotype writing: 0.300000 
		# For arrGenotype reading: 0.170000 
		# N=1000
		# For genotype writing: 1.280000 
		# For genotype reading: 1.780000 
		# For arrGenotype writing: 2.130000 
		# For arrGenotype reading: 1.280000 
		# N=10000
		# For genotype writing: 12.520000 
		# For genotype reading: 17.180000 
		# For arrGenotype writing: 21.120000 
		# For arrGenotype reading: 12.240000

if __name__ == '__main__':
    unittest.main()
