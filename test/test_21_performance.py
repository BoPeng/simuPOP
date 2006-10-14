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
#simuOpt.setOptions(quiet=True)

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
        # From ind: 0.180000
        # From pop, no rearrange: 0.190000
        # From pop, with rearrange: 0.270000
        # N=100000
        # From ind: 2.400000
        # From pop, no rearrange: 2.480000
        # From pop, with rearrange: 3.410000
        # N=1000000
        # From ind: 25.010000
        # From pop, no rearrange: 25.980000
        # From pop, with rearrange: 35.340000
        # 
        # This discourage the use of gappedIteartor
        # and enougrage the use of individual iterator
        # I am not sure why, but mating schemes have been 
        # rewitten according to this. 

    def TestRandomMating(self):
        'Test the performance of random mating '
        sel = maSelector(loci=[0], fitness=[1, 1-0.001/2, 1-0.001], wildtype=[0])
        r = 0.001
        migr = migrator(rate=[[1-r,r],[r,1-r]])
        p = 0.4
        for N in [10000, 100000, 1000000]:
            print "N=%d" % N
            pop = population(subPop=[N/2]*2, loci=[1], infoFields=['a', 'fitness'])
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
            TurnOnDebug(DBG_DEVEL)
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
            # a second test reveals less existing results
            #   4.24, 6.1, 7.74
            #   43.9, 65.63, 85.28
            #   ...
            # (kind of strange, system load difference?)
            # 
            # Tried to change maleIndex to IndIterator, but
            # the performance is slightly worse. (42,34, 64.35 etc)
            #
            # Change maSelector and get a little performance gain of
            #    4.26, 6.02, 7.65
            #    43.66, 64.6, 83.77

    def TestMatingAlgorithm(self):
        'Test the performance of mating algorithm'
        sel = maSelector(loci=[0], fitness=[1, 1-0.001/2, 1-0.001], wildtype=[0])
        for N in [40000, 400000]:
            print "N=%d" % N
            pop = population(subPop=[N/2]*2, loci=[1], infoFields=['a', 'fitness'])
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
            pop = population(subPop=[N/2]*2, loci=[1], infoFields=['a', 'fitness'])
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


        
if __name__ == '__main__':
    unittest.main()
