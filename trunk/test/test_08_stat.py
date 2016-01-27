#!/usr/bin/env python
#
# testing statistics calculation
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import math
import unittest, os, sys
from simuOpt import setOptions
from random import randint

setOptions(quiet=True) 
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
from simuPOP.utils import export
from simuPOP.gsl import gsl_cdf_chisq_Pinv

try:
    import rpy
    rpy.set_default_mode(rpy.DEFAULT_CONVERSION)
    has_rpy = True
except:
    has_rpy = False

class TestStat(unittest.TestCase):

    def testPopSize(self):
        'Testing calculation of population (subpopulation) size'
        pop = Population(size=[200,800])
        # do not calculate for subpopulations
        stat(pop, popSize=1, subPops=[])
        self.assertEqual(pop.dvars().subPopSize, [])
        self.assertEqual(pop.dvars().popSize, 0)
        self.assertRaises(ValueError, pop.dvars, 0)
        stat(pop, popSize=1, subPops=1, vars='popSize_sp')
        self.assertRaises(ValueError, pop.dvars, 0)
        self.assertEqual(pop.dvars(1).popSize, 800)
        # calculate for all subpopulations, using virtual subpopulation
        pop.setVirtualSplitter(SexSplitter())
        initSex(pop, sex=[MALE, FEMALE])
        stat(pop, popSize=1, subPops=[(0,0), (1,1), 1], vars=['subPopSize', 'popSize', 'popSize_sp'])
        self.assertEqual(pop.dvars().subPopSize, [100, 400, 800])
        self.assertEqual(pop.dvars([0,0]).popSize, 100)
        self.assertEqual(pop.dvars([1,1]).popSize, 400)
        # test the vars parameter
        pop = Population(size=[200,800])
        stat(pop, popSize=1, vars='popSize')
        self.assertEqual('subPopSize' in pop.vars(), False)
        self.assertEqual(pop.dvars().popSize, 1000)


    def testNumOfMales(self):
        'Testing counting number of male'
        pop = Population(size=[200, 800])
        for i in range(100):
            pop.individual(i,0).setSex(MALE)
            pop.individual(i,1).setSex(MALE)
        for i in range(100,200):
            pop.individual(i,0).setSex(FEMALE)
        for i in range(100,800):
            pop.individual(i,1).setSex(FEMALE)
        stat(pop, numOfMales=True, vars=['numOfMales', 'numOfFemales'])
        self.assertEqual(pop.dvars().numOfMales, 200)
        self.assertEqual(pop.dvars().numOfFemales, 800)
        self.assertRaises(ValueError, pop.dvars, 0)
        # all subpopulations
        stat(pop, numOfMales=True, vars=['numOfMales_sp', 'numOfFemales_sp', 'propOfMales_sp', 'propOfFemales_sp'])
        self.assertEqual(pop.dvars(0).numOfMales, 100)
        self.assertEqual(pop.dvars(0).numOfFemales, 100)
        self.assertEqual(pop.dvars(1).numOfMales, 100)
        self.assertEqual(pop.dvars(1).numOfFemales, 700)
        self.assertEqual(pop.dvars(1).propOfMales, 1./8)
        self.assertEqual(pop.dvars(1).propOfFemales, 7./8)
        # test virtual subpopulations
        pop.setVirtualSplitter(ProportionSplitter([0.4, 0.6]))
        stat(pop, numOfMales=True, subPops=[(0, 0), (1, 1)], vars=['numOfMales_sp', 'numOfFemales_sp', 'propOfFemales_sp'])
        self.assertRaises(ValueError, pop.dvars, (0, 1))
        self.assertEqual(pop.dvars([0, 0]).numOfMales, 80)
        self.assertEqual(pop.dvars([0, 0]).propOfFemales, 0)
        self.assertEqual(pop.dvars([1, 1]).numOfFemales, 480)

    def testNumOfAffected(self):
        'Testing counting number of affected individuals'
        pop = Population(size=[200, 800])
        initSex(pop, sex=[MALE, FEMALE])
        for i in range(100):
            pop.individual(i,0).setAffected(True)
            pop.individual(i,1).setAffected(True)
        for i in range(100,200):
            pop.individual(i,0).setAffected(False)
        for i in range(100,800):
            pop.individual(i,1).setAffected(False)
        stat(pop, numOfAffected=1, vars=['propOfAffected', 'propOfUnaffected', 'propOfAffected_sp', 'propOfUnaffected_sp',
            'numOfAffected', 'numOfUnaffected', 'numOfAffected_sp', 'numOfUnaffected_sp'])
        self.assertEqual(pop.dvars().numOfAffected, 200)
        self.assertEqual(pop.dvars().numOfUnaffected, 800)
        self.assertEqual(pop.dvars(0).numOfAffected, 100)
        self.assertEqual(pop.dvars(0).numOfUnaffected, 100)
        self.assertEqual(pop.dvars(1).numOfAffected, 100)
        self.assertEqual(pop.dvars(1).numOfUnaffected, 700)
        self.assertEqual(pop.dvars().propOfAffected, 0.2)
        self.assertEqual(pop.dvars().propOfUnaffected, 0.8)
        self.assertEqual(pop.dvars(0).propOfAffected, 0.5)
        self.assertEqual(pop.dvars(0).propOfUnaffected, 0.5)
        self.assertEqual(pop.dvars(1).propOfAffected, 1/8.)
        self.assertEqual(pop.dvars(1).propOfUnaffected, 7/8.)
        # virtual subpopulation?
        pop.setVirtualSplitter(SexSplitter())
        stat(pop, numOfAffected=True, subPops=[(0, 0), (1, 1)],
            vars=['numOfAffected_sp', 'propOfUnaffected_sp', 'numOfUnaffected_sp'])
        self.assertRaises(ValueError, pop.dvars, (0, 1))
        self.assertEqual(pop.dvars([0, 0]).numOfAffected, 50)
        self.assertEqual(pop.dvars([0, 0]).propOfUnaffected, 0.5)
        self.assertEqual(pop.dvars([1, 1]).numOfUnaffected, 350)

    def testNumOfSegSites(self):
        'Testing the number of segregating sites'
        pop = Population(size=1000, loci=[10]*10)
        # 100 mutations
        for i in range(100):
            pop.individual(randint(0, 999)).setAllele(i + 1, i)
        # number of segragation site should be 100
        stat(pop, numOfSegSites=ALL_AVAIL)
        self.assertEqual(pop.dvars().numOfSegSites, 100)
        #
        stat(pop, numOfSegSites=ALL_AVAIL, vars='segSites')
        self.assertEqual(len(pop.dvars().segSites), 100)
        #
        # male only?
        pop = Population(size=1000, loci=[10]*100)
        # 500 males and 500 females
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        # 100 mutations
        for i in range(100):
            # male ...
            pop.individual(i*2).setAllele(i + 1, i)
        # number of segragation site should be 100
        stat(pop, numOfSegSites=ALL_AVAIL)
        self.assertEqual(pop.dvars().numOfSegSites, 100)
        stat(pop, numOfSegSites=ALL_AVAIL, vars='segSites')
        self.assertEqual(len(pop.dvars().segSites), 100)
        stat(pop, numOfSegSites=ALL_AVAIL, subPops=[(0,0), (0,1)], vars=['numOfSegSites_sp', 'segSites_sp'])
        self.assertEqual(pop.dvars().numOfSegSites, 100)
        self.assertEqual(len(pop.dvars().segSites), 100)
        self.assertEqual(pop.dvars((0,0)).numOfSegSites, 100)
        self.assertEqual(pop.dvars((0,1)).numOfSegSites, 0)
        self.assertEqual(len(pop.dvars((0,0)).segSites), 100)
        self.assertEqual(len(pop.dvars((0,1)).segSites), 0)
        # selected loci
        stat(pop, numOfSegSites=range(100, 1000))
        self.assertEqual(pop.dvars().numOfSegSites, 0)
        stat(pop, numOfSegSites=range(100, 1000), vars='segSites')
        self.assertEqual(len(pop.dvars().segSites), 0)
        #
        # number of non-zero sites
        pop = Population(size=1000, loci=[10]*10)
        # 20 segreagating sites with 0
        initGenotype(pop, loci= range(20), freq=[0.5, 0.5])
        # 10 segragating sites without 0, fixed for binary module
        initGenotype(pop, loci= range(30, 40), freq=[0, 0.5, 0.5])
        # 10 fixed sites
        initGenotype(pop, loci= range(50, 60), freq=[0, 1])
        stat(pop, numOfSegSites=ALL_AVAIL, vars=['numOfSegSites', 'numOfFixedSites', 'segSites', 'fixedSites'])
        if moduleInfo()['alleleType'] == 'binary':
            self.assertEqual(pop.dvars().numOfSegSites, 20)
            self.assertEqual(pop.dvars().numOfFixedSites, 20)
            self.assertEqual(len(pop.dvars().segSites), 20)
            self.assertEqual(pop.dvars().segSites, list(range(20)))
            self.assertEqual(len(pop.dvars().fixedSites), 20)
            self.assertEqual(pop.dvars().fixedSites, list(range(30,40)) + list(range(50, 60)))
        else:
            self.assertEqual(pop.dvars().segSites, list(range(20)) + list(range(30,40)))
            self.assertEqual(pop.dvars().numOfSegSites, 30)
            self.assertEqual(pop.dvars().numOfFixedSites, 10)
            self.assertEqual(len(pop.dvars().segSites), 30)
            self.assertEqual(len(pop.dvars().fixedSites), 10)
            self.assertEqual(pop.dvars().fixedSites, list(range(50, 60)))
        #
        # there is a doubt if both ploidy are counted correctly
        pop = Population(size=10, loci=10)
        for idx in (1, 3, 5):
            pop.individual(idx).setAllele(1, idx, 0)
        for ind in pop.individuals():
            ind.setAllele(1, 7, 0)
            ind.setAllele(1, 7, 1)
        stat(pop, numOfSegSites=ALL_AVAIL, vars=['segSites', 'fixedSites'])
        self.assertEqual(pop.dvars().segSites, [1, 3, 5])
        self.assertEqual(pop.dvars().fixedSites, [7])
        #
        pop = Population(size=10, loci=10)
        for idx in (1, 3, 5):
            pop.individual(idx).setAllele(1, idx, 1)
        for ind in pop.individuals():
            ind.setAllele(1, 7, 0)
            ind.setAllele(1, 7, 1)
        stat(pop, numOfSegSites=ALL_AVAIL, vars=['segSites', 'fixedSites'])
        self.assertEqual(pop.dvars().segSites, [1, 3, 5])
        self.assertEqual(pop.dvars().fixedSites, [7])
        # there is however problem with sex chromosomes
        pop = Population(size=10, loci=[5,5], chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
        initSex(pop, sex=[MALE, FEMALE])
        for idx,ind in enumerate(pop.individuals()):
            if ind.sex() == MALE:
                # male has X, Y
                ind.setAllele(1, 2, 0)  # X
                ind.setAllele(1, 7, 1)  # Y
                if idx % 2 == 0:
                    ind.setAllele(1, 4, 0)  # X
                    ind.setAllele(1, 9, 1)  # Y
            else:
                ind.setAllele(1, 2, 0)  # X
                ind.setAllele(1, 2, 1)  # X
                if idx % 2 == 0:
                    ind.setAllele(1, 4, 0)  # X
                    ind.setAllele(1, 4, 1)  # Y
        #
        stat(pop, numOfSegSites=ALL_AVAIL, vars=['segSites', 'fixedSites'])
        self.assertEqual(pop.dvars().segSites, [4])
        self.assertEqual(pop.dvars().fixedSites, [2, 7, 9])



    def testNumOfMutants(self):
        'Testing the number of mutants'
        pop = Population(size=1000, loci=[100]*100)
        # 100 mutations
        for i in range(100):
            pop.individual(randint(0, 999)).setAllele(i+1, i*100)
        # number of Mutants should be 100
        stat(pop, numOfMutants=ALL_AVAIL)
        self.assertEqual(pop.dvars().numOfMutants, 100)
        pop = Population(size=10000, loci=[100]*5)
        for i in range(5):
            pop.individual(randint(0, 999)).setAllele(i+1, i*100)
        stat(pop, numOfMutants=[0,100,200,300,400])
        self.assertEqual(pop.dvars().numOfMutants, 5)
        #
        # if I set mutants back to zero?
        pop.setGenotype(0)
        stat(pop, numOfMutants=[0,100,200,300,400])
        self.assertEqual(pop.dvars().numOfMutants, 0)
        #
        #
        # male only?
        pop = Population(size=1000, loci=[10]*100)
        # 500 males and 500 females
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        # 100 mutations
        for i in range(100):
            # male ...
            pop.individual(i*2).setAllele(i + 1, i)
        stat(pop, numOfMutants=ALL_AVAIL)
        self.assertEqual(pop.dvars().numOfMutants, 100)
        stat(pop, numOfMutants=ALL_AVAIL, subPops=[(0,0), (0,1)], vars='numOfMutants_sp')
        self.assertEqual(pop.dvars().numOfMutants, 100)
        self.assertEqual(pop.dvars((0,0)).numOfMutants, 100)
        self.assertEqual(pop.dvars((0,1)).numOfMutants, 0)
        # selected loci
        stat(pop, numOfMutants=range(100, 1000))
        self.assertEqual(pop.dvars().numOfMutants, 0)


    def testDefDict(self):
        'Testing the default dictionary feature of statistics'
        pop = Population(size=1000, loci=[10])
        initGenotype(pop, freq=[0, 0.2, 0.8])
        stat(pop, alleleFreq=list(range(10)))
        d = pop.dvars().alleleFreq[0]
        if moduleInfo()['alleleType'] == 'binary':
            self.assertEqual(list(d.keys()), [1])
        else:
            self.assertEqual(list(d.keys()), [1, 2])
        self.assertEqual(d[0], 0)
        

    def testAlleleFreq(self):
        'Testing calculation of allele frequency and number of alleles'
        pop = Population(size=[500,100,1000], ploidy=2, loci = 1, lociNames='a')
        pop.setVirtualSplitter(RangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        for genotype, subPop in zip(
            [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
            initGenotype(pop, genotype=genotype, subPops=[subPop])
        stat(pop, alleleFreq=[0])
        self.assertEqual(pop.dvars().alleleNum[0], {0: 1230., 1:1970.})
        self.assertEqual(pop.dvars().alleleFreq[0], {0:1230./3200, 1:1970./3200})
        #
        stat(pop, alleleFreq=[0], subPops=(0, 1, 2), vars=['alleleFreq_sp', 'alleleNum_sp'])
        self.assertEqual(pop.dvars(0).alleleNum[0][0], 500)
        self.assertEqual(pop.dvars(0).alleleNum[0][1], 500)
        self.assertEqual(pop.dvars(0).alleleFreq[0][0], 0.5)
        self.assertEqual(pop.dvars(0).alleleFreq[0][1], 0.5)
        # can not use assertEqual for floating numbers in this case
        self.assertAlmostEqual(pop.dvars(1).alleleFreq[0][0], 13./20)
        self.assertAlmostEqual(pop.dvars(1).alleleFreq[0][1],  7./20)
        self.assertEqual(pop.dvars(1).alleleNum[0][0], 130)
        self.assertEqual(pop.dvars(1).alleleNum[0][1], 70)
        self.assertAlmostEqual(pop.dvars(2).alleleFreq[0][0], 0.3)
        self.assertAlmostEqual(pop.dvars(2).alleleFreq[0][1], 0.7)
        self.assertEqual(pop.dvars(2).alleleNum[0][0], 600)
        self.assertEqual(pop.dvars(2).alleleNum[0][1], 1400)
        pop.vars().clear()
        stat(pop, alleleFreq=[0], vars = ['alleleNum', 'alleleFreq_sp'])
        self.assertEqual('alleleNum' in pop.vars(), True)
        self.assertEqual('alleleFreq' in pop.vars(), False)
        self.assertEqual('alleleNum' in pop.vars(0), False)
        # Virtual subpopulation?
        stat(pop, alleleFreq=[0], vars = ['alleleNum', 'alleleFreq_sp'], subPops=[(0,0), (1,3), (1,4)])
        self.assertEqual(pop.dvars((0,0)).alleleFreq[0][0], 1)
        self.assertEqual(pop.dvars((1,3)).alleleFreq[0][0], 1)
        self.assertEqual(pop.dvars((1,4)).alleleFreq[0][0], 0.5)
        # does it accept ALL_AVAIL?
        pop.vars().clear()
        stat(pop, alleleFreq=ALL_AVAIL)
        self.assertEqual(pop.dvars().alleleNum[0], {0: 1230., 1:1970.})
        pop.vars().clear()
        stat(pop, alleleFreq='a')
        self.assertEqual(pop.dvars().alleleNum[0], {0: 1230., 1:1970.})
        # a special cast that fails alleleFreq
        pop = Population(size=[0,0,100,0], loci=5)
        initGenotype(pop, freq=[0.4, 0.6])
        stat(pop, alleleFreq=ALL_AVAIL)


    def testHeteroFreq(self):
        'Testing counting of heterozygote frequency'
        pop = Population(size=[500,100,1000],
            ploidy=2, loci = [1])
        pop.setVirtualSplitter(RangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        if moduleInfo()['alleleType'] == 'binary':
            for genotype, subPop in zip(
                    [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
                    [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
                initGenotype(pop, genotype=genotype, subPops=[subPop])
            stat(pop, heteroFreq=[0], vars=['heteroFreq', 'heteroNum',
                'heteroFreq_sp', 'heteroNum_sp', 'homoFreq', 'homoNum',
                'homoNum_sp', 'homoFreq_sp'])
            self.assertEqual(pop.dvars().heteroNum[0], 880)
            self.assertEqual(pop.dvars().heteroFreq[0], 880./1600)
            self.assertEqual(pop.dvars(0).heteroNum[0], 250)
            self.assertEqual(pop.dvars(0).heteroFreq[0], .5)
            self.assertEqual(pop.dvars(1).heteroNum[0], 30)
            self.assertEqual(pop.dvars(1).heteroFreq[0], 0.3)
            self.assertEqual(pop.dvars(2).heteroNum[0], 600)
            self.assertEqual(pop.dvars(2).heteroFreq[0], 0.6)
        else:
            for genotype, subPop in zip(
                [[1,1],[1,2],[2,3],[1,1],[3,2],[2,2],[1,2],[3,2],[2,2]],
                [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
                initGenotype(pop, genotype=genotype, subPops=[subPop])
            stat(pop, heteroFreq=[0], vars=['heteroFreq', 'heteroNum',
                'heteroFreq_sp', 'heteroNum_sp', 'homoFreq', 'homoNum',
                'homoNum_sp', 'homoFreq_sp'])
            self.assertEqual(pop.dvars().heteroNum[0], 1005)
            self.assertEqual(pop.dvars().heteroFreq[0], 1005./1600)
            self.assertEqual(pop.dvars(0).heteroNum[0], 375)
            self.assertEqual(pop.dvars(0).heteroFreq[0], 375/500.)
            self.assertEqual(pop.dvars(1).heteroNum[0], 30)
            self.assertEqual(pop.dvars(1).heteroFreq[0], 0.3)
            self.assertEqual(pop.dvars(2).heteroNum[0], 600)
            self.assertEqual(pop.dvars(2).heteroFreq[0], 0.6)
        # test coninuity of results
        pop = Population(size=1000, loci=2)
        initGenotype(pop, freq=[0.5, 0.5])
        stat(pop, heteroFreq=0)
        self.assertNotEqual(pop.dvars().heteroFreq[0], 0)
        self.assertRaises(KeyError, pop.dvars().heteroFreq.__getitem__, 1)
        stat(pop, heteroFreq=1)
        # calculating again at locus 1 should not override result at locus 0
        self.assertNotEqual(pop.dvars().heteroFreq[0], 0)
        self.assertNotEqual(pop.dvars().heteroFreq[1], 0)

    def testGenoFreq(self):
        'Testing the counting of genotype frequency'
        pop = Population(size=[500,100,1000], ploidy=2, loci = [1])
        pop.setVirtualSplitter(RangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        for genotype, subPop in zip(
            [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
            initGenotype(pop, genotype=genotype, subPops=[subPop])
        stat(pop, genoFreq=[0])
        #print pop.dvars().genoNum[0]
        self.assertEqual(pop.dvars().genoNum[0][(0, 0)], 175)
        self.assertEqual(pop.dvars().genoNum[0][(0, 1)], 880)
        self.assertEqual(pop.dvars().genoNum[0][(1, 1)], 545)
        self.assertEqual(pop.dvars().genoFreq[0][(0, 0)], 175./1600)
        self.assertEqual(pop.dvars().genoFreq[0][(0, 1)], 880./1600)
        self.assertEqual(pop.dvars().genoFreq[0][(1, 1)], 545./1600)
        self.assertRaises(ValueError, pop.dvars, 0)
        #
        stat(pop, genoFreq=[0], vars=['genoNum_sp', 'genoFreq_sp'])
        self.assertEqual(pop.dvars(0).genoNum[0][(0, 0)], 125)
        self.assertEqual(pop.dvars(0).genoNum[0][(0, 1)], 250)
        self.assertEqual(pop.dvars(0).genoNum[0][(1, 1)], 125)
        self.assertEqual(pop.dvars(0).genoFreq[0][(0, 0)], .25)
        self.assertEqual(pop.dvars(0).genoFreq[0][(0, 1)], .5)
        self.assertEqual(pop.dvars(0).genoFreq[0][(1, 1)], .25)
        #
        self.assertEqual(pop.dvars(1).genoNum[0][(0, 0)], 50)
        self.assertEqual(pop.dvars(1).genoNum[0][(0, 1)], 30)
        self.assertEqual(pop.dvars(1).genoNum[0][(1, 1)], 20)
        self.assertEqual(pop.dvars(1).genoFreq[0][(0, 0)], 0.5)
        self.assertEqual(pop.dvars(1).genoFreq[0][(0, 1)], 0.3)
        self.assertEqual(pop.dvars(1).genoFreq[0][(1, 1)], 0.2)
        # key does not exist
        self.assertEqual(pop.dvars(2).genoNum[0][(0, 1)], 600)
        self.assertEqual(pop.dvars(2).genoNum[0][(1, 1)], 400)
        # key does not exist
        self.assertEqual(pop.dvars(2).genoFreq[0][(0, 1)], 0.6)
        self.assertEqual(pop.dvars(2).genoFreq[0][(1, 1)], 0.4)

    def testInfoStat(self):
        'Testing summary statistics of information fields'
        import random
        pop = Population(size=[500, 1000, 1000], infoFields=['x', 'y', 'z'])
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        pop.setIndInfo([1], field='x', subPop=(0, 0))
        pop.setIndInfo([2], field='x', subPop=(0, 1))
        pop.setIndInfo([random.randint(4, 10) for x in range(500)],
            field='y', subPop=0)
        stat(pop, meanOfInfo='x', subPops=0)
        self.assertEqual(pop.dvars().meanOfInfo['x'], 1.5)
        stat(pop, sumOfInfo='x', subPops=[(0, 0)])
        self.assertEqual(pop.dvars().sumOfInfo['x'], 250)
        stat(pop, maxOfInfo='y', minOfInfo='y', subPops=[(0, 0), (0, 1)],
            vars=['maxOfInfo', 'minOfInfo', 'maxOfInfo_sp', 'minOfInfo_sp'])
        self.assertEqual(pop.dvars((0, 0)).maxOfInfo['y'], 10)
        self.assertEqual(pop.dvars((0, 1)).minOfInfo['y'], 4)
        self.assertEqual(pop.dvars().maxOfInfo['y'], 10)
        self.assertEqual(pop.dvars().minOfInfo['y'], 4)
        #
        # test cases with NO item, so mean, min, max etc should return None
        # 
        # Before 1.1.3, this results in seg dump
        pop = Population(size=[500, 0, 1000], infoFields=['x', 'y', 'z'])
        initSex(pop, sex=[MALE, FEMALE])
        stat(pop, meanOfInfo='x', minOfInfo='x', maxOfInfo='x', varOfInfo='x',
            subPops=1)
        self.assertEqual(pop.dvars().meanOfInfo['x'], None)
        self.assertEqual(pop.dvars().varOfInfo['x'], None)
        self.assertEqual(pop.dvars().minOfInfo['x'], None)
        self.assertEqual(pop.dvars().maxOfInfo['x'], None)


    def testFst(self):
        'Testing calculation of Fst value'
        pop = Population(size=[500,100,1000],
            ploidy=2, loci = [1])
        pop.setVirtualSplitter(RangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        for genotype, subPop in zip(
            [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
            initGenotype(pop, genotype=genotype, subPops=[subPop])
        #SaveFStat(simu.population(0), "p1.dat", maxAllele=2)
        # Fst is compaared with result from FStat.
        #
        #For locus : loc_0_
        #Allele    Capf	 Theta	Smallf	 Relat	Relatc	 Sig_a	 Sig_b	 Sig_w
        #        1	-0.103	 0.102	-0.229	 0.228
        #        2	-0.103	 0.102	-0.229	 0.228
        #    All	-0.103	 0.102	-0.229	 0.228	 0.373	 0.051	-0.102	 0.550
        stat(pop, structure=[0])
        self.assertTrue('F_st' in pop.vars())
        self.assertTrue('f_it' not in pop.vars())
        self.assertTrue('f_st' not in pop.vars())
        stat(pop, structure=[0], vars=['f_is', 'f_it', 'f_st', 'F_is', 'F_it'])
        self.assertAlmostEqual(pop.dvars().f_is[0], -0.2290202)
        self.assertAlmostEqual(pop.dvars().f_it[0], -0.1034594)
        self.assertAlmostEqual(pop.dvars().f_st[0],  0.1021633)
        self.assertAlmostEqual(pop.dvars().F_is, -0.2290202)
        self.assertAlmostEqual(pop.dvars().F_it, -0.1034594)
        self.assertAlmostEqual(pop.dvars().F_st,  0.1021633)
        #
        # test Fst with a subset of loci, and with unordered loci etc
        pop = Population(size=[500,100,1000],
            ploidy=2, loci = [10])
        pop.setVirtualSplitter(RangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        for genotype, subPop in zip(
            [[0,0,1],[0,1,1],[1,1,0],[0,1,0],[0,0,1],[1,1,1],[0,1,1],[0,1,1],[1,0, 1]],
            [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
            initGenotype(pop, genotype=genotype, subPops=[subPop])
        stat(pop, structure=ALL_AVAIL)
        self.assertAlmostEqual(pop.dvars().F_st,  0.0258261)
        stat(pop, structure=[2])
        self.assertAlmostEqual(pop.dvars().F_st,  0.0263055)
        stat(pop, structure=[2, 3, 5])
        self.assertAlmostEqual(pop.dvars().F_st,  0.0261665)
        # change order, but the result should be the same
        stat(pop, structure=[5, 3, 2])
        self.assertAlmostEqual(pop.dvars().F_st,  0.0261665)
        # or remove some loci
        pop.removeLoci(keep=[5,2,3])
        stat(pop, structure=ALL_AVAIL)
        self.assertAlmostEqual(pop.dvars().F_st,  0.0261665)

    def testHaploFreq(self):
        'Testing calculation of haplotype frequency'
        # test haplotype frequency
        pop = Population(size=[5000,1000], ploidy=2, loci = [10])
        if moduleInfo()['alleleType'] == 'binary':
            pop.setVirtualSplitter(ProportionSplitter([.3, .7]))
            initGenotype(pop, genotype=[0]*10, subPops=[(0,0), (1,0)])
            initGenotype(pop, genotype=[1]*10, subPops=[(0,1), (1,1)])
            stat(pop, haploFreq=[[0,1,5],[2,5]])
            self.assertEqual(pop.dvars().haploFreq[(0, 1, 5)][(0, 0, 0)], 0.3)
            self.assertEqual(pop.dvars().haploFreq[(0, 1, 5)][(1, 1, 1)], 0.7)
            self.assertEqual(pop.dvars().haploFreq[(2, 5)][(0, 0)], 0.3)
            self.assertEqual(pop.dvars().haploFreq[(2, 5)][(1, 1)], 0.7)
        else:
            pop.setVirtualSplitter(ProportionSplitter([.2, .3, .5]))
            initGenotype(pop, genotype=[1]*10, subPops=[(1,0), (0,0)])
            initGenotype(pop, genotype=[2]*10, subPops=[(1,1), (0,1)])
            initGenotype(pop, genotype=[3]*10, subPops=[(1,2), (0,2)])
            stat(pop, haploFreq=[[0,1,5], [2,5]])
            self.assertEqual(pop.dvars().haploFreq[(0, 1, 5)][(1, 1, 1)], 0.2)
            self.assertEqual(pop.dvars().haploFreq[(0, 1, 5)][(2, 2, 2)], 0.3)
            self.assertEqual(pop.dvars().haploFreq[(0, 1, 5)][(3, 3, 3)], 0.5)
            self.assertEqual(pop.dvars().haploFreq[(2, 5)][(1, 1)], 0.2)
            self.assertEqual(pop.dvars().haploFreq[(2, 5)][(2, 2)], 0.3)
            self.assertEqual(pop.dvars().haploFreq[(2, 5)][(3, 3)], 0.5) 


    def testHaploHomoFreq(self):
        'Testing calculation of haplotype homozygote frequency'
        # test haplotype frequency
        pop = Population(size=[5000,1000], ploidy=2, loci = [10])
        initSex(pop)
        initGenotype(pop, freq=[0.2, 0.3, 0.5])
        stat(pop, haploHomoFreq=[(0,1,4), (2,5)])
        for hap in [(0,1,4), (2,5)]:
            homo = sum([ [ind.allele(x, 0) for x in hap] == [ind.allele(x, 1) for x in hap] for ind in pop.individuals()])
            self.assertEqual(pop.dvars().haploHomoNum[hap], homo)
            self.assertEqual(pop.dvars().haploHomoFreq[hap], homo / 6000.)
        stat(pop, haploHeteroFreq=[(0,1,4), (2,5)])
        for hap in [(0,1,4), (2,5)]:
            homo = sum([[ind.allele(x, 0) for x in hap] == [ind.allele(x, 1) for x in hap] for ind in pop.individuals()])
            self.assertEqual(pop.dvars().haploHeteroNum[hap], 6000 - homo)
            self.assertAlmostEqual(pop.dvars().haploHeteroFreq[hap], 1 - homo / 6000., 6)
        pop.setVirtualSplitter(SexSplitter())
        stat(pop, haploHomoFreq=[(0,1,4), (2,5)], subPops=[(0,0)])
        for hap in [(0,1,4), (2,5)]:
            homo = sum([[ind.allele(x, 0) for x in hap] == [ind.allele(x, 1) for x in hap] for ind in pop.individuals([0,0])])
            self.assertEqual(pop.dvars().haploHomoNum[hap], homo)
            self.assertEqual(pop.dvars().haploHomoFreq[hap], float(homo) / pop.subPopSize([0,0]))
        # subpop ...
        stat(pop, haploHomoFreq=[(0,1,4), (2,5)], subPops=[(0,0)], vars=['haploHomoNum_sp'])
        for hap in [(0,1,4), (2,5)]:
            homo = sum([[ind.allele(x, 0) for x in hap] == [ind.allele(x, 1) for x in hap] for ind in pop.individuals([0,0])])
            self.assertEqual(pop.dvars([0,0]).haploHomoNum[hap], homo)


    def TestLD(self, freq):
        'Testing calculation of LD for a particular freq'
        #turnOnDebug('DBG_STATOR')
        #from simuUtil import viewVars
        #viewVars(pop.vars(), gui=False)
        def LD_single(var, loc1, loc2, allele1, allele2):
            p = var.alleleFreq[loc1][allele1]
            q = var.alleleFreq[loc2][allele2]
            pq = var.haploFreq[(loc1, loc2)][(allele1, allele2)]
            return pq-p*q
        def LD(var, loc1, loc2):
            LD = 0
            #allele1 is alleles in loc1
            for allele1 in range(len(var.alleleFreq[loc1])):
                for allele2 in range(len(var.alleleFreq[loc2])):
                    p = var.alleleFreq[loc1][allele1]
                    q = var.alleleFreq[loc2][allele2]
                    LD += p*q*abs(LD_single(var, loc1, loc2, allele1, allele2))
            return LD
        def LD_prime_single(var, loc1, loc2, allele1, allele2):
            p = var.alleleFreq[loc1][allele1]
            q = var.alleleFreq[loc2][allele2]
            pq = var.haploFreq[(loc1, loc2)][(allele1, allele2)]
            D = pq - p*q
            if D < 0:
                Dmax = min(p*q, (1 - p)*(1 - q))
            else:
                Dmax = min((1 - p)*q, p*(1 - q))
            return D/Dmax
        def LD_prime(var, loc1, loc2):
            LD_prime = 0
            #allele1 is alleles in loc1
            for allele1 in range(len(var.alleleFreq[loc1])):
                for allele2 in range(len(var.alleleFreq[loc2])):
                    p = var.alleleFreq[loc1][allele1]
                    q = var.alleleFreq[loc2][allele2]
                    LD_prime += p*q*abs(LD_prime_single(var, loc1, loc2, allele1, allele2))
            return LD_prime
        def R2_single(var, loc1, loc2, allele1, allele2):
            p = var.alleleFreq[loc1][allele1]
            q = var.alleleFreq[loc2][allele2]
            pq = var.haploFreq[(loc1, loc2)][(allele1, allele2)]
            return (pq-p*q)**2/(p*q*(1-p)*(1-q))
        def R2(var, loc1, loc2):
            R2 = 0
            #allele1 is alleles in loc1
            for allele1 in range(len(var.alleleFreq[loc1])):
                for allele2 in range(len(var.alleleFreq[loc2])):
                    p = var.alleleFreq[loc1][allele1]
                    q = var.alleleFreq[loc2][allele2]
                    R2 += p*q*R2_single(var, loc1, loc2, allele1, allele2)
            return R2
        pop = Population(size=[500, 100, 1000], ploidy=2, loci = [5])
        initGenotype(pop, freq=freq)
        # test case with primary alleles specified
        stat(pop, alleleFreq=[2,4], haploFreq=[2,4], LD=[2, 4, 0, 1])
        self.assertAlmostEqual(LD_single(pop.dvars(), 2, 4, 0, 1), pop.dvars().LD[2][4])
        self.assertAlmostEqual(LD_prime_single(pop.dvars(), 2, 4, 0, 1), pop.dvars().LD_prime[2][4])
        self.assertAlmostEqual(R2_single(pop.dvars(), 2, 4, 0, 1), pop.dvars().R2[2][4])
        #
        stat(pop, alleleFreq=[2,4], haploFreq=[2,4], LD=[2, 4])
        self.assertAlmostEqual(LD(pop.dvars(), 2, 4), pop.dvars().LD[2][4])
        self.assertAlmostEqual(LD_prime(pop.dvars(), 2, 4), pop.dvars().LD_prime[2][4])
        self.assertAlmostEqual(R2(pop.dvars(), 2, 4), pop.dvars().R2[2][4])
        stat(pop, alleleFreq=[2,4], haploFreq=[2,4], LD=[2, 4],
            vars=['alleleFreq_sp', 'haploFreq_sp', 'LD_sp', 'LD_prime_sp', 'R2_sp'])
        for sp in range(3):
            self.assertAlmostEqual(LD(pop.dvars(sp), 2, 4), pop.dvars(sp).LD[2][4])
            self.assertAlmostEqual(LD_prime(pop.dvars(sp), 2, 4), pop.dvars(sp).LD_prime[2][4])
            self.assertAlmostEqual(R2(pop.dvars(sp), 2, 4), pop.dvars(sp).R2[2][4])

    def testLD(self):
        '''Testing LD for both dialleleic and multi-allelic cases'''
        self.TestLD([.2, .8])
        self.TestLD([.2, .3, .5])
        #turnOnDebug('DBG_STATOR')
        pop = Population(size=[500,100,1000], ploidy=2, loci = [5])
        #
        # FIXME:
        #
        # In this test, we use and assume consecutive alleles, i.e. allele 0, 1 and 2
        # In simuPOP, these alleles can be discrete, i.e. something like
        #
        #    initGenotype(pop, freq=[0, 0, .2, .3, .5])
        #
        # This has not passed our test yet. (degree of freedom problem?)
        #
        initGenotype(pop, freq=[.2, .3, .5])
        stat(pop, alleleFreq=[2,4], haploFreq=[2,4], LD=[2,4], popSize=1,
            vars=['LD_ChiSq', 'CramerV', 'alleleNum', 'alleleFreq', 'haploNum', 'haploNum_sp'])
        def ChiSq(var, loc1, loc2):
            ChiSq = 0
            #allele1 is alleles in loc1
            N = var.popSize * 2.
            for allele1, p in var.alleleNum[loc1].items():
                for allele2, q in var.alleleNum[loc2].items():
                    pq = var.haploNum[(loc1, loc2)][(allele1, allele2)]
                    if p > 0 and q > 0:
                        ChiSq += (pq  - p * q / N) ** 2 / (p * q / N)
            return ChiSq
        def ChiSq_P(var, loc1, loc2):
            a = len(var.alleleFreq[loc1])
            b = len(var.alleleFreq[loc2])
            return 1 - rpy.r.pchisq(ChiSq(var, loc1, loc2), (a-1)*(b-1))
        def CramerV(var, loc1, loc2):
            r = len(var.alleleFreq[loc1])
            c = len(var.alleleFreq[loc2])
            CramerV = math.sqrt(ChiSq(var, loc1, loc2)/(2 * var.popSize * min(r - 1, c - 1)))
            return CramerV
        self.assertAlmostEqual(ChiSq(pop.dvars(), 2, 4), pop.dvars().LD_ChiSq[2][4])
        if has_rpy:
            self.assertAlmostEqual(ChiSq_P(pop.dvars(), 2, 4), pop.dvars().LD_ChiSq_P[2][4])
        self.assertAlmostEqual(CramerV(pop.dvars(), 2, 4), pop.dvars().CramerV[2][4])
        stat(pop, alleleFreq=[2,4], haploFreq=[2,4], LD=[2,4], popSize=1,
            vars=['alleleFreq_sp', 'alleleNum_sp', 'haploFreq_sp', 'LD_ChiSq_sp', 'CramerV_sp', 'popSize_sp'])
        for sp in range(3):
            self.assertAlmostEqual(ChiSq(pop.dvars(sp), 2, 4), pop.dvars(sp).LD_ChiSq[2][4])
            self.assertAlmostEqual(CramerV(pop.dvars(sp), 2, 4), pop.dvars(sp).CramerV[2][4])


    def testCombinedStats(self):
        '''Testing dependency of combined statistics'''
        pop = Population(size=[500,100,1000], ploidy=2, loci = [5])
        initGenotype(pop, freq=[.2, .3, .5])
        stat(pop, alleleFreq=[1,2], haploFreq=[[1,2], [1,3]], LD=[[1,2],[1,4]])
        self.assertTrue('alleleFreq' in pop.vars())
        self.assertTrue('LD' in pop.vars())
        self.assertTrue('haploFreq' in pop.vars())
        pop.dvars().LD[1][2]
        pop.dvars().LD[1][4]
        pop.dvars().alleleFreq[1][0]
        pop.dvars().alleleFreq[2][0]
        pop.dvars().haploFreq[(1, 2)]
        pop.dvars().haploFreq[(1, 3)]

    def pairwiseDiff(self, sample, loci):
        'Calculating pairwise difference'
        diff = []
        pd = sample.ploidy()
        for i in range(pd * sample.popSize() - 1):
            g1 = sample.individual(int(i / pd)).genotype(i % pd)
            for j in range(i + 1, pd * sample.popSize()):
                g2 = sample.individual(int(j / pd)).genotype(j % pd)
                diff.append(sum([x!=y for idx, x,y in zip(list(range(len(g1))), g1, g2) if idx in loci]) * 1.0)
        # return 0 if only one sequence
        if len(diff) == 0:
            return 0.0
        else:
            return sum(diff) / len(diff)

    def testNeutrality(self):
        '''Testing the calculation of Tajima's Pi'''
        pop = Population(size=1, ploidy=1, loci=[1, 1])
        initGenotype(pop, freq=[.3, .4, .3])
        pop1 = Population(size=[24, 31], ploidy=2, loci=[2, 1, 4])
        initGenotype(pop1, freq= [.3, .7])
        stat(pop, neutrality=[0, 1])
        self.assertEqual(pop.dvars().Pi, self.pairwiseDiff(pop, loci=[0, 1]))
        stat(pop1, neutrality=[0, 1, 2, 3, 4, 5, 6])
        self.assertEqual(pop1.dvars().Pi, self.pairwiseDiff(pop1, loci=[0, 1, 2, 3, 4, 5, 6]))
        stat(pop1, neutrality=[1, 4])
        self.assertEqual(pop1.dvars().Pi, self.pairwiseDiff(pop1, loci=[1, 4]))
        #stat(pop1, neutrality=[4, 4])
        #self.assertEqual(pop1.dvars().Pi, self.pairwiseDiff(pop1, loci=[4, 4]))
        stat(pop1, neutrality=[6, 1, 3])
        self.assertEqual(pop1.dvars().Pi, self.pairwiseDiff(pop1, loci=[6, 1, 3]))
        # for subpopulation
        pop1.dvars().clear()
        stat(pop1, neutrality=[1, 3, 4], vars=['Pi_sp'], suffix='_mt')
        pop1.removeSubPops(1)
        self.assertEqual(pop1.dvars(0).Pi_mt, self.pairwiseDiff(pop1, loci=[1, 3, 4]))

    def Waples89(self, S0, St, t, P0, Pt):
        # number of loci
        K_all = 0
        F_all = 0
        for loc in P0.keys():
            # number of alleles?
            alleles = set(list(P0[loc].keys()) + list(Pt[loc].keys()))
            K = len(alleles)
            if K == 1:
                continue
            # formula 9
            Fk = 0
            for allele in alleles:
                xi = P0[loc][allele]
                yi = Pt[loc][allele]
                Fk += (xi - yi) * (xi - yi) * 2. / (xi + yi)
            Fk /= K - 1
            # for multiple loci
            F_all += Fk * K
            K_all += K
        #
        F_all /= K_all
        # CI, n = Sum(K-1) is the total number of independent alleles surveyed
        n = K_all - len(P0)
        F_lower = n * F_all / gsl_cdf_chisq_Pinv(0.025, n)
        F_higher = n * F_all / gsl_cdf_chisq_Pinv(0.975, n)
        Ne = [t / (2 * (x - 0.5/S0 - 0.5/St)) for x in [F_all, F_lower, F_higher]]
        return [x if x > 0 else float('inf') for x in Ne]

    def saveToTempoFS(self, filename, S0, St, t, P0, Pt):
        with open(filename, 'w') as output:
            output.write('2 2 Test dataset\n')
            output.write('2 2\n')
            output.write('%.8f %.8f %d loc0 0\n' % (P0[0][0], P0[0][1], S0*2))
            output.write('%.8f %.8f %d loc0 %d\n' % (Pt[0][0], Pt[0][1], St*2, t))
            output.write('%.8f %.8f %d loc1 0\n' % (P0[1][0], P0[1][1], S0*2))
            output.write('%.8f %.8f %d loc1 %d\n' % (Pt[1][0], Pt[1][1], St*2, t))

    def demoNe(self, pop):
        #
        # get all lineage
        lin = [ind.alleleLineage(0, 0) for ind in pop.individuals()] + [ind.alleleLineage(0, 1) for ind in pop.individuals()]
        parent = {x: 0 for x in range(1, pop.popSize() + 1)}
        for x in lin:
            # get only the first locus
            parent[x] += 1
        # mean and variance
        N = pop.popSize()
        cnt = parent.values()
        k = sum(cnt) * 1.0 / N
        Vk = sum([(x-k)*(x-k) for x in cnt]) / (N * 1.)
        return (k * N - 1 ) / (k - 1 + Vk / k)

    def testEffectiveSize(self):
        '''Testing the effective population size estimated from genotype data'''
        setOptions(seed=1235)
        pop = Population(size=[500,100,1000], loci = 4)
        pop.setVirtualSplitter(RangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        for genotype, subPop in zip(
            [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
            initGenotype(pop, loci=0, genotype=genotype, subPops=[subPop])
            initGenotype(pop, loci=1, genotype=[genotype[1], genotype[0]], subPops=[subPop])
        stat(pop, effectiveSize=[0,1], vars='Ne_temporal_base')
        self.assertEqual(pop.dvars().Ne_temporal_base['freq'][0], {0:1230./3200, 1:1970./3200})
        self.assertEqual(pop.dvars().Ne_temporal_base['freq'][0], {0:1230./3200, 1:1970./3200})
        self.assertEqual(pop.dvars().Ne_temporal_base['size'], 1600)
        self.assertEqual(pop.dvars().Ne_temporal_base['gen'], 0)
        #
        stat(pop, effectiveSize=[0,1], subPops=(0, 1, 2), vars=['Ne_temporal_base_sp', 'Ne_waples89_sp'])
        self.assertEqual(pop.dvars(0).Ne_temporal_base['freq'][0][0], 0.5)
        self.assertEqual(pop.dvars(0).Ne_temporal_base['freq'][0][1], 0.5)
        self.assertEqual(pop.dvars(0).Ne_temporal_base['size'], 500)
        self.assertEqual(pop.dvars(1).Ne_temporal_base['size'], 100)
        self.assertEqual(pop.dvars(2).Ne_temporal_base['size'], 1000)
        self.assertEqual(pop.dvars(0).Ne_temporal_base['gen'], 0)
        self.assertEqual(pop.dvars(1).Ne_temporal_base['gen'], 0)
        self.assertEqual(pop.dvars(2).Ne_temporal_base['gen'], 0)
        #
        # export for testing using NeEstimator
        #export(pop, format='genepop', subPops=0, output='pre.txt')
        #
        t = 20
        S0 = 1600
        S00 = 500
        S01 = 100
        S02 = 1000
        P0 = pop.dvars().Ne_temporal_base['freq']
        P00 = pop.dvars(0).Ne_temporal_base['freq']
        P01 = pop.dvars(1).Ne_temporal_base['freq']
        P02 = pop.dvars(2).Ne_temporal_base['freq']
        # now, evolve population
        pop.evolve(
            initOps = [InitLineage(), InitSex()],
            matingScheme = RandomMating(subPopSize=(500, 100, 1000)),
            gen = t
        )
        # export for testing using NeEstimator
        #export(pop, format='genepop', subPops=0, output='post.txt')
        # get new allele frequency
        stat(pop, popSize=True, alleleFreq=[0,1], subPops=[0, 1, 2], vars=['popSize', 'popSize_sp', 'alleleFreq', 'alleleFreq_sp'])
        S1 = pop.dvars().popSize
        S10 = pop.dvars(0).popSize
        S11 = pop.dvars(1).popSize
        S12 = pop.dvars(2).popSize
        Pt = pop.dvars().alleleFreq
        Pt0 = pop.dvars(0).alleleFreq
        Pt1 = pop.dvars(1).alleleFreq
        Pt2 = pop.dvars(2).alleleFreq
        #
        # write to TempoFS input file for debugging
        #self.saveToTempoFS('all.txt', S0, S1, t, P0, Pt)
        #self.saveToTempoFS('sp0.txt', S00, S10, t, P00, Pt0)
        #self.saveToTempoFS('sp1.txt', S01, S11, t, P01, Pt1)
        #self.saveToTempoFS('sp2.txt', S02, S12, t, P02, Pt2)
        #
        # calculate Ne_temporal at the Python level
        #
        # for overall population and each subpopulation
        #
        #turnOnDebug('DBG_STATOR')
        stat(pop, effectiveSize=[0,1], vars=['Ne_waples89_P2', 'Ne_waples89_P2_sp', 'Ne_tempoFS_P2', 'Ne_tempoFS_P2_sp'])
        self.assertAlmostEqual(self.Waples89(S0, S1, t, P0, Pt), pop.dvars().Ne_waples89_P2)
        if moduleInfo()['threads'] == 1:
            for i in range(3):
                self.assertAlmostEqual(self.Waples89(S00, S10, t, P00, Pt0)[i], pop.dvars(0).Ne_waples89_P2[i])
                self.assertAlmostEqual(self.Waples89(S01, S11, t, P01, Pt1)[i], pop.dvars(1).Ne_waples89_P2[i])
                self.assertAlmostEqual(self.Waples89(S02, S12, t, P02, Pt2)[i], pop.dvars(2).Ne_waples89_P2[i])
            #
            self.assertAlmostEqual(pop.dvars().Ne_tempoFS_P2[0], 3026, 0)
            self.assertAlmostEqual(pop.dvars().Ne_tempoFS_P2[1], 1494, 0)
            self.assertAlmostEqual(pop.dvars(0).Ne_tempoFS_P2[0], 1620, 0)
            self.assertAlmostEqual(pop.dvars(0).Ne_tempoFS_P2[1], 543, 0)
            self.assertAlmostEqual(pop.dvars(1).Ne_tempoFS_P2[0], 94, 0)
            self.assertAlmostEqual(pop.dvars(1).Ne_tempoFS_P2[1], 28, 0)
            self.assertAlmostEqual(pop.dvars(2).Ne_tempoFS_P2[0], 2276, 0)
            self.assertAlmostEqual(pop.dvars(2).Ne_tempoFS_P2[1], 689, 0)
        #
        # 
        pop.vars().clear()
        pop.setLineage(0)
        stat(pop, effectiveSize=[0,1], vars=['Ne_demo_base'])
        if moduleInfo()['alleleType'] == 'lineage':
            l = []
            [l.extend([x,x,0,0,x,x,0,0]) for x in range(1, pop.popSize() + 1)]
            self.assertEqual(pop.lineage(), l)
            self.assertEqual(pop.dvars().Ne_demo_base, list(range(1, pop.popSize() + 1)))
        else:
            self.assertEqual(pop.dvars().Ne_demo_base, [])
        # 
        pop.evolve(matingScheme=RandomMating(), gen=1)
        # calculate Ne
        stat(pop, effectiveSize=[0,1], vars=['Ne_demo'])
        if moduleInfo()['alleleType'] == 'lineage':
            self.assertEqual(pop.dvars().Ne_demo[0], self.demoNe(pop))
            self.assertEqual(pop.dvars().Ne_demo[0], pop.dvars().Ne_demo[1])
        else:
            self.assertEqual(pop.dvars().Ne_demo[0], 0)
            self.assertEqual(pop.dvars().Ne_demo[0], pop.dvars().Ne_demo[1])
        #

    def testLDNe(self):
        # calculate LD Ne
        #turnOnDebug('DBG_STATOR')
        setOptions(seed=12345)
        pop = Population(size=[100], loci = [1]*10)
        pop.evolve(
            preOps=[
                InitSex(),
                InitGenotype(freq=[0.25]*4),
            ],
            matingScheme=RandomMating(),
            gen=10
        )
        #turnOnDebug('DBG_STATOR')
        #turnOnDebug('DBG_DEVEL')
        stat(pop, effectiveSize=ALL_AVAIL, vars=['Ne_LD', 'Ne_LD_mono'])
        if moduleInfo()['alleleType'] != 'binary' and moduleInfo()['threads'] == 1:
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.01][0], 137.6, 1)
            # parametric estimates 
            #self.assertAlmostEqual(pop.dvars().Ne_LD[1], 74.8, 1)
            #self.assertAlmostEqual(pop.dvars().Ne_LD[2], 413.4, 1)
            # Jackknife
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.01][1], 76.4, 0)
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.01][2], 383.9, 0)
            #
            # monopoly
            self.assertAlmostEqual(pop.dvars().Ne_LD_mono[0.01][0], 276.7, 1)
            # parametric estimate
            #self.assertAlmostEqual(pop.dvars().Ne_LD_mono[1], 151.1, 1)
            #self.assertAlmostEqual(pop.dvars().Ne_LD_mono[2], 828.3, 0)
            # jackknife
            self.assertAlmostEqual(pop.dvars().Ne_LD_mono[0.01][1], 154.4, 1)
            self.assertAlmostEqual(pop.dvars().Ne_LD_mono[0.01][2], 769.2, 0)
        #
        # calculate LD Ne
        #
        # calculate LD Ne
        #turnOnDebug('DBG_STATOR')
        pop = Population(size=[500], loci = [1]*10)
        pop.evolve(
            preOps=[
                InitSex(),
                InitGenotype(freq=[0.01]*5 + [0.02]*5 + [0.05]*5 + [0.2]*3),
            ],
            matingScheme=RandomMating(),
            gen=10
        )
        #export(pop, output='ldnenew.txt', format='GENEPOP', gui=False)
        #turnOnDebug('DBG_STATOR')
        #turnOnDebug('DBG_DEVEL')
        #turnOnDebug('DBG_WARNING')
        stat(pop, effectiveSize=range(10), vars=['Ne_LD', 'Ne_LD_mono'])
        #
        if moduleInfo()['alleleType'] != 'binary' and moduleInfo()['threads'] == 1:
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.05][0], 495.3, 1)
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.05][1], 377.9, 0)
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.05][2], 685.7, 0)
            #
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.02][0], 464.9, 1)
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.02][1], 397.3, 0)
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.02][2], 552.8, 0)
            #
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.01][0], 495.1, 1)
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.01][1], 423.0, 0)
            self.assertAlmostEqual(pop.dvars().Ne_LD[0.01][2], 589.5, 0)
            # monopoly
            self.assertAlmostEqual(pop.dvars().Ne_LD_mono[0.01][0], 991.7, 1)
            self.assertAlmostEqual(pop.dvars().Ne_LD_mono[0.01][1], 847.5, 1)
            self.assertAlmostEqual(pop.dvars().Ne_LD_mono[0.01][2], 1180.4, 0)
            #
        #  do not crash for fixed loci (will return nan, inf etc)
        pop = Population(size=[500], loci=[1]*10)
        stat(pop, effectiveSize=range(10), vars='Ne_LD')
 

if __name__ == '__main__':
    unittest.main()
