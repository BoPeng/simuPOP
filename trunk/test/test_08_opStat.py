#!/usr/bin/env python
#
# testing statistics calculation
#
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision$
# $LastChangedDate$
# 


import simuOpt
simuOpt.setOptions(quiet=False)

from simuPOP import *
import unittest, os, sys, exceptions

class TestStat(unittest.TestCase):

    def testPopSize(self):
        'Testing calculation of population (subPopulation) size'
        pop = population(subPop=[200,800])
        Stat(pop, popSize=1)
        self.assertEqual(pop.dvars().numSubPop, 2)
        self.assertEqual(pop.dvars().popSize, 1000)
        self.assertEqual(pop.dvars(0).popSize, 200)
        self.assertEqual(pop.dvars(1).popSize, 800)
        self.assertEqual(pop.dvars().subPopSize, [200,800])

    def testNumOfMale(self):
        'Testing counting number of male'
        pop = population(subPop=[200, 800])
        for i in range(100):
            pop.individual(i,0).setSex(Male)
            pop.individual(i,1).setSex(Male)
        for i in range(100,200):
            pop.individual(i,0).setSex(Female)
        for i in range(100,800):
            pop.individual(i,1).setSex(Female)
        Stat(pop, numOfMale=1)
        self.assertEqual(pop.dvars().numOfMale, 200)
        self.assertEqual(pop.dvars().numOfFemale, 800)
        self.assertEqual(pop.dvars(0).numOfMale, 100)
        self.assertEqual(pop.dvars(0).numOfFemale, 100)
        self.assertEqual(pop.dvars(1).numOfMale, 100)
        self.assertEqual(pop.dvars(1).numOfFemale, 700)
            
    def testNumOfAffected(self):
        'Testing counting number of affected individuals'
        pop = population(subPop=[200, 800])
        for i in range(100):
            pop.individual(i,0).setAffected(True)
            pop.individual(i,1).setAffected(True)
        for i in range(100,200):
            pop.individual(i,0).setAffected(False)
        for i in range(100,800):
            pop.individual(i,1).setAffected(False)
        Stat(pop, numOfAffected=1)
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
        
    def testAlleleFreq(self):
        'Testing calculation of allele frequency and number of alleles'
        pop = population(subPop=[500,100,1000], 
            ploidy=2, loci = [1])
        InitByValue(pop, 
            value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            indRange = [[0,125], [125,375],[375,500],[500,550],
                [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
        Stat(pop, alleleFreq=[0], numOfAlleles=[0])
        self.assertEqual(pop.dvars().alleleFreq[0], [1230./3200, 1970./3200])
        self.assertEqual(pop.dvars().alleleNum[0], [1230, 1970])
        self.assertEqual(pop.dvars(0).alleleFreq[0], [.5, .5])
        self.assertEqual(pop.dvars(0).alleleNum[0], [500, 500])
        # can not use assertEqual for floating numbers in this case
        assert abs(pop.dvars(1).alleleFreq[0][0] - 13./20) < 1e-5
        assert abs(pop.dvars(1).alleleFreq[0][1] -  7./20) < 1e-5
        self.assertEqual(pop.dvars(1).alleleNum[0], [130, 70])
        assert abs(pop.dvars(2).alleleFreq[0][0] - 0.3) < 1e-5
        assert abs(pop.dvars(2).alleleFreq[0][1] - 0.7) < 1e-5
        self.assertEqual(pop.dvars(2).alleleNum[0], [600, 1400])
        self.assertEqual(pop.dvars().numOfAlleles[0], 2)
        self.assertEqual(pop.dvars(0).numOfAlleles[0], 2)
        self.assertEqual(pop.dvars(1).numOfAlleles[0], 2)

    def testNumOfAlleles(self):
        'Testing calculation of number of alleles'
        pop = population(subPop=[5000,15000], loci=[1])
        InitByFreq(pop, [.2, 0, .5, .3])
        Stat(pop, numOfAlleles=[0])
        if alleleType() == 'binary':
            self.assertEqual(pop.dvars().numOfAlleles[0], 2)
            self.assertEqual(pop.dvars(0).numOfAlleles[0], 2)
            self.assertEqual(pop.dvars(1).numOfAlleles[0], 2)
        else:
            self.assertEqual(pop.dvars().numOfAlleles[0], 3)
            self.assertEqual(pop.dvars(0).numOfAlleles[0], 3)
            self.assertEqual(pop.dvars(1).numOfAlleles[0], 3)
        
    def testHeteroFreq(self):
        'Testing counting of heterozygote frequency'
        pop = population(subPop=[500,100,1000], 
            ploidy=2, loci = [1])
        if alleleType() == 'binary':
            InitByValue(pop, 
                value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
                indRange = [[0,125], [125,375],[375,500],[500,550],
                    [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
            Stat(pop, heteroFreq=[0])
            self.assertEqual(pop.dvars().HeteroNum[0], 880)
            self.assertEqual(pop.dvars().heteroNum[0][0], 880)
            self.assertEqual(pop.dvars().heteroNum[0][1], 880)
            self.assertEqual(pop.dvars().HeteroFreq[0], 880./1600)
            self.assertEqual(pop.dvars().heteroFreq[0][0], 880./1600)
            self.assertEqual(pop.dvars().heteroFreq[0][1], 880./1600)
            self.assertEqual(pop.dvars(0).HeteroNum[0], 250)
            self.assertEqual(pop.dvars(0).heteroNum[0][0], 250)
            self.assertEqual(pop.dvars(0).heteroNum[0][1], 250)
            self.assertEqual(pop.dvars(0).HeteroFreq[0], .5)
            self.assertEqual(pop.dvars(0).heteroFreq[0][0], .5)
            self.assertEqual(pop.dvars(0).heteroFreq[0][1], .5)
            self.assertEqual(pop.dvars(1).HeteroNum[0], 30)
            self.assertEqual(pop.dvars(1).heteroNum[0][0], 30)
            self.assertEqual(pop.dvars(1).heteroNum[0][1], 30)
            self.assertEqual(pop.dvars(1).HeteroFreq[0], 0.3)
            self.assertEqual(pop.dvars(1).heteroFreq[0][0], 0.3)
            self.assertEqual(pop.dvars(1).heteroFreq[0][1], 0.3)
            self.assertEqual(pop.dvars(2).HeteroNum[0], 600)
            self.assertEqual(pop.dvars(2).heteroNum[0][0], 600)
            self.assertEqual(pop.dvars(2).heteroNum[0][1], 600)
            self.assertEqual(pop.dvars(2).HeteroFreq[0], 0.6)
            self.assertEqual(pop.dvars(2).heteroFreq[0][0], 0.6)
            self.assertEqual(pop.dvars(2).heteroFreq[0][1], 0.6)
        else:
            InitByValue(pop, 
                value = [[1,1],[1,2],[2,3],[1,1],[3,2],[2,2],[1,2],[3,2],[2,2]],
                indRange = [[0,125], [125,375],[375,500],[500,550],
                    [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
            Stat(pop, heteroFreq=[0])
            self.assertEqual(pop.dvars().HeteroNum[0], 1005)
            self.assertEqual(pop.dvars().heteroNum[0][1], 350)
            self.assertEqual(pop.dvars().heteroNum[0][2], 1005)
            self.assertEqual(pop.dvars().heteroNum[0][3], 655)
            self.assertEqual(pop.dvars().HeteroFreq[0], 1005./1600)
            self.assertEqual(pop.dvars().heteroFreq[0][1], 350./1600)
            self.assertEqual(pop.dvars().heteroFreq[0][2], 1005./1600)
            self.assertEqual(pop.dvars().heteroFreq[0][3], 655./1600)
            self.assertEqual(pop.dvars(0).HeteroNum[0], 375)
            self.assertEqual(pop.dvars(0).heteroNum[0][1], 250)
            self.assertEqual(pop.dvars(0).heteroNum[0][2], 375)
            self.assertEqual(pop.dvars(0).heteroNum[0][3], 125)
            self.assertEqual(pop.dvars(0).HeteroFreq[0], 375/500.)
            self.assertEqual(pop.dvars(0).heteroFreq[0][1], .5)
            self.assertEqual(pop.dvars(0).heteroFreq[0][2], 375/500.)
            self.assertEqual(pop.dvars(0).heteroFreq[0][3], 125/500.)
            self.assertEqual(pop.dvars(1).HeteroNum[0], 30)
            self.assertEqual(pop.dvars(1).heteroNum[0][1], 0)
            self.assertEqual(pop.dvars(1).heteroNum[0][2], 30)
            self.assertEqual(pop.dvars(1).heteroNum[0][3], 30)
            self.assertEqual(pop.dvars(1).HeteroFreq[0], 0.3)
            self.assertEqual(pop.dvars(1).heteroFreq[0][1], 0)
            self.assertEqual(pop.dvars(1).heteroFreq[0][2], 0.3)
            self.assertEqual(pop.dvars(1).heteroFreq[0][3], 0.3)
            self.assertEqual(pop.dvars(2).HeteroNum[0], 600)
            self.assertEqual(pop.dvars(2).heteroNum[0][1], 100)
            self.assertEqual(pop.dvars(2).heteroNum[0][2], 600)
            self.assertEqual(pop.dvars(2).heteroNum[0][3], 500)
            self.assertEqual(pop.dvars(2).HeteroFreq[0], 0.6)
            self.assertEqual(pop.dvars(2).heteroFreq[0][1], 0.1)
            self.assertEqual(pop.dvars(2).heteroFreq[0][2], 0.6)
            self.assertEqual(pop.dvars(2).heteroFreq[0][3], 0.5)
        
    def testExpHetero(self):
        'Testing expected heterozygosity 1-sum p_i2'
        pop = population(subPop=[500,100,1000], 
            ploidy=2, loci = [1])
        InitByValue(pop, 
            value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            indRange = [[0,125], [125,375],[375,500],[500,550],
                [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
        Stat(pop, expHetero=[0])
        #
        assert abs(pop.dvars().expHetero[0] - (1-(123./320)**2-(197./320)**2)) < 0.00001
        assert abs(pop.dvars(0).expHetero[0] - (1-0.5**2-0.5**2)) < 0.00001
        assert abs(pop.dvars(1).expHetero[0] - (1-(13./20)**2-(7./20)**2)) < 0.00001
        assert abs(pop.dvars(2).expHetero[0] - (1-0.3**2-0.7**2)) < 0.00001

    def testGenoFreq(self):
        'Testing the counting of genotype frequency'
        pop = population(subPop=[500,100,1000], 
            ploidy=2, loci = [1])
        InitByValue(pop, 
            value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            indRange = [[0,125], [125,375],[375,500],[500,550],
                [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
        Stat(pop, genoFreq=[0])
        self.assertEqual(pop.dvars().genoNum[0][0][0], 175)
        self.assertEqual(pop.dvars().genoNum[0][0][1], 880)
        self.assertEqual(pop.dvars().genoNum[0][1][1], 545)
        self.assertEqual(pop.dvars().genoFreq[0][0][0], 175./1600)
        self.assertEqual(pop.dvars().genoFreq[0][0][1], 880./1600)
        self.assertEqual(pop.dvars().genoFreq[0][1][1], 545./1600)
        #
        self.assertEqual(pop.dvars(0).genoNum[0][0][0], 125)
        self.assertEqual(pop.dvars(0).genoNum[0][0][1], 250)
        self.assertEqual(pop.dvars(0).genoNum[0][1][1], 125)
        self.assertEqual(pop.dvars(0).genoFreq[0][0][0], .25)
        self.assertEqual(pop.dvars(0).genoFreq[0][0][1], .5)
        self.assertEqual(pop.dvars(0).genoFreq[0][1][1], .25)
        #
        self.assertEqual(pop.dvars(1).genoNum[0][0][0], 50)
        self.assertEqual(pop.dvars(1).genoNum[0][0][1], 30)
        self.assertEqual(pop.dvars(1).genoNum[0][1][1], 20)
        self.assertEqual(pop.dvars(1).genoFreq[0][0][0], 0.5)
        self.assertEqual(pop.dvars(1).genoFreq[0][0][1], 0.3)
        self.assertEqual(pop.dvars(1).genoFreq[0][1][1], 0.2)
        # key does not exist
        self.assertEqual(pop.dvars(2).genoNum[0][0].setdefault(0,0), 0)
        self.assertEqual(pop.dvars(2).genoNum[0][0][1], 600)
        self.assertEqual(pop.dvars(2).genoNum[0][1][1], 400)
        # key does not exist
        self.assertEqual(pop.dvars(2).genoFreq[0][0].setdefault(0,0), 0)
        self.assertEqual(pop.dvars(2).genoFreq[0][0][1], 0.6)
        self.assertEqual(pop.dvars(2).genoFreq[0][1][1], 0.4)

    def testFst(self):
        'Testing calculation of Fst value'    
        pop = population(subPop=[500,100,1000], 
            ploidy=2, loci = [1])
        InitByValue(pop, 
            value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            indRange = [[0,125], [125,375],[375,500],[500,550],
                [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
        #SaveFstat(simu.population(0), "p1.dat", maxAllele=2)
        # Fst is compaared with result from Fstat.
        #
        #For locus : loc_0_
        #Allele    Capf	 Theta	Smallf	 Relat	Relatc	 Sig_a	 Sig_b	 Sig_w
        #        1	-0.103	 0.102	-0.229	 0.228
        #        2	-0.103	 0.102	-0.229	 0.228
        #    All	-0.103	 0.102	-0.229	 0.228	 0.373	 0.051	-0.102	 0.550
        Stat(pop, Fst=[0])
        assert abs(pop.dvars().Fis[0] - (-0.229)) < 0.01
        assert abs(pop.dvars().Fit[0] - (-0.103)) < 0.01
        assert abs(pop.dvars().Fst[0] - (0.102)) < 0.01

    def testHaploFreq(self):
        'Testing calculation of haplotype frequency'
        # test haplotype frequency
        pop = population(subPop=[5000,1000], ploidy=2, loci = [10])
        if alleleType() == 'binary':
            InitByValue(pop, value=[[0]*10,[1]*10], proportions=[.3,.7])
            Stat(pop, haploFreq=[[0,1,5],[2,5]])
            assert abs(pop.dvars().haploFreq['0-1-5']['0-0-0'] - 0.3) < 0.05
            assert abs(pop.dvars().haploFreq['0-1-5']['1-1-1'] - 0.7) < 0.05
            assert abs(pop.dvars().haploFreq['2-5']['0-0'] - 0.3) < 0.05
            assert abs(pop.dvars().haploFreq['2-5']['1-1'] - 0.7) < 0.05
        else:
            InitByValue(pop, value=[[1]*10,[2]*10,[3]*10], 
                proportions=[.2,.3,.5])
            Stat(pop, haploFreq=[[0,1,5],[2,5]])
            assert abs(pop.dvars().haploFreq['0-1-5']['1-1-1'] - 0.2) < 0.05
            assert abs(pop.dvars().haploFreq['0-1-5']['2-2-2'] - 0.3) < 0.05
            assert abs(pop.dvars().haploFreq['0-1-5']['3-3-3'] - 0.5) < 0.05
            assert abs(pop.dvars().haploFreq['2-5']['1-1'] - 0.2) < 0.05
            assert abs(pop.dvars().haploFreq['2-5']['2-2'] - 0.3) < 0.05
            assert abs(pop.dvars().haploFreq['2-5']['3-3'] - 0.5) < 0.05

    def TestLD(self, freq):
        'Testing calculation of LD for a particular freq'
        #TurnOnDebug(DBG_STATOR)
        pop = population(subPop=[500,100,1000], 
            ploidy=2, loci = [5])
        InitByFreq(pop, freq)
        if alleleType() == 'binary':
            Stat(pop, LD=[2,4], haploFreq=[2,4], alleleFreq=[2,4], numOfAlleles=[2,4])
        else:
            Stat(pop, LD=[2,4], haploFreq=[2,4], alleleFreq=[2,4], numOfAlleles=[2,4])
        def LD_single(var, loc1, loc2, allele1, allele2):
            p = var.alleleFreq[loc1][allele1]
            q = var.alleleFreq[loc2][allele2]
            pq = var.haploFreq['%d-%d' % (loc1, loc2)]['%d-%d' % (allele1, allele2)]
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
            pq = var.haploFreq['%d-%d' % (loc1, loc2)]['%d-%d' % (allele1, allele2)]
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
            pq = var.haploFreq['%d-%d' % (loc1, loc2)]['%d-%d' % (allele1, allele2)]
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
        def Delta2(var, loc1, loc2):
            P11 = var.haploFreq['%d-%d' % (loc1, loc2)]['0-0']
            P12 = var.haploFreq['%d-%d' % (loc1, loc2)]['0-1']
            P21 = var.haploFreq['%d-%d' % (loc1, loc2)]['1-0']
            P22 = var.haploFreq['%d-%d' % (loc1, loc2)]['1-1']
            p1 = var.alleleFreq[loc1][0]
            q1 = var.alleleFreq[loc2][0]
            return (P11*P22-P12*P21)**2/(p1*(1-p1)*q1*(1-q1))
        #import simuUtil
        #simuUtil.ListVars(pop.dvars())
        assert (LD(pop.dvars(), 2, 4) - pop.dvars().LD[2][4]) < 1e-6
        assert (LD_prime(pop.dvars(), 2, 4) - pop.dvars().LD_prime[2][4]) < 1e-6
        assert (R2(pop.dvars(), 2, 4) - pop.dvars().R2[2][4]) < 1e-6
        if pop.dvars().numOfAlleles[2] > 2 or pop.dvars().numOfAlleles[4] > 2 :
            assert not pop.vars().has_key('Delta2')
        else:
            assert (Delta2(pop.dvars(), 2, 4) - pop.dvars().Delta2[2][4]) < 1e-6
        for sp in range(3):
            assert (LD(pop.dvars(sp), 2, 4) - pop.dvars(sp).LD[2][4]) < 1e-6
            assert (LD_prime(pop.dvars(sp), 2, 4) - pop.dvars(sp).LD_prime[2][4]) < 1e-6
            assert (R2(pop.dvars(sp), 2, 4) - pop.dvars(sp).R2[2][4]) < 1e-6
            if pop.dvars(sp).numOfAlleles[2] > 2 or pop.dvars(sp).numOfAlleles[4] > 2 :
                assert not pop.vars(sp).has_key('Delta2')
                #self.assertRaises(exceptions.AttributeError, pop.dvars(sp).Delta2)
            else:
                assert (Delta2(pop.dvars(sp), 2, 4) - pop.dvars(sp).Delta2[2][4]) < 1e-6
        #
        # test for single allele cases
        # for binary alleles, LD should be the same 
        # for standard or long alleles, LD should be different from average LD
        Stat(pop, LD=[1,3,0,1], haploFreq=[0,1], midValues=True)
        assert (abs(LD_single(pop.dvars(), 1, 3, 0, 1)) - abs(pop.dvars().ld['1-3']['0-1'])) < 1e-6
        assert (abs(LD_prime_single(pop.dvars(), 1, 3, 0, 1)) - abs(pop.dvars().ld_prime['1-3']['0-1'])) < 1e-6
        assert (abs(R2_single(pop.dvars(), 1, 3, 0, 1)) - abs(pop.dvars().r2['1-3']['0-1'])) < 1e-6
        if pop.dvars().numOfAlleles[1] > 2 or pop.dvars().numOfAlleles[3] > 2 :
            assert not pop.vars().has_key('delta2')
        else:
            assert (abs(Delta2(pop.dvars(), 1, 3)) - abs(pop.dvars().delta2['1-3']['0-1'])) < 1e-6
        for sp in range(3):
            assert (abs(LD_single(pop.dvars(sp), 1, 3, 0, 1)) - abs(pop.dvars(sp).ld['1-3']['0-1'])) < 1e-6
            assert (abs(LD_prime_single(pop.dvars(sp), 1, 3, 0, 1)) - abs(pop.dvars(sp).ld_prime['1-3']['0-1'])) < 1e-6
            assert (abs(R2_single(pop.dvars(sp), 1, 3, 0, 1)) - abs(pop.dvars(sp).r2['1-3']['0-1'])) < 1e-6
            if pop.dvars(sp).numOfAlleles[1] > 2 or pop.dvars(sp).numOfAlleles[3] > 2:
                assert not pop.vars().has_key('delta2')
            else:
                assert (abs(Delta2(pop.dvars(sp), 1, 3)) - abs(pop.dvars(sp).delta2['1-3']['0-1'])) < 1e-6
                
        
    def testLD(self):
        '''Testing LD for both dialleleic and multi-allelic cases'''
        self.TestLD([.2, .8])
        self.TestLD([.2, .3, .5])
                
                
if __name__ == '__main__':
    unittest.main()
