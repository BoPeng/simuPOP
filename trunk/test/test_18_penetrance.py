#!/usr/bin/env python
#
# Purpose:
#     Testing penetrance.
#
# Author:
#     Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

class TestPenetrance(unittest.TestCase):

    def setUp(self):
        self.pop = population(size=[500,100,1000],
            ploidy=2, loci = [1])
        self.pop.setVirtualSplitter(rangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        InitByValue(self.pop, value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            subPops = [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)])

    def testMapPenetrance(self):
        'Testing map penetrance'
        MapPenetrance(self.pop, loci = 0,
            penetrance={'0-0':0, '0-1':1, '1-1':1})
        Stat(self.pop, numOfAffected=1, vars=['numOfAffected', 'numOfAffected_sp'])
        self.assertEqual(self.pop.dvars().numOfAffected, 1425)
        self.assertEqual(self.pop.dvars(0).numOfAffected, 375)
        self.assertEqual(self.pop.dvars(1).numOfAffected, 50)
        self.assertEqual(self.pop.dvars(2).numOfAffected, 1000)
        #
        # imcomlete penetrance
        MapPenetrance(self.pop, loci = 0,
            penetrance={'0-0':0, '0-1':.3, '1-1':.5})
        Stat(self.pop, numOfAffected=1, vars=['numOfAffected', 'numOfAffected_sp'])
        assert abs(self.pop.dvars().numOfAffected - 880*0.3 - 545*0.5) < 100
        assert abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5) < 30
        assert abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5) < 15
        assert abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5) < 50



    def testNoInfoField(self):
        'Testing info field for penetrance opeartors'
        pop = population(size=[500,100,1000],
            ploidy=2, loci = [1], infoFields=['penetrance'])
        pop.setVirtualSplitter(rangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        InitByValue(pop, value = [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            subPops = [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)])
        #
        MapPenetrance(pop, loci = 0,
            penetrance={'0-0':0, '0-1':1, '1-1':1},
            infoFields=['penetrance'])
        Stat(pop, numOfAffected=1)
        self.assertEqual(pop.dvars().numOfAffected, 1425)
        self.assertEqual(sum(pop.indInfo('penetrance')), 1425)



    def testMaPenetrance(self):
        'Testing multi-allele penetrance'
        MaPenetrance(self.pop, loci = 0, wildtype=0, penetrance=[0, 1, 1])
        Stat(self.pop, numOfAffected=True, vars=['numOfAffected', 'numOfAffected_sp'])
        self.assertEqual(self.pop.dvars().numOfAffected, 1425)
        self.assertEqual(self.pop.dvars(0).numOfAffected, 375)
        self.assertEqual(self.pop.dvars(1).numOfAffected, 50)
        self.assertEqual(self.pop.dvars(2).numOfAffected, 1000)
        #
        # imcomlete penetrance
        self.pop.dvars().clear()
        MaPenetrance(self.pop, loci = 0,    wildtype=0,
            penetrance=[0, .3, .5])
        Stat(self.pop, numOfAffected=1, vars=['numOfAffected', 'numOfAffected_sp'])
        assert abs(self.pop.dvars().numOfAffected -    880*0.3 - 545*0.5) < 100
        assert abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5) < 30
        assert abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5) < 15
        assert abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5) < 50


    def testMultiLocusMaPenetrance(self):
        'Testing the multi-locus version of maPenetrance'
        pop = population(1000, loci=[3,5], infoFields=['penetrance'])
        InitByFreq(pop, [.3, .7])
        #
        MaPenetrance(pop, loci=[3,5], wildtype=0,
            penetrance=[0, .3, .5, 0.3, 0.6, 0.8, 0.1, 1, 0.8])


    def testMlPenetrance(self):
        'Testing multi-locus penetrance'
        pop = population(1000, loci=[3,5], infoFields=['penetrance'])
        InitByFreq(pop, [.3, .7])
        #
        MlPenetrance(pop, [
            maPenetrance(loci = 0,    wildtype=0,
                penetrance=[0, .3, .5]),
            mapPenetrance(loci = 1,
                penetrance={'0-0':0, '0-1':1, '1-1':1})
            ],
            mode=Additive
        )
        #
        MlPenetrance(pop, [
            maPenetrance(loci = 2,    wildtype=0,
                penetrance=[0, .3, .5]),
            mapPenetrance(loci = 4,
                penetrance={'0-0':0, '0-1':1, '1-1':1})
            ],
            mode=Multiplicative
        )

    def testPyPenetrance(self):
        'Testing python penetrance operator'
        def pen(geno):
            if geno == [0, 0]:
                return 0
            elif geno == [0, 1]:
                return 0.5
            elif geno == [1, 0]:
                return 0.5
            else:
                return 1
        PyPenetrance(self.pop, loci = 0, func=pen)
        Stat(self.pop, numOfAffected=1, vars=['numOfAffected', 'numOfAffected_sp'])
        assert abs(self.pop.dvars().numOfAffected -  880*0.5 - 545) < 100
        assert abs(self.pop.dvars(0).numOfAffected - 250*0.5 - 125) < 30
        assert abs(self.pop.dvars(1).numOfAffected - 30*0.5 - 20) < 15
        assert abs(self.pop.dvars(2).numOfAffected - 600*0.5 - 400) < 50

    def testAncestralPenetrance(self):
        'Testing the ancestralGen parameter... (FIXME)'
        # test the ancestralGen parameter
        # 0: set affection status for the current generation
        # -1: for all generation
        # otherwise: up to this level of ancestral generation

if __name__ == '__main__':
    unittest.main()
