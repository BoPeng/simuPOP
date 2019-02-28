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

import unittest, os, sys
from simuOpt import setOptions
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

class TestPenetrance(unittest.TestCase):

    def setUp(self):
        self.pop = Population(size=[500,100,1000],
            ploidy=2, loci = [1])
        self.pop.setVirtualSplitter(RangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        for genotype, subPop in zip(
            [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
            initGenotype(self.pop, genotype=genotype, subPops=[subPop])

    def testMapPenetrance(self):
        'Testing map penetrance'
        mapPenetrance(self.pop, loci = 0,
            penetrance={(0,0):0, (0,1):1, (1,1):1})
        stat(self.pop, numOfAffected=1, vars=['numOfAffected', 'numOfAffected_sp'])
        self.assertEqual(self.pop.dvars().numOfAffected, 1425)
        self.assertEqual(self.pop.dvars(0).numOfAffected, 375)
        self.assertEqual(self.pop.dvars(1).numOfAffected, 50)
        self.assertEqual(self.pop.dvars(2).numOfAffected, 1000)
        #
        # imcomlete penetrance
        mapPenetrance(self.pop, loci = 0,
            penetrance={(0,0):0, (0,1):.3, (1,1):.5})
        stat(self.pop, numOfAffected=1, vars=['numOfAffected', 'numOfAffected_sp'])
        # self.assertTrue(abs(self.pop.dvars().numOfAffected - 880*0.3 - 545*0.5) < 100, 
        #     "Expression abs(self.pop.dvars().numOfAffected - 880*0.3 - 545*0.5) (test value %f) be less than 100. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars().numOfAffected - 880*0.3 - 545*0.5)))
        # self.assertTrue(abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5) < 30, 
        #     "Expression abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5) (test value %f) be less than 30. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5)))
        # self.assertTrue(abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5) < 15, 
        #     "Expression abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5) (test value %f) be less than 15. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5)))
        # self.assertTrue(abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5) < 50, 
        #     "Expression abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5) (test value %f) be less than 50. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5)))



    def testNoInfoField(self):
        'Testing info field for penetrance opeartors'
        pop = Population(size=[500,100,1000],
            ploidy=2, loci = [1], infoFields=['penetrance'])
        pop.setVirtualSplitter(RangeSplitter([[0,125], [125, 375], [375, 500],
            [0, 50], [50, 80], [80, 100],
            [0, 100],[100, 600], [600, 1000]]))
        for genotype, subPop in zip(
            [[0,0],[0,1],[1,1],[0,0],[0,1],[1,1],[0,1],[0,1],[1,1]],
            [(0, 0), (0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]):
            initGenotype(pop, genotype=genotype, subPops=[subPop])
        #
        mapPenetrance(pop, loci = 0,
            penetrance={(0,0):0, (0,1):1, (1,1):1},
            infoFields=['penetrance'])
        stat(pop, numOfAffected=1)
        self.assertEqual(pop.dvars().numOfAffected, 1425)
        self.assertEqual(sum(pop.indInfo('penetrance')), 1425)



    def testMaPenetrance(self):
        'Testing multi-allele penetrance'
        maPenetrance(self.pop, loci = 0, wildtype=0, penetrance=[0, 1, 1])
        stat(self.pop, numOfAffected=True, vars=['numOfAffected', 'numOfAffected_sp'])
        self.assertEqual(self.pop.dvars().numOfAffected, 1425)
        self.assertEqual(self.pop.dvars(0).numOfAffected, 375)
        self.assertEqual(self.pop.dvars(1).numOfAffected, 50)
        self.assertEqual(self.pop.dvars(2).numOfAffected, 1000)
        #
        # imcomlete penetrance
        self.pop.dvars().clear()
        maPenetrance(self.pop, loci = 0,    wildtype=0,
            penetrance=[0, .3, .5])
        stat(self.pop, numOfAffected=1, vars=['numOfAffected', 'numOfAffected_sp'])
        # self.assertTrue(abs(self.pop.dvars().numOfAffected - 880*0.3 - 545*0.5) < 100, 
        #     "Expression abs(self.pop.dvars().numOfAffected - 880*0.3 - 545*0.5) (test value %f) be less than 100. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars().numOfAffected - 880*0.3 - 545*0.5)))
        # self.assertTrue(abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5) < 30, 
        #     "Expression abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5) (test value %f) be less than 30. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(0).numOfAffected - 250*0.3 - 125*0.5)))
        # self.assertTrue(abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5) < 15, 
        #     "Expression abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5) (test value %f) be less than 15. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(1).numOfAffected - 30*0.3 - 20*0.5)))
        # self.assertTrue(abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5) < 50, 
        #     "Expression abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5) (test value %f) be less than 50. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(2).numOfAffected - 600*0.3 - 400*0.5)))


    def testMultiLocusmaPenetrance(self):
        'Testing the multi-locus version of MaPenetrance'
        pop = Population(1000, loci=[3,5], infoFields=['penetrance'])
        initGenotype(pop, freq=[.3, .7])
        #
        maPenetrance(pop, loci=[3,5], wildtype=0,
            penetrance=[0, .3, .5, 0.3, 0.6, 0.8, 0.1, 1, 0.8])

    def testMlPenetrance(self):
        'Testing multi-locus penetrance'
        pop = Population(1000, loci=[3,5], infoFields=['penetrance'])
        initGenotype(pop, freq=[.3, .7])
        #
        mlPenetrance(pop, [
            MaPenetrance(loci = 0,    wildtype=0,
                penetrance=[0, .3, .5]),
            MapPenetrance(loci = 1,
                penetrance={(0,0):0, (0,1):1, (1,1):1})
            ],
            mode=ADDITIVE
        )
        #
        mlPenetrance(pop, [
            MaPenetrance(loci = 2,    wildtype=0,
                penetrance=[0, .3, .5]),
            MapPenetrance(loci = 4,
                penetrance={(0,0):0, (0,1):1, (1,1):1})
            ],
            mode=MULTIPLICATIVE
        )

    def testEvolveMaPenetrance(self):
        'Testing using MaPenetrance in evolve function'
        pop = Population(1000, loci=[3,5])
        pop.evolve(
            initOps=[
                InitSex(),
                InitGenotype(freq=[.3, .7]),
            ],
            #
            preOps= MaPenetrance(loci = 0,    wildtype=0,
                    penetrance=[0, .3, .5]),
            matingScheme=RandomMating(),
            gen=4
        )

    def testEvolveMlPenetrance(self):
        'Testing using MlPenetrance in evolve function'
        pop = Population(1000, loci=[3,5])
        pop.evolve(
            initOps=[
                InitSex(),
                InitGenotype(freq=[.3, .7]),
            ],
            #
            preOps=MlPenetrance(ops=[
                MaPenetrance(loci = 0,    wildtype=0,
                    penetrance=[0, .3, .5]),
                MapPenetrance(loci = 1,
                    penetrance={(0,0):0, (0,1):1, (1,1):1})
                ],
                mode=ADDITIVE
            ),
            matingScheme=RandomMating(),
            gen=4
        )

    def testPyPenetrance(self):
        'Testing python penetrance operator'
        def pen(geno):
            if geno == (0, 0):
                return 0
            elif geno == (0, 1):
                return 0.5
            elif geno == (1, 0):
                return 0.5
            else:
                return 1
        pyPenetrance(self.pop, loci = 0, func=pen)
        stat(self.pop, numOfAffected=1, vars=['numOfAffected', 'numOfAffected_sp'])
        # self.assertTrue(abs(self.pop.dvars().numOfAffected -  880*0.5 - 545) < 100, 
        #     "Expression abs(self.pop.dvars().numOfAffected -  880*0.5 - 545) (test value %f) be less than 100. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars().numOfAffected -  880*0.5 - 545)))
        # self.assertTrue(abs(self.pop.dvars(0).numOfAffected - 250*0.5 - 125) < 30, 
        #     "Expression abs(self.pop.dvars(0).numOfAffected - 250*0.5 - 125) (test value %f) be less than 30. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(0).numOfAffected - 250*0.5 - 125)))
        # self.assertTrue(abs(self.pop.dvars(1).numOfAffected - 30*0.5 - 20) < 15, 
        #     "Expression abs(self.pop.dvars(1).numOfAffected - 30*0.5 - 20) (test value %f) be less than 15. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(1).numOfAffected - 30*0.5 - 20)))
        # self.assertTrue(abs(self.pop.dvars(2).numOfAffected - 600*0.5 - 400) < 50, 
        #     "Expression abs(self.pop.dvars(2).numOfAffected - 600*0.5 - 400) (test value %f) be less than 50. This test may occasionally fail due to the randomness of outcome." % (abs(self.pop.dvars(2).numOfAffected - 600*0.5 - 400)))

    def testAncestralPenetrance(self):
        'Testing the ancestralGen parameter... '
        # test the ancestralGen parameter
        # 0: set affection status for the current generation
        # -1: for all generation
        # otherwise: up to this level of ancestral generation

if __name__ == '__main__':
    unittest.main()
