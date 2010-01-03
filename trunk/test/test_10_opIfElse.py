#!/usr/bin/env python
#
# This is a unittest file for ifElse operator
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

class TestIfElseOperator(unittest.TestCase):

    def testIfElseOperator(self):
        'Testing opeartor ifElse'
        pop = population(1000, loci=[2])
        # now if we want to flip a lot of alleles whenever it reaches 0.2
        pop.evolve(
            # init to 1,2 (or 0, 1 in the binary case)
            initOps = [ initSex(), initByFreq([.5,.5]) ],
            matingScheme = randomMating(),
            postOps = [
                # count number of allels at this locus
                stat(alleleFreq=[0]),
                # inject 50% of allele 2 if this allele get low freq
                ifElse('alleleFreq[0][0]<0.2',
                    KamMutator(rates=.6, k=2, loci=[0]) ),
                # the other way around?
                ifElse('alleleFreq[0][0]>0.8',
                    KamMutator(rates=.6, k=2, loci=[0]) ),
                # terminate if
                TerminateIf('alleleFreq[0][0]<0.1 or alleleFreq[0][0]>0.9')
            ],
            gen = 1000
        )
        self.assertEqual(pop.dvars().gen, 1000)

    def testIfElseOperators(self):
        'Testing opeartor ifElse with multiple operators'
        simu = simulator(population(1000, loci=[2]))
        # now if we want to flip a lot of alleles whenever it reaches 0.2
        simu.evolve(
            # init to 1,2 (or 0, 1 in the binary case)
            initOps = [ initSex(), initByFreq([.5,.5]) ],
            matingScheme = randomMating(),
            postOps = [
                # count number of allels at this locus
                stat(alleleFreq=[0]),
                # inject 50% of allele 2 if this allele get low freq
                ifElse('alleleFreq[0][0]<0.2', ifOps=[
                    KamMutator(rates=.6, k=2, loci=[0]),
                    NoneOp(),
                    ]),
                # the other way around?
                ifElse('alleleFreq[0][0]>0.8', ifOps = [
                    KamMutator(rates=.6, k=2, loci=[0]),
                    NoneOp()
                    ]),
                # terminate if
                TerminateIf('alleleFreq[0][0]<0.1 or alleleFreq[0][0]>0.9')
            ],
            gen = 1000
        )
        self.assertEqual(simu.dvars(0).gen, 1000)


if __name__ == '__main__':
    unittest.main()



