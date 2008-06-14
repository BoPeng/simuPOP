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
        simu = simulator(population(1000, loci=[2]),
            randomMating())
        # now if we want to flip a lot of alleles whenever it reaches 0.2
        simu.evolve(
            # init to 1,2 (or 0, 1 in the binary case)
            preOps = [ initByFreq([.5,.5]) ],
            ops = [
                # count number of allels at this locus
                stat(alleleFreq=[0]),
                # inject 50% of allele 2 if this allele get low freq
                ifElse('alleleFreq[0][0]<0.2',
                    kamMutator(rate=.6, maxAllele=1, loci=[0]) ),
                # the other way around?
                ifElse('alleleFreq[0][0]>0.8',
                    kamMutator(rate=.6, maxAllele=1, loci=[0]) ),
                # terminate if
                terminateIf('alleleFreq[0][0]<0.1 or alleleFreq[0][0]>0.9')
            ],
            end=1000
        )
        self.assertEqual(simu.gen(), 1001)

if __name__ == '__main__':
    unittest.main()



