#!/usr/bin/env python
#
# Purpose:
#  testing pause operator
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

class TestTerminator(unittest.TestCase):

    def testTerminator(self):
        'Testing operator terminateIf that terminates a replicate'
        pop = population(size=100, loci=[2])
        simu = simulator(pop, randomMating(), rep=5)
        gens = simu.evolve(
            preOps = [initByFreq([0.3, 0.7])],
            ops = [
                stat(alleleFreq=[0]),
                terminateIf('alleleNum[0][0] == 0 or alleleNum[0][0] == 200')
            ]
        )
        # it is very unlikely that gens are equal
        self.assertEqual(gens[0] == gens[1] == gens[2] == gens[3] == gens[4], False)
        for pop in simu.populations():
            Stat(pop, alleleFreq=[0])
            self.assertEqual(pop.dvars().alleleFreq[0][0] == 0 or \
                pop.dvars().alleleFreq[0][0] == 1, True)


    def testTerminateAll(self):
        'Testing operator terminateIf that terminates all replicates'
        pop = population(size=100, loci=[2])
        simu = simulator(pop, randomMating(), rep=5)
        gens = simu.evolve(
            preOps = [initByFreq([0.3, 0.7])],
            ops = [
                stat(alleleFreq=[0]),
                terminateIf('alleleNum[0][0] == 0 or alleleNum[0][0] == 200', stopAll=True)
            ]
        )
        # If a previous replicate stops, all others stop evolving.
        self.assertEqual(gens[0] >= gens[1] >= gens[2] >= gens[3] >= gens[4], True)
        # which one stopped?
        pop = simu.population(4)
        for rep in range(4):
            if gens[rep+1] < gens[rep]:
                pop = simu.population(rep)
                break
        Stat(pop, alleleFreq=[0])
        self.assertEqual(pop.dvars().alleleFreq[0][0] == 0 or \
            pop.dvars().alleleFreq[0][0] == 1, True)


if __name__ == '__main__':
  unittest.main()
