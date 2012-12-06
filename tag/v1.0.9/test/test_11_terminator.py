#!/usr/bin/env python
#
# Purpose:
#  testing Pause operator
#
# Author:
#  Bo Peng (bpeng@rice.edu)
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

class TestTerminator(unittest.TestCase):

    def testTerminator(self):
        'Testing operator TerminateIf that terminates a replicate'
        pop = Population(size=100, loci=[2])
        simu = Simulator(pop, rep=5)
        gens = simu.evolve(
            initOps = [InitSex(), InitGenotype(freq=[0.3, 0.7])],
            matingScheme = RandomMating(),
            postOps = [
                Stat(alleleFreq=[0]),
                TerminateIf('alleleNum[0][0] == 0 or alleleNum[0][0] == 200')
            ]
        )
        # it is very unlikely that gens are equal
        self.assertEqual(gens[0] == gens[1] == gens[2] == gens[3] == gens[4], False)
        for pop in simu.populations():
            stat(pop, alleleFreq=[0])
            self.assertEqual(pop.dvars().alleleFreq[0][0] == 0 or \
                pop.dvars().alleleFreq[0][0] == 1, True)


    def testTerminateAll(self):
        'Testing operator TerminateIf that terminates all replicates'
        pop = Population(size=100, loci=[2])
        simu = Simulator(pop, rep=5)
        gens = simu.evolve(
            initOps = [InitSex(), InitGenotype(freq=[0.3, 0.7])],
            matingScheme = RandomMating(),
            postOps = [
                Stat(alleleFreq=[0]),
                TerminateIf('alleleNum[0][0] == 0 or alleleNum[0][0] == 200', stopAll=True)
            ]
        )
        return
        # If a previous replicate stops, all others stop evolving.
        self.assertEqual(gens[0] >= gens[1] >= gens[2] >= gens[3] >= gens[4], True)
        # which one stopped?
        pop = simu.population(4)
        for rep in range(4):
            if gens[rep+1] < gens[rep]:
                pop = simu.population(rep)
                break
        stat(pop, alleleFreq=[0])
        self.assertEqual(pop.dvars().alleleFreq[0][0] == 0 or \
            pop.dvars().alleleFreq[0][0] == 1, True)


if __name__ == '__main__':
    unittest.main()
