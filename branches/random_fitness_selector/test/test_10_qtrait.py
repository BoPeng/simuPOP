#!/usr/bin/env python
#
# Purpose:
#     Testing quantitative trait.
#
# Author:
#     Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 146 $
# $LastChangedDate: 2006-02-03 00:18:04 -0600 (Fri, 03 Feb 2006) $
#

import unittest, os, sys, math, random
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

class TestQuanTrait(unittest.TestCase):

    def setUp(self):
        self.pop = Population(size=[5000],
            ploidy=2, loci = [1], infoFields=['qtrait'])
        self.pop.setVirtualSplitter(RangeSplitter([[0,1250], [1250,3750],[3750,5000]]))
        initGenotype(self.pop, genotype = [0,0], subPops = [(0, 0)])
        initGenotype(self.pop, genotype = [0,1], subPops = [(0, 1)])
        initGenotype(self.pop, genotype = [1,1], subPops = [(0, 2)])

    def stdev(self, x):
        'calculate standard devisaion'
        mean = sum(x)/len(x)
        return math.sqrt( sum([(i-mean)**2 for i in x])/(len(x)-1) )

    def testPyQuanTrait(self):
        'Testing the hybrid quantitative trait operator'
        def qt(geno):
            if geno == [0, 0]:
                return random.normalvariate(0, 0.5)
            elif geno == [0, 1]:
                return random.normalvariate(0.5, 1)
            elif geno == [1, 0]:
                return random.normalvariate(0.5, 1)
            else:
                return random.normalvariate(1, 2)
        pyQuanTrait(self.pop, loci=[0], func=qt, infoFields='qtrait')
        #
        # multi-locus
        pop = Population(1000, loci=[3,5], infoFields=['qtrait'])
        initGenotype(pop, freq=[.3, .7])
        def qt1(geno):
            self.assertEqual(len(geno), 4)
            return random.normalvariate(0, 0.5*sum(geno) )
        pyQuanTrait(pop, loci=[2,6], func=qt1, infoFields='qtrait')
        # multi-fields
        pop = Population(1000, loci=[3,5], infoFields=['qtrait1', 'qtrait2'])
        initGenotype(pop, freq=[.3, .7])
        def qt1(geno):
            self.assertEqual(len(geno), 4)
            return random.normalvariate(0, 0.5*sum(geno) ), 1
        pyQuanTrait(pop, loci=[2,6], func=qt1, infoFields=['qtrait1', 'qtrait2'])

    def testAncestralGen(self):
        'Testing parameter ancestralGen of qtrait... (FIXME)'
        # test the ancestralGen parameter of qtrait


if __name__ == '__main__':
    unittest.main()
