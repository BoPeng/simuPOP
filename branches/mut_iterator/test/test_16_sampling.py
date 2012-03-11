#!/usr/bin/env python
#
# Purpose:
#
#     Unittest testcases for ascertainment operators
#
# Bo Peng (bpeng@rice.edu)
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
from simuPOP.utils import *
from simuPOP.sampling import *

class TestSampling(unittest.TestCase):

    def setUp(self):
        self.pop = Population(size=[1000,2000], loci=[5,10],
                ancGen=1,
                infoFields=['father_idx', 'mother_idx', 'father_id', 'mother_id', 'ind_id', 'oldindex'])
        self.pop.evolve(
            initOps = [
                 InitSex(),
                 InitGenotype(freq=[.2, .8]),
                 IdTagger(),
            ],
            matingScheme=RandomMating(
                numOffspring=(UNIFORM_DISTRIBUTION, 2, 4),
                ops=[MendelianGenoTransmitter(),
                IdTagger(),
                PedigreeTagger(),
                ParentsTagger(),
            ]),
            postOps = [
                MapPenetrance(loci=0,
                    penetrance={(0,0):0.1,(0,1):.7,(1,1):1}),
            ],
            gen = 4
        )
        for gen in range(self.pop.ancestralGens(), -1, -1):
            self.pop.useAncestralGen(gen)
            self.pop.setIndInfo(list(range(self.pop.popSize())), 'oldindex')
        # more complicated one
        self.largepop = Population(size=[5000, 20000], ploidy=2, loci=[5,10],
                ancGen=2, infoFields=['fitness', 'father_idx', 'mother_idx', 'migrate_to', 'oldindex', 'father_id', 'mother_id', 'ind_id'])
        self.largepop.evolve(
            initOps=[
                 InitSex(),
                 InitGenotype(freq=[.2, .8], loci=[0]),
                 InitGenotype(freq=[.2]*5, loci=list(range(1, self.largepop.totNumLoci()))),
                 IdTagger(),
            ],
            #preOps = Migrator(rate=[[0.1,0.1],[0.1,0.1]]),
            matingScheme=RandomMating(
                numOffspring=(UNIFORM_DISTRIBUTION, 2, 5),
                ops=[MendelianGenoTransmitter(),
                IdTagger(),
                PedigreeTagger(),
                ParentsTagger(),
                ]),
            postOps = [
                Stat( alleleFreq=[0,1], genoFreq=[0,1]),
                MapPenetrance(loci=0,
                    penetrance={(0,0):0.1, (0,1):.7, (1,1):1}),
            ],
            gen = 10
        )
        for gen in range(self.largepop.ancestralGens(), -1, -1):
            self.largepop.useAncestralGen(gen)
            self.largepop.setIndInfo(list(range(self.largepop.popSize())), 'oldindex')

    def testRandomSample(self):
        'Testing random sampling (incomplete)'
        #
        s = drawRandomSample(self.pop, 10)
        self.assertEqual(s.popSize(), 10)
        for ind in s.individuals():
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        # 
        s = drawRandomSample(self.pop, [2, 8])
        self.assertEqual(s.subPopSize(0), 2)
        self.assertEqual(s.subPopSize(1), 8)
        #
        for ind in s.individuals():
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        #
        self.pop.setVirtualSplitter(SexSplitter())
        s = drawRandomSample(self.pop, 10, subPops=[(0,0), (1,0)])
        # all samples should be MALE
        self.assertEqual(s.popSize(), 10)
        for ind in s.individuals():
            self.assertEqual(ind.sex(), MALE)
        #
        s = drawRandomSample(self.pop, [2, 8], subPops=[(0,0), (1,0)])
        # all samples should be MALE
        self.assertEqual(s.subPopSizes(), (2, 8))
        for ind in s.individuals():
            self.assertEqual(ind.sex(), MALE)
        #
        samples = drawRandomSamples(self.pop, [2, 8], subPops=[(0,0), (1,0)], numOfSamples=10)
        self.assertEqual(len(samples), 10)
        for s in samples:
            # all samples should be MALE
            self.assertEqual(s.subPopSizes(), (2, 8))
            for ind in s.individuals():
                self.assertEqual(ind.sex(), MALE)

    def testCaseControlSample(self):
        'Testing case control sampling (incomplete)'
        # case control sampling.
        s = drawCaseControlSample(self.pop, 10, 10)
        self.assertEqual(s.popSize(), 20)
        #
        s = drawCaseControlSample(self.pop, cases=[1,2], controls=[5,4])
        self.assertEqual(s.popSize(), 12)
        self.assertEqual(s.subPopSize(0), 6)
        self.assertEqual(s.subPopSize(1), 6)
        # # old index
        self.assertEqual(s.hasInfoField('oldindex'), True)
        #
        for ind in s.individuals():
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        #
        # draw from vsp?


    def testAffectedSibpairSample(self):
        'Testing affected sibpair sampling (incomplete)'
		# FIXME: testing sharing of parents
		# (father_idx and mother_idx of original and sample Population,
		# and if the parents are the same.)
        s = drawAffectedSibpairSample(self.pop, families = [2, 3])
        self.assertEqual(s.subPopSize(0), 4)
        self.assertEqual(s.subPopSize(1), 6)
        for ind in s.individuals(0):
            self.assertEqual(ind.affected(), True)
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        for ind in s.individuals(1):
            self.assertEqual(ind.affected(), True)
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        #
        s = drawAffectedSibpairSample(self.pop, 2)
        s.useAncestralGen(1)
        self.assertEqual(s.popSize(), 4)
        s.useAncestralGen(0)
        self.assertEqual(s.popSize(), 4)
        for ind in s.individuals():
            self.assertEqual(ind.affected(), True)
            inpop = self.pop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)


    def testThreeGenFamilySample(self):
        'Testing large pedigree sampling (FIXME)'
        s = drawThreeGenFamilySample(self.largepop, families=10, pedSize=(3, 20),
            numOffspring=(1,5), numOfAffected=(0, 5))
        for ind in s.individuals():
            #old index?
            inpop = self.largepop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)
        s = drawThreeGenFamilySample(self.largepop, families=50, pedSize=(3, 20),
            numOffspring=(1, 5), numOfAffected=(0, 5))
        for ind in s.individuals():
            #old index?
            inpop = self.largepop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)


    def testNuclearFamilySample(self):
        'Testing nuclear family sampling (incomplete)'
        s = drawNuclearFamilySample(self.largepop, families=50, numOffspring=(1,5),
            affectedParents=(0, 2), affectedOffspring=(1,5))
        for ind in s.individuals():
            #old index?
            inpop = self.largepop.individual(ind.oldindex)
            self.assertEqual(ind, inpop)
        #
        s = drawNuclearFamilySample(self.largepop, families=[2, 3],
            numOffspring=(1,5), affectedParents=(0,2), affectedOffspring=(1,5))
        for ind in s.individuals():
            #old index?
            inpop = self.largepop.individual(int(ind.oldindex))
            self.assertEqual(ind, inpop)



if __name__ == '__main__':
    unittest.main()
