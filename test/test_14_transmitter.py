#!/usr/bin/env python
#
#  This is a unittest file for operator Recombinator
#
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

class TestTransmitters(unittest.TestCase):
    def getPop(self, *args, **kwargs):
        'Create a population for testing.'
        pop = Population(*args, **kwargs)
        initSex(pop)
        initGenotype(pop, freq=[0.4] + [0.1]*6)
        return pop

    def testCloneGenoTransmitter(self):
        'Testing operator CloneGenoTransmitter()'
        pop = self.getPop(size=100, loci=[10, 20])
        applyDuringMatingOperator(CloneGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        # check if 1 is copied to 2.
        self.assertEqual(pop.individual(1).genotype(),
            pop.individual(2).genotype())
        # and 0 is not copied to 2
        self.assertNotEqual(pop.individual(0).genotype(),
            pop.individual(2).genotype())
        #
        # customized chromosomes are NOT copied
        pop = self.getPop(size=100, loci=[10, 20, 30, 30],
            chromTypes=[AUTOSOME, AUTOSOME, CUSTOMIZED, CUSTOMIZED])
        applyDuringMatingOperator(CloneGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        # check if 1 is copied to 2.
        for idx in range(2, pop.popSize()):
            for p in range(2):
                for ch in range(2):
                    self.assertEqual(pop.individual(1).genotype(p, ch),
                        pop.individual(idx).genotype(p, ch))
                for ch in range(2,4):
                    self.assertNotEqual(pop.individual(1).genotype(p, ch),
                        pop.individual(idx).genotype(p, ch))
                for ch in range(4):
                    # and 0 is not copied to 2
                    self.assertNotEqual(pop.individual(0).genotype(p, ch),
                        pop.individual(idx).genotype(p, ch))


    def testMendelianGenoTransmitter(self):
        'Testing operator MendelianGenoTransmitter()'
        pop = self.getPop(size=1000, loci=[20]*5)
        applyDuringMatingOperator(MendelianGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for idx in range(2, pop.popSize()):
            for ch in range(5):
                # check if 1 is copied to 2.
                g1 = pop.individual(idx).genotype(0, ch)
                g2 = pop.individual(idx).genotype(1, ch)
                p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
                p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
        #
        # customized chromosomes are NOT copied
        pop = self.getPop(size=1000, loci=[20]*7,
            chromTypes=[AUTOSOME]*5 + [CUSTOMIZED]*2)
        applyDuringMatingOperator(MendelianGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for idx in range(2, pop.popSize()):
            for ch in range(7):
                # check if 1 is copied to 2.
                g1 = pop.individual(idx).genotype(0, ch)
                g2 = pop.individual(idx).genotype(1, ch)
                p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
                p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
                if ch < 5:
                    self.assertEqual(g1 in p1, True)
                    self.assertEqual(g2 in p2, True)
                else:
                    self.assertNotEqual(g1 in p1, True)
                    self.assertNotEqual(g2 in p2, True)
        # MALE...
        pop = self.getPop(size=1000, loci=[20]*9,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*2)
        for idx in range(2, pop.popSize()):
            pop.individual(idx).setSex(MALE)
        applyDuringMatingOperator(MendelianGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for idx in range(2, pop.popSize()):
            for ch in range(9):
                # check if 1 is copied to 2.
                g1 = pop.individual(idx).genotype(0, ch)
                g2 = pop.individual(idx).genotype(1, ch)
                p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
                p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
                if ch < 5:
                    self.assertEqual(g1 in p1, True)
                    self.assertEqual(g2 in p2, True)
                elif ch == 5:
                    # get chrom X from mother
                    self.assertEqual(g1 in p1, True)
                elif ch == 6:
                    # get chrom Y from father
                    self.assertEqual(g2, pop.individual(0).genotype(1, ch))
                else:
                    self.assertNotEqual(g1 in p1, True)
                    self.assertNotEqual(g2 in p2, True)
        # FEMALE ...
        pop = self.getPop(size=1000, loci=[20]*9,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*2)
        for idx in range(2, pop.popSize()):
            pop.individual(idx).setSex(FEMALE)
        applyDuringMatingOperator(MendelianGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for idx in range(2, pop.popSize()):
            for ch in range(9):
                # check if 1 is copied to 2.
                g1 = pop.individual(idx).genotype(0, ch)
                g2 = pop.individual(idx).genotype(1, ch)
                p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
                p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
                if ch < 5:
                    self.assertEqual(g1 in p1, True)
                    self.assertEqual(g2 in p2, True)
                elif ch == 5:
                    # get chrom X from mother
                    self.assertEqual(g1 in p1, True)
                    # get chrom X from father
                    self.assertEqual(g2, pop.individual(0).genotype(0, ch))
                elif ch == 6:
                    # unused.
                    self.assertEqual(g1 in p1, False)
                    self.assertEqual(g2 in p2, False)
                else:
                    self.assertNotEqual(g1 in p1, True)
                    self.assertNotEqual(g2 in p2, True)


    def testSelfingGenoTransmitter(self):
        'Testing operator SelfingGenoTransmitter()'
        pop = self.getPop(size=100, loci=[20]*5)
        applyDuringMatingOperator(SelfingGenoTransmitter(),
            pop, pop, dad = 0, mom = -1, off=(2, pop.popSize()))
        for ch in range(5):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            self.assertEqual(g1 in p1, True)
            self.assertEqual(g2 in p1, True)
        #
        # customized chromosomes are NOT copied
        pop = self.getPop(size=100, loci=[20]*7,
            chromTypes=[AUTOSOME]*5 + [CUSTOMIZED]*2)
        applyDuringMatingOperator(SelfingGenoTransmitter(),
            pop, pop, dad = 0, mom = -1, off=(2, pop.popSize()))
        for ch in range(7):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p1, True)
            else:
                self.assertNotEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p1, True)
        # MALE...
        pop = self.getPop(size=100, loci=[20]*9,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*2)
        pop.individual(2).setSex(MALE)
        applyDuringMatingOperator(SelfingGenoTransmitter(),
            pop, pop, dad = 0, mom = -1, off=(2, pop.popSize()))
        for ch in range(9):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p1, True)
            elif ch == 5:
                # get chrom X from mother
                self.assertEqual(g1 in p1, True)
            elif ch == 6:
                # get chrom Y from father
                self.assertEqual(g2, pop.individual(0).genotype(1, ch))
            else:
                self.assertNotEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p1, True)
        # FEMALE ...
        pop = self.getPop(size=100, loci=[20]*9,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*2)
        pop.individual(2).setSex(FEMALE)
        applyDuringMatingOperator(SelfingGenoTransmitter(),
            pop, pop, dad = 0, mom = -1, off=(2, pop.popSize()))
        for ch in range(9):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p1, True)
            elif ch == 5:
                # get chrom X from mother
                self.assertEqual(g1 in p1, True)
                # get chrom X from father
                self.assertEqual(g2, pop.individual(0).genotype(0, ch))
            elif ch == 6:
                # unused.
                self.assertEqual(g1 in p1, False)
                self.assertEqual(g2 in p1, False)
            else:
                self.assertNotEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p1, True)


    def testHaplodiploidGenoTransmitter(self):
        'Testing operator HaplodiploidGenoTransmitter()'
        pop = self.getPop(size=100, loci=[20]*5)
        # MALE...
        pop = self.getPop(size=100, loci=[20]*9,
            chromTypes=[AUTOSOME]*7 + [CUSTOMIZED]*2)
        pop.individual(2).setSex(MALE)
        applyDuringMatingOperator(HaplodiploidGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(9):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            if ch < 7:
                self.assertEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p1, True)
            else:
                self.assertNotEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p1, True)
        # FEMALE ...
        pop = self.getPop(size=100, loci=[20]*9,
            chromTypes=[AUTOSOME]*7  + [CUSTOMIZED]*2)
        pop.individual(2).setSex(FEMALE)
        applyDuringMatingOperator(HaplodiploidGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(9):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 7:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            else:
                self.assertNotEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p2, True)
    
    def testMitochondrialGenoTransmitter(self):
        'Testing operator MitochondrialGenoTransmitter()'
        #
        pop = self.getPop(size=100, loci=[10, 20] + [20] + [30]*4,
            chromTypes=[CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*5)
        self.assertRaises(ValueError,
            applyDuringMatingOperator, MitochondrialGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        #
        pop = self.getPop(size=100, loci=[10, 20] + [30]*5,
            chromTypes=[CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*5)
        applyDuringMatingOperator(MitochondrialGenoTransmitter(),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        # 
        src = [pop.individual(1).genotype(0, ch) for ch in range(2, 7)]
        for ch in range(2, 7):
            self.assertEqual(pop.individual(2).genotype(0, ch) in src,
                    True)


    def testRecombinatorAsGenoTransmitter(self):
        'Testing operator Recombinator as a genotype transmitter.'
        # test recombine on empty population
        pop = self.getPop(size=100, loci=[0])
        # this should be allowed.
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        #
        pop = self.getPop(size=100, loci=[20]*5)
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(5):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            self.assertEqual(g1 in p1, True)
            self.assertEqual(g2 in p2, True)
        #
        # customized chromosomes are NOT copied
        pop = self.getPop(size=100, loci=[20]*7,
            chromTypes=[AUTOSOME]*5 + [CUSTOMIZED]*2)
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(7):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            else:
                self.assertNotEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p2, True)
        # MALE...
        pop = self.getPop(size=100, loci=[20]*9,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*2)
        pop.individual(2).setSex(MALE)
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(9):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            elif ch == 5:
                # get chrom X from mother
                self.assertEqual(g1 in p1, True)
            elif ch == 6:
                # get chrom Y from father
                self.assertEqual(g2, pop.individual(0).genotype(1, ch))
            else:
                self.assertNotEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p2, True)
        # FEMALE ...
        pop = self.getPop(size=100, loci=[20]*9,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*2)
        pop.individual(2).setSex(FEMALE)
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(9):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            elif ch == 5:
                # get chrom X from mother
                self.assertEqual(g1 in p1, True)
                # get chrom X from father
                self.assertEqual(g2, pop.individual(0).genotype(0, ch))
            elif ch == 6:
                # unused.
                self.assertEqual(g1 in p1, False)
                self.assertEqual(g2 in p2, False)
            else:
                self.assertNotEqual(g1 in p1, True)
                self.assertNotEqual(g2 in p2, True)
        #
        #
        # With non-zero recombination rate
        #
        pop = self.getPop(size=100, loci=[20]*5)
        applyDuringMatingOperator(Recombinator(rates=0.1),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for index in range(pop.chromEnd(4)):
            ch,loc = pop.chromLocusPair(index)
            # check if 1 is copied to 2.
            g1 = pop.individual(2).allele(loc, 0, ch)
            g2 = pop.individual(2).allele(loc, 1, ch)
            p1 = [pop.individual(1).allele(loc, p, ch) for p in range(2)]
            p2 = [pop.individual(0).allele(loc, p, ch) for p in range(2)]
            self.assertEqual(g1 in p1, True)
            self.assertEqual(g2 in p2, True)
        #
        # customized chromosomes are NOT copied
        pop = self.getPop(size=100, loci=[20]*7,
            chromTypes=[AUTOSOME]*5 + [CUSTOMIZED]*2)
        applyDuringMatingOperator(Recombinator(rates=0.1),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for index in range(pop.chromEnd(6)):
            # check if 1 is copied to 2.
            ch,loc = pop.chromLocusPair(index)
            g1 = pop.individual(2).allele(loc, 0, ch)
            g2 = pop.individual(2).allele(loc, 1, ch)
            p1 = [pop.individual(1).allele(loc, p, ch) for p in range(2)]
            p2 = [pop.individual(0).allele(loc, p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            #else:
            #    self.assertNotEqual(g1 in p1, True)
            #    self.assertNotEqual(g2 in p2, True)
        # MALE...
        pop = self.getPop(size=100, loci=[20]*9,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*2)
        pop.individual(2).setSex(MALE)
        applyDuringMatingOperator(Recombinator(rates=0.1),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for index in range(pop.chromEnd(8)):
            # check if 1 is copied to 2.
            ch,loc = pop.chromLocusPair(index)
            g1 = pop.individual(2).allele(loc, 0, ch)
            g2 = pop.individual(2).allele(loc, 1, ch)
            p1 = [pop.individual(1).allele(loc, p, ch) for p in range(2)]
            p2 = [pop.individual(0).allele(loc, p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            elif ch == 5:
                # get chrom X from mother
                self.assertEqual(g1 in p1, True)
            elif ch == 6:
                # get chrom Y from father
                self.assertEqual(g2, pop.individual(0).allele(loc, 1, ch))
            #else:
                #self.assertNotEqual(g1 in p1, True)
                #self.assertNotEqual(g2 in p2, True)
        # FEMALE ...
        pop = self.getPop(size=100, loci=[20]*9,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X, CHROMOSOME_Y] + [CUSTOMIZED]*2)
        pop.individual(2).setSex(FEMALE)
        applyDuringMatingOperator(Recombinator(rates=0.1),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for index in range(pop.chromEnd(8)):
            # check if 1 is copied to 2.
            ch,loc = pop.chromLocusPair(index)
            g1 = pop.individual(2).allele(loc, 0, ch)
            g2 = pop.individual(2).allele(loc, 1, ch)
            p1 = [pop.individual(1).allele(loc, p, ch) for p in range(2)]
            p2 = [pop.individual(0).allele(loc, p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            elif ch == 5:
                # get chrom X from mother
                self.assertEqual(g1 in p1, True)
                # get chrom X from father
                self.assertEqual(g2, pop.individual(0).allele(loc, 0, ch))
            #elif ch == 6:
                #self.assertEqual(g1 in p1, False)
                #self.assertEqual(g2 in p2, False)
            #else:
                #self.assertNotEqual(g1 in p1, True)
                #self.assertNotEqual(g2 in p2, True)
        # 
        # with only chromosome X
        #
        # male
        pop = self.getPop(size=100, loci=[20]*6,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X])
        pop.individual(2).setSex(MALE)
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        return
        for ch in range(6):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            elif ch == 5:
                # get chrom X from mother
                self.assertEqual(g1 in p1, True)
        # FEMALE ...
        pop = self.getPop(size=100, loci=[20]*6,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_X])
        pop.individual(2).setSex(FEMALE)
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(6):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            elif ch == 5:
                # get chrom X from mother
                self.assertEqual(g1 in p1, True)
                # get chrom X from father
                self.assertEqual(g2, pop.individual(0).genotype(0, ch))
        # with only chromosome Y
        # MALE...
        pop = self.getPop(size=100, loci=[20]*6,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_Y])
        pop.individual(2).setSex(MALE)
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(9):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            elif ch == 5:
                # get chrom Y from father
                self.assertEqual(g2, pop.individual(0).genotype(1, ch))
        # FEMALE ...
        pop = self.getPop(size=100, loci=[20]*6,
            chromTypes=[AUTOSOME]*5 + [CHROMOSOME_Y])
        pop.individual(2).setSex(FEMALE)
        applyDuringMatingOperator(Recombinator(rates=0),
            pop, pop, dad = 0, mom = 1, off=(2, pop.popSize()))
        for ch in range(9):
            # check if 1 is copied to 2.
            g1 = pop.individual(2).genotype(0, ch)
            g2 = pop.individual(2).genotype(1, ch)
            p1 = [pop.individual(1).genotype(p, ch) for p in range(2)]
            p2 = [pop.individual(0).genotype(p, ch) for p in range(2)]
            if ch < 5:
                self.assertEqual(g1 in p1, True)
                self.assertEqual(g2 in p2, True)
            elif ch == 5:
                # unused.
                self.assertEqual(g1 in p1, False)
                self.assertEqual(g2 in p2, False)


    def testRecRate(self):
        'Testing to see if we actually recombine at this rate '
        a1, a2 = 0, 1

        pop = Population(10000, loci=[2,3,2])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        simu = Simulator(pop)
        simu.evolve( postOps = Stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating(ops = Recombinator(rates = 0.1)),
            gen=1 )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        #print simu.dvars(0).haploFreq
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.05) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.05) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.05)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(a1,a2)] - 0.05) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(a1,a2)] - 0.05) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(a1,a2)] - 0.05)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.05) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.05) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.05)))
        # compare to the next test
        pop = Population(10000, loci=[3,4])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        simu = Simulator(pop)
        simu.evolve( postOps= Stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating(ops=Recombinator(rates = 0.4, loci=[1,3])),
            gen=1 )
        # 0.5
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30)))
        # # 0.25
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30)))
        # # no recombination.
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5)))
        #
        # alrorithm 0?
        #
        pop = Population(10000, loci=[3,10])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*13+[a2]*13)
        simu = Simulator(pop)
        simu.evolve( postOps =Stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating(ops=Recombinator(rates = 0.4, loci=[1,3, 8])),
            gen=1 )
        # 0.5
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30)))
        # # 0.25
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30)))
        # # no recombination.
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5)))
        # algorithm 2 (uniform rare), just see if it crashes
        pop = Population(10000, loci=[100])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*13+[a2]*13)
        simu = Simulator(pop)
        simu.evolve( postOps =Stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating(ops=Recombinator(rates = 1e-6)),
            gen=1 )
 
    def testConversionRate(self):
        'Testing to see if we actually convert at this rate '
        a1, a2 = 0, 1
        pop = Population(10000, loci=[3,4])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        rec = Recombinator(rates = 0.4, convMode = (NUM_MARKERS, 1, 1), loci=[1,3])
        simu = Simulator(pop)
        simu.evolve( postOps = Stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating(ops=rec),
            gen=1 )
        #print simu.dvars(0).haploFreq
        #print rec.convCounts()
        # 0.5
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30)))
        # # 0.25
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30)))
        # # r/2 recombination induced recombination
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.3) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.3) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.3)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(0,1)] - 0.2) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(0,1)] - 0.2) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(0,1)] - 0.2)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.3) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.3) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.3)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(1,0)] - 0.2) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(1,0)] - 0.2) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(1,0)] - 0.2)))
        # # copied
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.5)))
        #
        # length 2
        pop = Population(100000, loci=[3,4])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        simu = Simulator(pop)
        simu.evolve( postOps = Stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating( ops=Recombinator(rates = 0.4, convMode = (NUM_MARKERS, 1, 2), loci=[1,3])),
            gen=1 )
        #
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30)))
        # # 0.25
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30)))
        # # copied
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.5)))
        # # r/2 recombination induced recombination
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.3) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.3) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.3)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(0,1)] - 0.2) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(0,1)] - 0.2) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(0,1)] - 0.2)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.3) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.3) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.3)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(1,0)] - 0.2) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(1,0)] - 0.2) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(1,0)] - 0.2)))
        #
        # algorithm 0??
        pop = Population(10000, loci=[3,10])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*13+[a2]*13)
        simu = Simulator(pop)
        simu.evolve(
            postOps = Stat( haploFreq = [[0,1], [1,2], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating(ops=Recombinator(rates = 0.4, convMode=(NUM_MARKERS, 1, 2), loci=[1,3,8]) ),
            gen=1 )
        # #
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,1)][(1,1)] - 0.5)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(1,2)][(1,1)] - 0.30)))
        # # 0.25
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(0,1)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,0)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(2,3)][(1,1)] - 0.25)))
        # # r/2
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,1)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,0)] - 0.20)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(0,0)] - 0.30)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,4)][(1,1)] - 0.30)))
        # # copied
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(0,0)] - 0.5)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.5) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.5) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(4,5)][(1,1)] - 0.5)))
        # # r/2 recombination induced recombination
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.3) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.3) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(0,0)] - 0.3)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(0,1)] - 0.2) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(0,1)] - 0.2) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(0,1)] - 0.2)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.3) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.3) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(1,1)] - 0.3)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(1,0)] - 0.2) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(5,6)][(1,0)] - 0.2) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(5,6)][(1,0)] - 0.2)))


    def testAtLociRecRates(self):
        'Testing loci parameter'
        if moduleInfo()['alleleType'] == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        pop = Population(10000, loci=[2,3,2])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        simu = Simulator(pop)
        simu.evolve(
            postOps = Stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating(ops=Recombinator(rates = [0.1,0.15,0.3], loci=[0,2,5] )),
            gen=1 )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        #print simu.dvars(0).haploFreq
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05) < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.075) < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.075) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.075)))
        try:    # do not have this haplotype
            simu.dvars(0).haploFreq[(3,4)][(a1,a2)]
        except KeyError:
            pass
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25) < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.15) < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.15) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.15)))
        #
        # test the case with loci= ALL_AVAIL
        pop = Population(10000, loci=[2,2])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*4+[a2]*4)
        # for ALL_AVAIL, rates should have a length of totNumLoci
        pop.evolve(
            postOps = Stat( haploFreq = [[0,1], [2,3]]),
            matingScheme = RandomMating(ops=Recombinator(rates = [0, 0.8, 0.5, 0.8])),
            gen=1 )
        # the first two no recombine, the second two free recombine
        self.assertEqual(pop.dvars().haploFreq[(0,1)][(a1,a2)], 0)
        self.assertEqual(pop.dvars().haploFreq[(0,1)][(a2,a1)], 0)
        # self.assertTrue(abs(pop.dvars().haploFreq[(2,3)][(a1,a2)] - 0.25) < 0.01)
        # self.assertTrue(abs(pop.dvars().haploFreq[(2,3)][(a2,a1)] - 0.25) < 0.01)
        #
        # another test, using another order of rates
        initGenotype(pop, genotype=[a1]*4+[a2]*4)
        # for ALL_AVAIL, rates should have a length of totNumLoci
        pop.evolve(
            postOps = Stat( haploFreq = [[0,1], [2,3]]),
            matingScheme = RandomMating(ops=Recombinator(rates = [0.5, 0.8, 0, 0.8])),
            gen=1 )
        # the first two no recombine, the second two free recombine
        # self.assertTrue(abs(pop.dvars().haploFreq[(0,1)][(a1,a2)] - 0.25) < 0.01)
        # self.assertTrue(abs(pop.dvars().haploFreq[(0,1)][(a2,a1)] - 0.25) < 0.01)
        self.assertEqual(pop.dvars().haploFreq[(2,3)][(a1,a2)], 0)
        self.assertEqual(pop.dvars().haploFreq[(2,3)][(a2,a1)], 0)
        # a little bit more complicated
        pop = Population(10000, loci=[5,10])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*15+[a2]*15)
        # for ALL_AVAIL, rates should have a length of totNumLoci
        pop.evolve(
            postOps = Stat( haploFreq = [[1,2], [3,4], [8,9], [13, 14]]),
            matingScheme = RandomMating(ops=Recombinator(
                rates = [0]*2 + [0.5]*3 + [0]*5 + [0.5]* 5)),
            gen=1 )
        self.assertEqual(pop.dvars().haploFreq[(1,2)][(a1,a2)], 0)
        self.assertEqual(pop.dvars().haploFreq[(1,2)][(a2,a1)], 0)
        self.assertEqual(pop.dvars().haploFreq[(8,9)][(a1,a2)], 0)
        self.assertEqual(pop.dvars().haploFreq[(8,9)][(a2,a1)], 0)
        # self.assertTrue(abs(pop.dvars().haploFreq[(3,4)][(a1,a2)] - 0.25) < 0.01)
        # self.assertTrue(abs(pop.dvars().haploFreq[(3,4)][(a2,a1)] - 0.25) < 0.01)
        self.assertTrue(abs(pop.dvars().haploFreq[(13,14)][(a1,a2)] - 0.25) < 0.01)
        self.assertTrue(abs(pop.dvars().haploFreq[(13,14)][(a2,a1)] - 0.25) < 0.01)




    def testRecIntensity(self):
        'Testing recombination intensity'
        if moduleInfo()['alleleType'] == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        pop = Population(100000, loci=[2,3,2], lociPos=[0,1,0,2,4,0,4] )
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        simu = Simulator(pop)
        simu.evolve( postOps = Stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating( ops=Recombinator(intensity = 0.1) ),
            gen=1 )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        self.assertNotEqual(simu.dvars(0).haploFreq[(0,1)][(a1,a2)], 0)
        self.assertNotEqual(simu.dvars(0).haploFreq[(2,3)][(a1,a2)], 0)
        self.assertNotEqual(simu.dvars(0).haploFreq[(3,4)][(a1,a2)], 0) 
        self.assertNotEqual(simu.dvars(0).haploFreq[(4,5)][(a1,a2)], 0)
        #print simu.dvars(0).haploFreq
        # self.assertTrue((simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05) < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(0,1)][(a1,a2)] - 0.05)))
        # self.assertTrue((simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.1)  < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.1)  (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.1)   ))
        # self.assertTrue((simu.dvars(0).haploFreq[(3,4)][(a1,a2)] - 0.1)  < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(3,4)][(a1,a2)] - 0.1)  (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(3,4)][(a1,a2)] - 0.1)   ))
        # self.assertTrue((simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25) < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25)))
        # self.assertTrue((simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.2)  < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.2)  (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.2)   ))
        # lower intensity rate
        pop = Population(10000, loci=[40])
        pop.evolve( 
            matingScheme = SelfMating(ops=Recombinator(intensity = 0.00001) ),
            gen=1 )
        # test for uniform rare case
        pop = Population(100000, loci=[7], lociPos=[0,1,2,3,4,5,6], infoFields='parent_id')
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        #turnOnDebug('DBG_TRANSMITTER')
        pop.evolve( postOps = Stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme = RandomMating(ops=Recombinator(intensity = 0.001, infoFields='parent_id', output='>>test.log') ),
            gen=1 )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        #print simu.dvars(0).haploFreq
        self.assertNotEqual(pop.dvars().haploFreq[(0,1)][(a1,a2)], 0)
        self.assertNotEqual(pop.dvars().haploFreq[(2,3)][(a1,a2)], 0)
        self.assertNotEqual(pop.dvars().haploFreq[(3,4)][(a1,a2)], 0) 
        self.assertNotEqual(pop.dvars().haploFreq[(4,5)][(a1,a2)], 0)
        # self.assertTrue(abs(pop.dvars().haploFreq[(0,1)][(a1,a2)] - 0.0005) < 0.0003, 
        #     "Expression (pop.dvars().haploFreq[(0,1)][(a1,a2)] - 0.0005) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((pop.dvars().haploFreq[(0,1)][(a1,a2)])))
        # self.assertTrue(abs(pop.dvars().haploFreq[(2,3)][(a1,a2)] - 0.0005) < 0.0003, 
        #     "Expression (pop.dvars().haploFreq[(2,3)][(a1,a2)] - 0.0005)  (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((pop.dvars().haploFreq[(2,3)][(a1,a2)])   ))
        # self.assertTrue(abs(pop.dvars().haploFreq[(3,4)][(a1,a2)] - 0.0005) < 0.0003, 
        #     "Expression (pop.dvars().haploFreq[(3,4)][(a1,a2)] - 0.0005)  (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((pop.dvars().haploFreq[(3,4)][(a1,a2)])   ))
        # self.assertTrue(abs(pop.dvars().haploFreq[(4,5)][(a1,a2)] - 0.0005) < 0.0003, 
        #     "Expression (pop.dvars().haploFreq[(4,5)][(a1,a2)] - 0.0005) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((pop.dvars().haploFreq[(4,5)][(a1,a2)])))
        # self.assertTrue(abs(pop.dvars().haploFreq[(5,6)][(a1,a2)] - 0.0005) < 0.0003, 
        #     "Expression (pop.dvars().haploFreq[(5,6)][(a1,a2)] - 0.0005)  (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((pop.dvars().haploFreq[(5,6)][(a1,a2)])   ))
  


    def testRecIntensityAfterLoci(self):
        'Testing RecIntensity after loci'
        if moduleInfo()['alleleType'] == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        pop = Population(10000, loci=[2,3,2], lociPos=[0,1,0,2,4,0,4] )
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        simu = Simulator(pop)
        simu.evolve( postOps = Stat( haploFreq = [[0,1], [2,3], [3,4], [4,5], [5,6]]),
            matingScheme =  RandomMating(ops=Recombinator(intensity = 0.1, loci=[2,5]) ),
            gen=1 )
        # the supposed proportions are 1-1: 0.5-r/2, 1-2: r/2, 2-1: r/2, 2-2: 0.5-r/2
        #print simu.dvars(0).haploFreq
        self.assertEqual(simu.dvars(0).haploFreq[(0,1)].setdefault((a1,a2),0), 0)
        # self.assertTrue((simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.1)    < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.1)    (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(2,3)][(a1,a2)] - 0.1)   ))
        self.assertEqual(simu.dvars(0).haploFreq[(3,4)].setdefault((a1,a2),0), 0)
        # self.assertTrue((simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25) < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(4,5)][(a1,a2)] - 0.25)))
        # self.assertTrue((simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.2) < 0.01, 
        #     "Expression (simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.2) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % ((simu.dvars(0).haploFreq[(5,6)][(a1,a2)] - 0.2)))


    def testNoNewAllele(self):
        'Testing that no new allele appear because of recombination'
        if moduleInfo()['alleleType'] == 'binary':
            geno = [x%2 for x in [1,2,3,4,5,6,7] ]
        else:
            geno = [1,2,3,4,5,6,7]
        pop = Population(1000, loci=[2,3,2])
        initSex(pop)
        initGenotype(pop, genotype=geno)
        simu = Simulator(pop)
        simu.evolve(  gen=100, matingScheme=RandomMating(ops=Recombinator(rates = 0.4)))
        stat(simu.population(0), alleleFreq=list(range(0,7)))
        for i in range(7):
            self.assertEqual(simu.dvars(0).alleleFreq[i][geno[i]], 1.)


    def testCrossBetweenChrom(self):
        'Testing if chromsomes are crossed by default'
        if moduleInfo()['alleleType'] == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        pop = Population(100000, loci=[2,3,2])
        initSex(pop)
        initGenotype(pop, genotype=[a1]*7+[a2]*7)
        simu = Simulator(pop)
        #turnOnDebug(DBG_RECOMBINATOR)
        simu.evolve( postOps = Stat( haploFreq = [[0,1], [2,3], [2,4], [5,6], [0,2], [0,6], [3,6]]),
            matingScheme = RandomMating(ops=Recombinator(rates = 0) ),
            gen=1 )
        self.assertEqual(simu.dvars(0).haploFreq[(0,1)].setdefault((a1,a2), 0), 0.)
        self.assertEqual(simu.dvars(0).haploFreq[(2,3)].setdefault((a1,a2), 0), 0.)
        self.assertEqual(simu.dvars(0).haploFreq[(2,4)].setdefault((a1,a2), 0), 0.)
        self.assertEqual(simu.dvars(0).haploFreq[(5,6)].setdefault((a1,a2), 0), 0.)
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,2)].setdefault((a1,a2), 0) - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,2)].setdefault((a1,a2), 0) - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,2)].setdefault((a1,a2), 0) - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(0,6)].setdefault((a1,a2), 0) - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(0,6)].setdefault((a1,a2), 0) - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(0,6)].setdefault((a1,a2), 0) - 0.25)))
        # self.assertTrue(abs(simu.dvars(0).haploFreq[(3,6)].setdefault((a1,a2), 0) - 0.25) < 0.01, 
        #     "Expression abs(simu.dvars(0).haploFreq[(3,6)].setdefault((a1,a2), 0) - 0.25) (test value %f) be less than 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).haploFreq[(3,6)].setdefault((a1,a2), 0) - 0.25)))


    def testRecProportion(self):
        'Testing table 4 of H&C 3nd edition P49 '
        N = 100000
        r = 0.1
        genoDad = [[0,0],[0,0],[0,0],[0,0],[0,1],[0,1],[0,1],[1,0],[1,0],[1,1]]
        genoMom = [[0,0],[0,1],[1,0],[1,1],[0,1],[1,0],[1,1],[1,0],[1,1],[1,1]]
        prop = [ [1, 0, 0, 0], [.5, .5, 0, 0], [.5, 0, 0.5, 0],
            [0.5-r/2, r/2, r/2, 0.5-r/2], [0, 1, 0, 0], [r/2, .5-r/2, .5-r/2, r/2],
            [0, .5, 0, .5], [0, 0, 1, 0], [0, 0, .5, .5], [0, 0, 0, 1] ]
        for i in range(0, len(genoDad)):
            pop = Population(size=N, loci=[2])
            initSex(pop)
            initGenotype(pop, genotype=genoDad[i]+genoMom[i])
            simu = Simulator(pop)
            simu.evolve(postOps = Stat(haploFreq=[0,1]),
                matingScheme = RandomMating(ops=Recombinator(rates=r)),
                    gen=1)
            hf = simu.dvars(0).haploFreq[(0,1)]
            # self.assertTrue(abs(hf[(0,0)] - prop[i][0]) <= 0.01, 
            # "Expression abs(hf[(0,0)]) - prop[i][0] (test value %f) be less than or equal to 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(hf[(0,0)]) - prop[i][0]))
            # self.assertTrue(abs(hf[(0,1)] - prop[i][1]) <= 0.01, 
            # "Expression abs(hf[(0,1)]) - prop[i][1] (test value %f) be less than or equal to 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(hf[(0,1)]) - prop[i][1]))
            # self.assertTrue(abs(hf[(1,0)] - prop[i][2]) <= 0.01, 
            # "Expression abs(hf[(1,0)]) - prop[i][2] (test value %f) be less than or equal to 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(hf[(1,0)]) - prop[i][2]))
            # self.assertTrue(abs(hf[(1,1)] - prop[i][3]) <= 0.01, 
            # "Expression abs(hf[(1,1)]) - prop[i][3] (test value %f) be less than or equal to 0.01. This test may occasionally fail due to the randomness of outcome." % (abs(hf[(1,1)]) - prop[i][3]))


    def testLDDecay(self):
        'Testing formula Dn=(1-r)^n D0 '
        r = 0.1
        N = 10000
        pop = Population(size=N, loci=[2])
        # genotype 11/22, with LD=1
        if moduleInfo()['alleleType'] == 'binary':
            a1, a2 = 0, 1
        else:
            a1, a2 = 1, 2
        initSex(pop)
        initGenotype(pop, genotype=[a1,a1,a2,a2])
        simu = Simulator(pop)
        simu.evolve(
            finalOps = Stat(LD=[0,1]),
            matingScheme = RandomMating(ops=Recombinator(rates=r) ),
            gen=9)
        # check the change of LD, hopefully, the variation is not too high.
        # self.assertTrue(abs(simu.dvars(0).LD[0][1] - 0.25*(1-r)**10) < 0.02 , 
        #     "Expression abs(simu.dvars(0).LD[0][1] - 0.25*(1-r)**10) (test value %f) be less than 0.02 . This test may occasionally fail due to the randomness of outcome." % (abs(simu.dvars(0).LD[0][1] - 0.25*(1-r)**10)))

#     def testNoMaleRec(self):
#         'Testing recombination of male chromosome. This is currently wrong.'
#         # create such an situation where female has 11111/22222, male has 1111/33333
#         # we will test if 3333 is untouched under recombination
#         r = 0.1
#         N = 100
#         if moduleInfo()['alleleType'] == 'binary':
#             a1, a2, a3 = 0, 1, 1
#         else:
#             a1, a2, a3 = 1, 2, 3
#         pop = Population(size=N, loci=[2,5], sexChrom=True)
#         # male     1 3
#         # female 1 2
#         initGenotype(pop, genotype=indRange=[0,N/2], sex=[MALE]*(N/2), atPloidy=0, value=[a1]*7)
#         initGenotype(pop, genotype=indRange=[0,N/2], sex=[MALE]*(N/2), atPloidy=1, value=[a1]*2+[a3]*5)
#         initGenotype(pop, genotype=indRange=[N/2,N], sex=[FEMALE]*(N/2), value=[a1]*7+[a2]*7)
#         # now let us recombine
#         simu = Simulator(pop, RandomMating())
#         simu.evolve( [ Recombinator(rates=r) ], gen=100)
#         pop = simu.population(0)
#         #
#         for i in range( pop.popSize() ):
#             ind = pop.individual(i)
#             if ind.sex() == MALE:
#                 # check the second chromosome
#                 # arrGenotype(ploidy, chrom), no dict parameter for efficiency purpose
#                 #print ind.arrGenotype()
#                 self.assertEqual(ind.arrGenotype(1, 1), [a3]*5)
#             else:
#                 # there is no allele 3 anywhere (only non-binary...)
#                 if moduleInfo()['alleleType'] != 'binary':
#                     self.assertEqual(ind.arrGenotype().count(a3), 0)


    def testHaplodiploid(self):
        'Testing recombination in haplodiploid populations'
        pop = Population(size=[20, 20], ploidy=HAPLODIPLOID, loci=[3,5])
        simu = Simulator(pop)
        simu.evolve(
            initOps = [
                InitSex(), InitGenotype(genotype=[0]*8 + [1]*8)
                ],
            matingScheme = HaplodiploidMating(),
            gen = 1)
        # all individuals get the second copy from the first copy of male parents
        # which are all zero
        for ind in simu.population(0).individuals():
            self.assertEqual(ind.genotype(1), [0]*8)

    def testLineage(self):
        'Testing the transmission of lineage information'
        # pretend that we advance a generation
        if moduleInfo()['alleleType'] != 'lineage':
            return
        # set lingeage with ind_id field
        pop = Population(100, infoFields='ind_id', loci=10)
        tagID(pop)
        mom = pop.individual(0)
        dad = pop.individual(1)
        child = pop.individual(2)
        MendelianGenoTransmitter().transmitGenotype(mom, child, 0)
        MendelianGenoTransmitter().transmitGenotype(dad, child, 1)
        # test lineage assignment in genotype transition
        self.assertEqual(child.lineage(), (list(mom.lineage(0)) + list(dad.lineage(1))))
       
        
if __name__ == '__main__':
    unittest.main()
