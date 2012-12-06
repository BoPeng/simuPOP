#!/usr/bin/env python
#
# Purpose:
#
# This is a unittest file for individual object
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

class TestIndividual(unittest.TestCase):

    # define a function to create basic populations
    def getPop(self):
        if moduleInfo()['alleleType'] != 'binary':
            pop = Population(size=[20, 80], ploidy=2, loci=[5, 7],
                lociPos=[2, 3, 4, 5, 6, 2, 4, 6, 8, 10, 12, 14],
                alleleNames=['_', 'A', 'C', 'T', 'G'])
        else: # binary
            pop = Population(size=[20, 80], ploidy=2, loci=[5, 7],
                lociPos=[2, 3, 4, 5, 6, 2, 4, 6, 8, 10, 12, 14],
                alleleNames=['1', '2'])
        return pop

    def testAllele(self):
        'Testing individual::Allele(idx), Allele(idx, p), Allele(idx, p, chrom), '
        'setAllele(allele, idx), setAllele(allele, idx, p), setAllele(allele, idx, p, chrom)'
        pop = self.getPop()
        ind = pop.individual(0)
        gt = ind.genotype()
        if moduleInfo()['alleleType'] == 'binary':
            gt[:] = [0, 1]*12
        else:
            gt[:] = [2, 3, 4]*8
        if moduleInfo()['alleleType'] == 'binary':
            # layout 0 1 0 1 0 | 1 0 1 0 1 0 1 || 0 1 0 1 0 | 1 0 1 0 1 0 1
            self.assertEqual(ind.allele(0), 0)
            self.assertEqual(ind.allele(1), 1)
            self.assertEqual(ind.allele(5), 1)
            self.assertEqual(ind.allele(1, 1), 1)
            self.assertEqual(ind.allele(2, 1), 0)
            self.assertEqual(ind.allele(2, 1, 1), 1)
            self.assertEqual(ind.allele(2, 0, 1), 1)
            self.assertRaises(IndexError, ind.allele, 24)
            self.assertRaises(IndexError, ind.allele, 12, 0)
            self.assertRaises(IndexError, ind.allele, 5, 0, 0)
            self.assertRaises(IndexError, ind.allele, 0, 2)
            self.assertRaises(IndexError, ind.allele, 0, 0, 2)
        else:
            # layout 2 3 4 2 3 | 4 2 3 4 2 3 4 || 2 3 4 2 3 | 4 2 3 4 2 3 4
            self.assertEqual(ind.allele(0), 2)
            self.assertEqual(ind.allele(1), 3)
            self.assertEqual(ind.allele(5), 4)
            self.assertEqual(ind.allele(1, 1), 3)
            self.assertEqual(ind.allele(2, 1), 4)
            self.assertEqual(ind.allele(2, 1, 1), 3)
            self.assertEqual(ind.allele(2, 0, 1), 3)
            self.assertRaises(IndexError, ind.allele, 24)
            self.assertRaises(IndexError, ind.allele, 12, 0)
            self.assertRaises(IndexError, ind.allele, 5, 0, 0)
            self.assertRaises(IndexError, ind.allele, 0, 2)
            self.assertRaises(IndexError, ind.allele, 0, 0, 2)
        if moduleInfo()['alleleType'] == 'binary':
            # set allele
            # layout 0 1 0 1 0 | 1 0 1 0 1 0 1 || 0 1 0 1 0 | 1 0 1 0 1 0 1
            ind.setAllele(0, 1)
            self.assertEqual(ind.allele(0), 0)
            ind.setAllele(0, 5)
            self.assertEqual(ind.allele(5), 0)
            ind.setAllele(0, 1, 1)
            self.assertEqual(ind.allele(1, 1), 0)
            ind.setAllele(1, 2, 1)
            self.assertEqual(ind.allele(2, 1), 1)
            ind.setAllele(0, 2, 1, 1)
            self.assertEqual(ind.allele(2, 1, 1), 0)
            ind.setAllele(0, 2, 0, 1)
            self.assertEqual(ind.allele(2, 0, 1), 0)
            self.assertRaises(IndexError, ind.setAllele, 0, 24)
            self.assertRaises(IndexError, ind.setAllele, 0, 12, 0)
            self.assertRaises(IndexError, ind.setAllele, 0, 5, 0, 0)
            self.assertRaises(IndexError, ind.setAllele, 0, 0, 2)
            self.assertRaises(IndexError, ind.setAllele, 0, 0, 0, 2)
        else:
            # layout 2 3 4 2 3 | 4 2 3 4 2 3 4 || 2 3 4 2 3 | 4 2 3 4 2 3 4
            ind.setAllele(1, 1)
            self.assertEqual(ind.allele(1), 1)
            ind.setAllele(2, 1)
            self.assertEqual(ind.allele(1), 2)
            ind.setAllele(3, 5)
            self.assertEqual(ind.allele(5), 3)
            ind.setAllele(1, 1, 1)
            self.assertEqual(ind.allele(1, 1), 1)
            ind.setAllele(2, 2, 1)
            self.assertEqual(ind.allele(2, 1), 2)
            ind.setAllele(1, 2, 1, 1)
            self.assertEqual(ind.allele(2, 1, 1), 1)
            ind.setAllele(1, 2, 0, 1)
            self.assertEqual(ind.allele(2, 0, 1), 1)
            self.assertRaises(IndexError, ind.setAllele, 0, 24)
            self.assertRaises(IndexError, ind.setAllele, 0, 12, 0)
            self.assertRaises(IndexError, ind.setAllele, 0, 5, 0, 0)
            self.assertRaises(IndexError, ind.setAllele, 0, 0, 2)
            self.assertRaises(IndexError, ind.setAllele, 0, 0, 0, 2)

    def testGenotype(self):
        'Testing individual::genotype(), genotype(p), genotype(p, chrom)'
        pop = self.getPop()
        ind = pop.individual(0)
        gt = ind.genotype()
        return
        if moduleInfo()['alleleType'] == 'binary':
            ind.setGenotype([0, 1])
            self.assertEqual(gt, [0, 1]*12)
        else:
            ind.setGenotype([2, 3, 4])
            self.assertEqual(gt, [2, 3, 4]*8)
        # ploidy 1
        gt = ind.genotype(1)
        if moduleInfo()['alleleType'] == 'binary':
            self.assertEqual(gt, [0, 1]*6)
        else:
            self.assertEqual(gt, [2, 3, 4]*4)
        # ploidy 1, ch 1
        gt = ind.genotype(1, 1)
        if moduleInfo()['alleleType'] == 'binary':
            self.assertEqual(gt, [1, 0, 1, 0, 1, 0, 1])
        else:
            self.assertEqual(gt, [4, 2, 3, 4, 2, 3, 4])
        #
        self.assertRaises(IndexError, ind.genotype, 2)
        self.assertRaises(IndexError, ind.genotype, 0, 2)
        # accept single form
        pop.setGenotype(1)
        self.assertEqual(pop.individual(0).allele(0), 1)
        pop.individual(0).setGenotype(0)
        self.assertEqual(pop.individual(0).allele(0), 0)

    def testIndMutants(self):
        'Testing individual::mutants()'
        if moduleInfo()['alleleType'] == 'mutant':
            pop = self.getPop()
            ind = pop.individual(0)
            ind.setGenotype([2, 3, 4])
            gt = tuple(ind.mutants())
            self.assertEqual(gt, ((0,2),(1,3),(2,4),(3,2),(4,3),(5,4),(6,2),(7,3),(8,4),(9,2),(10,3),(11,4),(12,2),(13,3),(14,4),(15,2),(16,3),(17,4),(18,2),(19,3),(20,4),(21,2),(22,3),(23,4)))
            # ploidy 1
            gt = tuple(ind.mutants(1))
            self.assertEqual(gt, ((12,2),(13,3),(14,4),(15,2),(16,3),(17,4),(18,2),(19,3),(20,4),(21,2),(22,3),(23,4)))
            # ploidy 1, ch 1
            gt = tuple(ind.mutants(1, 1))
            self.assertEqual(gt, ((17,4),(18,2),(19,3),(20,4),(21,2),(22,3),(23,4)))
            #
            self.assertRaises(IndexError, ind.mutants, 2)
            self.assertRaises(IndexError, ind.mutants, 0, 2)
            # accept single form
            pop.setGenotype(1)
            self.assertEqual(pop.individual(0).allele(0), 1)
            pop.individual(0).setGenotype(0)
            self.assertEqual(pop.individual(0).allele(0), 0)

            #  case: Genotype has zero value
            ind.setGenotype([2, 0, 4, 0, 6])
            gt = tuple(ind.mutants())
            self.assertEqual(gt,((0, 2), (2, 4), (4, 6), (5, 2), (7, 4), (9, 6), (10, 2), (12, 4), (14, 6), (15, 2), (17, 4), (19, 6), (20, 2), (22, 4)))
            gt = tuple(ind.mutants(1))
            self.assertEqual(gt,((12, 4), (14, 6), (15, 2), (17, 4), (19, 6), (20, 2), (22, 4)))
            gt = tuple(ind.mutants(1, 1))
            self.assertEqual(gt,((17, 4), (19, 6), (20, 2), (22, 4)))


    def testSetGenotype(self):
        'Testing individual::setGenotype(geno), setgenotype(geno, p), '
        'setGenotype(geno, p, chrom)'
        pop = self.getPop()
        ind = pop.individual(0)
        gt = ind.genotype()
        if moduleInfo()['alleleType'] == 'binary':
            ind.setGenotype([0, 1])
            self.assertEqual(gt, [0, 1]*12)
        else:
            ind.setGenotype([2, 3, 4])
            self.assertEqual(gt, [2, 3, 4]*8)
        # ploidy 1
        gt = ind.genotype(1)
        if moduleInfo()['alleleType'] == 'binary':
            ind.setGenotype([0, 1], 1)
            self.assertEqual(gt, [0, 1]*6)
        else:
            ind.setGenotype([2, 3, 4], 1)
            self.assertEqual(gt, [2, 3, 4]*4)
        # ploidy 1, ch 1
        gt = ind.genotype(1, 1)
        if moduleInfo()['alleleType'] == 'binary':
            ind.setGenotype([1], 1, 1)
            self.assertEqual(gt, [1, 1, 1, 1, 1, 1, 1])
        else:
            ind.setGenotype([4, 2, 3], 1, 1)
            self.assertEqual(gt, [4, 2, 3, 4, 2, 3, 4])
        # FIXME: test the list version of parameter ploidy and chroms

    def testSex(self):
        'Testing individual::sex(), setSex(sex)'
        pop = self.getPop()
        ind = pop.individual(0)
        self.assertEqual(ind.sex(), MALE)
        ind.setSex(FEMALE)
        self.assertEqual(ind.sex(), FEMALE)

    def testAffected(self):
        'Testing individual::affected(), setAffected(affected)'
        pop = self.getPop()
        ind = pop.individual(0)
        self.assertEqual(ind.affected(), False)
        ind.setAffected(True)
        self.assertEqual(ind.affected(), True)

    def testInfo(self):
        'Testing individual::info(idx), info(name)'
        'setInfo(value, idx), setInfo(value, name)'
        pop = Population(10, infoFields=['age', 'fitness', 'trait1'])
        ind = pop.individual(0)
        ind.setInfo(2.5, 0)
        self.assertEqual(ind.info(0), 2.5)
        self.assertEqual(ind.info('age'), 2.5)
        ind.setInfo(2, 0)
        ind.setInfo(1, 'fitness')
        self.assertEqual(ind.info('fitness'), 1)
        self.assertEqual(ind.info(1), 1)

    def testMutants(self):
        'Testing indivudal.mutants()'
        pop = Population(10, loci=20)
        ind = pop.individual(5)
        ind.setGenotype([0,0,1])
        # mutants are at
        # 0 0 1=2 0 0 1=5 0 0 1 0 0 1 0 0 1 0 0 1=17 0 0
        # 1=0 0 0 1=3 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1=18 0
        mutants = list(ind.mutants())
        self.assertEqual(len(mutants), 13)
        self.assertEqual([x[0] for x in mutants], [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38])
        #
        mutants = list(ind.mutants(ploidy=1))
        self.assertEqual(len(mutants), 7)
        self.assertEqual([x[0] for x in mutants], [20, 23, 26, 29, 32, 35, 38])
        #
        ind.setAllele(0, 6, 1)
        self.assertEqual([x[0] for x in ind.mutants(ploidy=1)], [20, 23, 29, 32, 35, 38])

if __name__ == '__main__':
    unittest.main()
