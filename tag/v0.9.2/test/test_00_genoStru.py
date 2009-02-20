#!/usr/bin/env python
#
# Purpose:
#
# This is a unittest file for carray datatype, and class genoStruTrait.
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

class TestGenoStru(unittest.TestCase):
    # define a function to create basic populations
    def getPop(self):
        pop = population(size=[20, 80], ploidy=2, loci=[5, 7],
            lociPos=[ [2, 3, 4, 5, 6], [2, 4, 6, 8, 10, 12, 14]],
            alleleNames=['_', 'A', 'C', 'T', 'G'],
            infoFields=['a', 'b'])
        InitSex(pop)
        return pop

    def testLociPos(self):
        'Testing the specification of loci position'
        pop = population(loci=[4, 5], lociPos=[1, 3, 2, 4, 1, 5, 4, 3, 2])
        self.assertEqual(pop.lociPos(), (1, 2, 3, 4, 1, 2, 3, 4, 5))
        self.assertRaises(exceptions.ValueError, population, loci=[4, 5], lociPos=[1, 3, 2, 4, 1, 3, 4, 3, 2])
        #
        pop = population(loci=[4, 2], lociPos=[1, 3, 2, 4, 10, 5], lociNames=['1', '2', '3', '4', '5', '6'])
        self.assertEqual(pop.lociPos(), (1, 2, 3, 4, 5, 10))
        self.assertEqual(pop.lociNames(), ('1', '3', '2', '4', '6', '5'))

    def testGenotypeCarray(self):
        'Testing allele carray type returned by genotype'
        pop = population(size=2, loci=[2, 1])
        InitByValue(pop, [1, 2, 3])
        arr = pop.genotype()
        arr[:] = [0, 1, 2]*4
        # can print
        # expression
        if AlleleType() != 'binary':
            self.assertEqual(str(arr), "[0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]")
        else:
            self.assertEqual(str(arr), "[0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]")
        # count
        if AlleleType() != 'binary':
            self.assertEqual(arr[5], 2)
            self.assertEqual(arr.count(2), 4)
        else:
            self.assertEqual(arr[5], 1)
            self.assertEqual(arr.count(2), 0)
            self.assertEqual(arr.count(1), 8)
        arr[0] = 1
        if AlleleType() != 'binary':
            self.assertEqual(arr.count(1), 5)
        else:
            self.assertEqual(arr.count(1), 9)
        # index
        self.assertRaises(exceptions.ValueError,
            arr.index, 10)
        self.assertEqual(arr.index(1), 0)
        # can read write
        arr[3] = 3
        if AlleleType() != 'binary':
            self.assertEqual(arr[3], 3)
        else:
            self.assertEqual(arr[3], 1)
        # convert to list
        arr[:] = [0, 1, 2]*4
        if AlleleType() != 'binary':
            self.assertEqual(arr.tolist(), [0, 1, 2]*4)
        else:
            self.assertEqual(arr.tolist(), [0, 1, 1]*4)
            self.assertNotEqual(arr.tolist(), [0, 1, 2]*4)
        # convert to list
        if AlleleType() != 'binary':
            self.assertEqual(arr, [0, 1, 2]*4)
        else:
            self.assertEqual(arr, [0, 1, 1]*4)
        # slice
        arr[:] = [0, 1, 2]*4
        arr1 = arr[:3]
        if AlleleType() != 'binary':
            self.assertEqual(arr1, [0, 1, 2])
        else:
            self.assertEqual(arr1, [0, 1, 1])
        arr1 = arr[3:5]
        if AlleleType() != 'binary':
            self.assertEqual(arr1, [0, 1])
        else:
            self.assertEqual(arr1, [0, 1])
        # assign slice
        arr1[:] = 5
        # IMPORTANT NOTE that arr will also be affected
        if AlleleType() != 'binary':
            self.assertEqual(arr1, [5, 5] )
            self.assertEqual(arr, [0, 1, 2, 5, 5, 2, 0, 1, 2, 0, 1, 2])
        else:
            self.assertEqual(arr1, [1, 1] )
            self.assertEqual(arr, [0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1])
        # assign vector
        arr1[:] = [0, 0]
        self.assertEqual(arr1, [0, 0] )
        if AlleleType() != 'binary':
            self.assertEqual(arr, [0, 1, 2, 0, 0, 2, 0, 1, 2, 0, 1, 2])
        else:
            self.assertEqual(arr, [0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1])
        # assign from another part
        arr[:6] = arr[6:12]
        if AlleleType() != 'binary':
            self.assertEqual(arr, [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2])
        else:
            self.assertEqual(arr, [0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1])

    def testPloidy(self):
        'Testing genoStruTrait::Ploidy(), PloidyName()'
        pop = self.getPop()
        self.assertEqual(pop.ploidy(), 2)
        self.assertEqual(pop.ploidyName(), 'diploid')
        pop = population(size=100, ploidy=Haplodiploid, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), 'haplodiploid')
        pop = population(size=100, ploidy=1, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), 'haploid')
        pop = population(size=100, ploidy=3, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), 'triploid')
        pop = population(size=100, ploidy=4, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), 'tetraploid')
        pop = population(size=100, ploidy=5, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), '5-ploid')
        self.assertRaises(exceptions.ValueError, population, size=[20, 20], ploidy=0)

    def testChromBeginEnd(self):
        'Testing genoStruTrait::ChromBegin(chrom), chromEnd(chrom)'
        pop = self.getPop()
        self.assertEqual(pop.chromBegin(0), 0)
        self.assertEqual(pop.chromBegin(1), 5)
        self.assertEqual(pop.chromEnd(0), 5)
        self.assertEqual(pop.chromEnd(1), 12)
        self.assertEqual(pop.numChrom(), 2)
        self.assertRaises(exceptions.IndexError, pop.chromBegin, 2 )
        self.assertRaises(exceptions.IndexError, pop.chromEnd, 2 )

    def testChromName(self):
        'Testing genoStruTrait::chromByName(name), chromName(chrom), chromNames()'
        pop = population(size=100, ploidy=2, loci=[5, 7])
        self.assertEqual(pop.chromName(0), 'chrom1')
        self.assertEqual(pop.chromName(1), 'chrom2')
        pop = population(size=100, ploidy=2, loci=[5, 7], chromNames=["c1", "c2"])
        self.assertEqual(pop.chromName(0), 'c1')
        self.assertEqual(pop.chromName(1), 'c2')
        self.assertEqual(pop.chromNames(), ('c1', 'c2'))
        self.assertEqual(pop.chromByName("c2"), 1)
        self.assertRaises(exceptions.ValueError, pop.chromByName, 'c3')
        self.assertRaises(exceptions.IndexError, pop.chromName, 2)

    def testChromType(self):
        'Testing genoStruTrait::chromType(chron), chromTypes()'
        pop = population(size=100, ploidy=2, loci=[2, 3, 2, 4],
        chromTypes=[Autosome, ChromosomeX, ChromosomeY, Customized])
        self.assertEqual(pop.chromType(0), Autosome)
        self.assertEqual(pop.chromType(1), ChromosomeX)
        self.assertEqual(pop.chromType(2), ChromosomeY)
        self.assertEqual(pop.chromType(3), Customized)
        self.assertRaises(exceptions.ValueError, population, ploidy=4,
            chromTypes=[Autosome, ChromosomeX, ChromosomeY, Customized])

    def testNumChrom(self):
        'Testing genoStruTrait::numChrom()'
        pop = self.getPop()
        self.assertEqual(pop.numChrom(), 2)

    def testAbsoluteLocusIndex(self):
        'Testing genoStruTrait::absLocusIndex(chrom, locus)'
        pop = self.getPop()
        self.assertEqual(pop.absLocusIndex(1, 5), 10)
        self.assertEqual(pop.locusPos(pop.absLocusIndex(1, 2) ), 6)
        self.assertRaises(exceptions.IndexError, pop.absLocusIndex, 2, 5 )

    def testChromLocusPair(self):
        'Testing genoStruTrait::chromLocusPair(locus)'
        pop = self.getPop()
        self.assertEqual(pop.chromLocusPair(10), (1, 5) )
        self.assertRaises(exceptions.IndexError, pop.chromLocusPair, 50 )

    def testLociName(self):
        'Testing genoStruTrait::lociByNames(names), lociNames(), locusByName(name), locusName(loc)'
        pop = self.getPop()
        self.assertEqual(pop.locusName(0), 'loc1-1')
        self.assertEqual(pop.locusName(1), 'loc1-2')
        self.assertEqual(pop.locusName(2), 'loc1-3')
        pop = population(loci=[1, 2], lociNames=['la', 'lb', 'lc'])
        self.assertEqual(pop.locusName(0), 'la')
        self.assertEqual(pop.locusName(1), 'lb')
        self.assertEqual(pop.locusName(2), 'lc')
        self.assertEqual(pop.locusByName("la"), 0)
        self.assertEqual(pop.locusByName("lb"), 1)
        self.assertEqual(pop.locusByName("lc"), 2)
        self.assertRaises(exceptions.IndexError, pop.locusName, 5)
        self.assertRaises(exceptions.ValueError, pop.locusByName, 'somename')
        self.assertRaises(exceptions.ValueError, pop.lociByNames, ['somename', 'other'])
        self.assertEqual(pop.lociByNames(['lb', 'lc']), (1, 2))
        self.assertEqual(pop.lociByNames(['lb', 'la']), (1, 0))

    def testLociDist(self):
        'Testing genoStruTrait::LociDist(loc1, loc2)'
        pop = self.getPop()
        self.assertEqual(pop.lociDist(0, 3), 3)
        self.assertEqual(pop.lociDist(2, 4), 2)
        self.assertRaises(exceptions.ValueError, pop.lociDist, 2, 8)

    def testLocusPos(self):
        'Testing genoStruTrait::lociPos(), locusPos(loc)'
        pop = self.getPop()
        self.assertEqual(pop.locusPos(10), 12)
        self.assertRaises(exceptions.IndexError, pop.locusPos, 20 )

    def testNumLoci(self):
        'Testing genoStruTrait::numLoci(chrom), numLoci(), totNumLoci()'
        pop = self.getPop()
        self.assertEqual(pop.numLoci(0), 5)
        self.assertEqual(pop.numLoci(1), 7)
        self.assertEqual(pop.totNumLoci(), 12)
        self.assertEqual(pop.numLoci(), (5, 7))
        self.assertRaises(exceptions.IndexError, pop.numLoci, 2 )

    def testAlleleName(self):
        'Testing genoStruTrait::AlleleName(allele), alleleNames()'
        pop = population(size=[20, 80], ploidy=2, loci=[5, 7])
        self.assertEqual(pop.alleleName(0), '0')
        self.assertEqual(pop.alleleName(1), '1')
        pop = self.getPop()
        self.assertEqual(pop.alleleName(0), '_')
        self.assertEqual(pop.alleleName(1), 'A')
        if AlleleType() != 'binary':
            self.assertEqual(pop.alleleName(2), 'C')
            self.assertEqual(pop.alleleName(3), 'T')
            self.assertEqual(pop.alleleName(4), 'G')
            self.assertEqual(pop.alleleNames(), ('_', 'A', 'C', 'T', 'G'))
        else:
            self.assertEqual(pop.alleleNames(), ('_', 'A'))
        self.assertRaises(exceptions.IndexError, pop.alleleName, MaxAllele()+1)

    def testInfoField(self):
        'Testing genoStruTrait::infoField(idx), infoFields(), infoIdx(name)'
        pop = population(10, infoFields=['age', 'fitness', 'trait1'])
        self.assertEqual(pop.infoField(0), 'age')
        self.assertEqual(pop.infoField(2), 'trait1')
        self.assertEqual(pop.infoIdx('age'), 0)
        self.assertEqual(pop.infoIdx('fitness'), 1)
        self.assertRaises(exceptions.IndexError, pop.infoField, 3)
        self.assertEqual(pop.infoFields(), ('age', 'fitness', 'trait1'))


if __name__ == '__main__':
    unittest.main()

