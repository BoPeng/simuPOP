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

class TestGenoStru(unittest.TestCase):
    # define a function to create basic populations
    def getPop(self):
        pop = Population(size=[20, 80], ploidy=2, loci=[5, 7],
            lociPos=[2, 3, 4, 5, 6, 2, 4, 6, 8, 10, 12, 14],
            alleleNames=['_', 'A', 'C', 'T', 'G'],
            infoFields=['a', 'b'])
        initSex(pop)
        return pop

    def testLociPos(self):
        'Testing the specification of loci position'
        pop = Population(loci=[4, 5], lociPos=[1, 3, 2, 4, 1, 5, 4, 3, 2])
        self.assertEqual(pop.lociPos(), (1, 2, 3, 4, 1, 2, 3, 4, 5))
        self.assertRaises(ValueError, Population, loci=[4, 5], lociPos=[1, 3, 2, 4, 1, 3, 4, 3, 2])
        #
        pop = Population(loci=[4, 2], lociPos=[1, 3, 2, 4, 10, 5], lociNames=['1', '2', '3', '4', '5', '6'])
        self.assertEqual(pop.lociPos(), (1, 2, 3, 4, 5, 10))
        self.assertEqual(pop.lociNames(), ('1', '3', '2', '4', '6', '5'))
        #
        pop = Population(loci=[4, 2], lociPos=[1, 3, 2, 4, 10, 5],
            alleleNames=[['1', 'x'], ['2', 'a'], ['3', 'y', 'z'], ['4', 'b'], ['5', 'd'], ['6', 'c']])
        self.assertEqual(pop.lociPos(), (1, 2, 3, 4, 5, 10))
        self.assertEqual(pop.alleleNames(0), ('1', 'x'))
        if moduleInfo()['alleleType'] == 'binary':
            self.assertEqual(pop.alleleNames(1), ('3', 'y'))
        else:
            self.assertEqual(pop.alleleNames(1), ('3', 'y', 'z'))
        self.assertEqual(pop.alleleNames(2), ('2', 'a'))
        self.assertEqual(pop.alleleNames(3), ('4', 'b'))
        self.assertEqual(pop.alleleNames(4), ('6', 'c'))
        self.assertEqual(pop.alleleNames(5), ('5', 'd'))

    def testGenotypeCarray(self):
        'Testing allele carray type returned by genotype'
        pop = Population(size=2, loci=[2, 1])
        initGenotype(pop, genotype=[1, 2, 3])
        arr = pop.genotype()
        arr[:] = [0, 1, 2]*4
        # can print
        # expression
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(str(arr), "[0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]")
        else:
            self.assertEqual(str(arr), "[0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]")
        # count
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr[5], 2)
            self.assertEqual(arr.count(2), 4)
        else:
            self.assertEqual(arr[5], 1)
            self.assertEqual(arr.count(2), 0)
            self.assertEqual(arr.count(1), 8)
        arr[0] = 1
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr.count(1), 5)
        else:
            self.assertEqual(arr.count(1), 9)
        # index
        self.assertRaises(ValueError,
            arr.index, 10)
        self.assertEqual(arr.index(1), 0)
        # can read write
        arr[3] = 3
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr[3], 3)
        else:
            self.assertEqual(arr[3], 1)
        # convert to list
        arr[:] = [0, 1, 2]*4
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr.tolist(), [0, 1, 2]*4)
        else:
            self.assertEqual(arr.tolist(), [0, 1, 1]*4)
            self.assertNotEqual(arr.tolist(), [0, 1, 2]*4)
        # convert to list
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr, [0, 1, 2]*4)
        else:
            self.assertEqual(arr, [0, 1, 1]*4)
        # slice
        arr[:] = [0, 1, 2]*4
        arr1 = arr[:3]
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr1, [0, 1, 2])
        else:
            self.assertEqual(arr1, [0, 1, 1])
        arr1 = arr[3:5]
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr1, [0, 1])
        else:
            self.assertEqual(arr1, [0, 1])
        # assign slice
        arr1[:] = 5
        # IMPORTANT NOTE that arr will also be affected
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr1, [5, 5] )
            self.assertEqual(arr, [0, 1, 2, 5, 5, 2, 0, 1, 2, 0, 1, 2])
        else:
            self.assertEqual(arr1, [1, 1] )
            self.assertEqual(arr, [0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1])
        # assign vector
        arr1[:] = [0, 0]
        self.assertEqual(arr1, [0, 0] )
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr, [0, 1, 2, 0, 0, 2, 0, 1, 2, 0, 1, 2])
        else:
            self.assertEqual(arr, [0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1])
        # assign from another part
        arr[:6] = arr[6:12]
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(arr, [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2])
        else:
            self.assertEqual(arr, [0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1])

    def testPloidy(self):
        'Testing genoStruTrait::Ploidy(), PloidyName()'
        pop = self.getPop()
        self.assertEqual(pop.ploidy(), 2)
        self.assertEqual(pop.ploidyName(), 'diploid')
        pop = Population(size=100, ploidy=HAPLODIPLOID, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), 'haplodiploid')
        pop = Population(size=100, ploidy=1, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), 'haploid')
        pop = Population(size=100, ploidy=3, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), 'triploid')
        pop = Population(size=100, ploidy=4, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), 'tetraploid')
        pop = Population(size=100, ploidy=5, loci=[5, 7])
        self.assertEqual(pop.ploidyName(), '5-ploid')
        self.assertRaises(ValueError, Population, size=[20, 20], ploidy=0)

    def testChromBeginEnd(self):
        'Testing genoStruTrait::ChromBegin(chrom), chromEnd(chrom)'
        pop = self.getPop()
        self.assertEqual(pop.chromBegin(0), 0)
        self.assertEqual(pop.chromBegin(1), 5)
        self.assertEqual(pop.chromEnd(0), 5)
        self.assertEqual(pop.chromEnd(1), 12)
        self.assertEqual(pop.numChrom(), 2)
        self.assertRaises(IndexError, pop.chromBegin, 2 )
        self.assertRaises(IndexError, pop.chromEnd, 2 )

    def testChromName(self):
        'Testing genoStruTrait::chromByName(name), chromName(chrom), chromNames()'
        pop = Population(size=100, ploidy=2, loci=[5, 7])
        self.assertEqual(pop.chromName(0), '')
        self.assertEqual(pop.chromName(1), '')
        pop = Population(size=100, ploidy=2, loci=[5, 7], chromNames=["c1", "c2"])
        self.assertEqual(pop.chromName(0), 'c1')
        self.assertEqual(pop.chromName(1), 'c2')
        self.assertEqual(pop.chromNames(), ('c1', 'c2'))
        self.assertEqual(pop.chromByName("c2"), 1)
        self.assertRaises(ValueError, pop.chromByName, 'c3')
        self.assertRaises(IndexError, pop.chromName, 2)

    def testChromType(self):
        'Testing genoStruTrait::chromType(chron), chromTypes()'
        pop = Population(size=100, ploidy=2, loci=[2, 3, 2, 4, 4],
            chromTypes=[AUTOSOME, CHROMOSOME_X, CHROMOSOME_Y, MITOCHONDRIAL, CUSTOMIZED])
        self.assertEqual(pop.chromType(0), AUTOSOME)
        self.assertEqual(pop.chromType(1), CHROMOSOME_X)
        self.assertEqual(pop.chromType(2), CHROMOSOME_Y)
        self.assertEqual(pop.chromType(3), MITOCHONDRIAL)
        self.assertEqual(pop.chromType(4), CUSTOMIZED)
        self.assertRaises(ValueError, Population, ploidy=4,
            chromTypes=[AUTOSOME, CHROMOSOME_X, CHROMOSOME_Y, CUSTOMIZED])
        #
        self.assertRaises(ValueError, Population, loci=[3, 4],
            chromTypes=[MITOCHONDRIAL]*2)

    def testNumChrom(self):
        'Testing genoStruTrait::numChrom()'
        pop = self.getPop()
        self.assertEqual(pop.numChrom(), 2)

    def testAbsoluteLocusIndex(self):
        'Testing genoStruTrait::absLocusIndex(chrom, locus)'
        pop = self.getPop()
        self.assertEqual(pop.absLocusIndex(1, 5), 10)
        self.assertEqual(pop.locusPos(pop.absLocusIndex(1, 2) ), 6)
        self.assertRaises(IndexError, pop.absLocusIndex, 2, 5 )

    def testChromLocusPair(self):
        'Testing genoStruTrait::chromLocusPair(locus)'
        pop = self.getPop()
        self.assertEqual(pop.chromLocusPair(10), (1, 5) )
        self.assertRaises(IndexError, pop.chromLocusPair, 50 )

    def testLociName(self):
        'Testing genoStruTrait::lociByNames(names), lociNames(), locusByName(name), locusName(loc)'
        pop = self.getPop()
        self.assertEqual(pop.locusName(0), '')
        self.assertEqual(pop.locusName(1), '')
        self.assertEqual(pop.locusName(2), '')
        pop = Population(loci=[1, 2], lociNames=['la', 'lb', 'lc'])
        self.assertEqual(pop.locusName(0), 'la')
        self.assertEqual(pop.locusName(1), 'lb')
        self.assertEqual(pop.locusName(2), 'lc')
        self.assertEqual(pop.locusByName("la"), 0)
        self.assertEqual(pop.locusByName("lb"), 1)
        self.assertEqual(pop.locusByName("lc"), 2)
        self.assertRaises(IndexError, pop.locusName, 5)
        self.assertRaises(ValueError, pop.locusByName, 'somename')
        self.assertRaises(ValueError, pop.lociByNames, ['somename', 'other'])
        self.assertEqual(pop.lociByNames(['lb', 'lc']), (1, 2))
        self.assertEqual(pop.lociByNames(['lb', 'la']), (1, 0))

    def testIndexesOfLoci(self):
        'Testing genoStruTrait::indexesOfLoci'
        pop = Population(size=10, loci=10)
        self.assertEqual(pop.indexesOfLoci([('', 1), ('', 3)]), (0, 2))
        self.assertRaises(ValueError, pop.indexesOfLoci, [('', 1), ('', 3.5)])
        pop = Population(size=10, loci=[5,5])
        self.assertEqual(pop.indexesOfLoci(loci=[('', 1), ('', 3)]), (5, 7))
        #
        pop = Population(size=10, loci=[5,5], chromNames=['a', 'b'], lociPos=[x*0.1 for x in range(10)])
        self.assertEqual(pop.indexesOfLoci([('a', 0.3), ('b', 0.9)]), (3, 9))
        self.assertRaises(ValueError, pop.indexesOfLoci, loci=[('a', .1), ('c', .3)])
        # single number?
        pop = Population(size=10, loci=[5,5], chromNames=['a', 'b'], lociPos=[x*0.1 for x in range(10)])
        self.assertEqual(pop.indexesOfLoci(loci=('a', 0.3)), (3,))
        self.assertRaises(ValueError, pop.indexesOfLoci, loci=('a', 1.9))

    def testLociDist(self):
        'Testing genoStruTrait::LociDist(loc1, loc2)'
        pop = self.getPop()
        self.assertEqual(pop.lociDist(0, 3), 3)
        self.assertEqual(pop.lociDist(2, 4), 2)
        self.assertRaises(ValueError, pop.lociDist, 2, 8)

    def testLocusPos(self):
        'Testing genoStruTrait::lociPos(), locusPos(loc)'
        pop = self.getPop()
        self.assertEqual(pop.locusPos(10), 12)
        self.assertRaises(IndexError, pop.locusPos, 20 )

    def testNumLoci(self):
        'Testing genoStruTrait::numLoci(chrom), numLoci(), totNumLoci()'
        pop = self.getPop()
        self.assertEqual(pop.numLoci(0), 5)
        self.assertEqual(pop.numLoci(1), 7)
        self.assertEqual(pop.totNumLoci(), 12)
        self.assertEqual(pop.numLoci(), (5, 7))
        self.assertRaises(IndexError, pop.numLoci, 2 )

    def testAlleleName(self):
        'Testing genoStruTrait::AlleleName(allele), alleleNames()'
        pop = Population(size=[20, 80], ploidy=2, loci=[5, 7])
        self.assertEqual(pop.alleleName(0), '0')
        self.assertEqual(pop.alleleName(1), '1')
        pop = self.getPop()
        self.assertEqual(pop.alleleName(0), '_')
        self.assertEqual(pop.alleleName(1), 'A')
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(pop.alleleName(2), 'C')
            self.assertEqual(pop.alleleName(3), 'T')
            self.assertEqual(pop.alleleName(4), 'G')
            self.assertEqual(pop.alleleNames(), ('_', 'A', 'C', 'T', 'G'))
        else:
            self.assertEqual(pop.alleleNames(), ('_', 'A'))
        self.assertRaises((IndexError, TypeError, OverflowError), pop.alleleName, 
            moduleInfo()['maxAllele']+1)
        # test locus-specific allele names
        pop = Population(size=[20, 80], ploidy=2, loci=[1, 2],
            alleleNames = [['A1', 'A2'], ['B1', 'B2', 'B3'], ['C1', 'C2', 'C3', 'C4']])
        self.assertEqual(pop.alleleName(0), 'A1')
        self.assertEqual(pop.alleleName(1), 'A2')
        self.assertEqual(pop.alleleName(0, 1), 'B1')
        self.assertEqual(pop.alleleName(1, 2), 'C2')
        if moduleInfo()['alleleType'] != 'binary':
            self.assertEqual(pop.alleleName(2), '2')
            self.assertEqual(pop.alleleNames(1), ('B1', 'B2', 'B3'))
        else:
            self.assertEqual(pop.alleleName(1), 'A2')
            self.assertEqual(pop.alleleNames(1), ('B1', 'B2'))
        self.assertRaises((IndexError, TypeError, OverflowError), pop.alleleName, 
            moduleInfo()['maxAllele']+1)


    def testInfoField(self):
        'Testing genoStruTrait::infoField(idx), infoFields(), infoIdx(name)'
        pop = Population(10, infoFields=['age', 'fitness', 'trait1'])
        self.assertEqual(pop.infoField(0), 'age')
        self.assertEqual(pop.infoField(2), 'trait1')
        self.assertEqual(pop.infoIdx('age'), 0)
        self.assertEqual(pop.infoIdx('fitness'), 1)
        self.assertRaises(IndexError, pop.infoField, 3)
        self.assertEqual(pop.infoFields(), ('age', 'fitness', 'trait1'))


if __name__ == '__main__':
    unittest.main()

