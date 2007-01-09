#!/usr/bin/env python
#
# Purpose:
#    testing the mpi version of simuPOP
#
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision: 475 $
# $LastChangedDate: 2006-10-04 01:00:38 -0500 (Wed, 04 Oct 2006) $
# 


import simuOpt
simuOpt.setOptions(quiet=True, optimized=False, mpi=True,
    chromMap=[2, 1, 3])

from simuPOP import *
import unittest, os, sys, exceptions, time

class TestMPI(unittest.TestCase):

    def testMPIStart(self):
        '''Start mpi'''
        if not mpi():
            return True
        testMPI()

    def testChromMap(self):
        '''Testing chromosome map'''
        if not mpi():
            return True
        # for all nodes
        pop = population(10, loci=[1,2,3,4,5,1], chromMap=[1,3,2])
        # chromosomes
        # 0 | 1 2 | 3 4 5 | 6 7 8 9 | 10 11 12 13 14 | 15
        # rank
        # chrom: 0 | 1 2 3 | 4 5
        # locus: 0 | 1 2 3 4 5 6 7 8 9 | 10 11 12 13 14 15
        self.assertEqual(pop.chromMap(), (1, 3, 2))
        self.assertEqual(pop.rankOfChrom(0), 1)
        self.assertEqual(pop.rankOfChrom(1), 2)
        self.assertEqual(pop.rankOfChrom(2), 2)
        self.assertEqual(pop.rankOfChrom(3), 2)
        self.assertEqual(pop.rankOfChrom(4), 3)
        self.assertEqual(pop.rankOfChrom(5), 3)
        self.assertRaises(exceptions.IndexError, pop.rankOfChrom, 6)
        # rank of locus
        self.assertEqual(pop.rankOfLocus(0), 1)
        self.assertEqual(pop.rankOfLocus(1), 2)
        self.assertEqual(pop.rankOfLocus(3), 2)
        self.assertEqual(pop.rankOfLocus(7), 2)
        self.assertEqual(pop.rankOfLocus(9), 2)
        self.assertEqual(pop.rankOfLocus(10), 3)
        self.assertEqual(pop.rankOfLocus(14), 3)
        self.assertEqual(pop.rankOfLocus(15), 3)
        self.assertRaises(exceptions.IndexError, pop.rankOfLocus, 16)
        # begin chrom of rank
        self.assertRaises(exceptions.IndexError, pop.beginChromOfRank, 0)
        self.assertEqual(pop.beginChromOfRank(1), 0)
        self.assertEqual(pop.beginChromOfRank(2), 1)
        self.assertEqual(pop.beginChromOfRank(3), 4)
        self.assertRaises(exceptions.IndexError, pop.beginChromOfRank, 4)
        # end chrom of rank
        self.assertRaises(exceptions.IndexError, pop.endChromOfRank, 0)
        self.assertEqual(pop.endChromOfRank(1), 1)
        self.assertEqual(pop.endChromOfRank(2), 4)
        self.assertEqual(pop.endChromOfRank(3), 6)
        self.assertRaises(exceptions.IndexError, pop.endChromOfRank, 4)
        # begin locus of rank
        self.assertRaises(exceptions.IndexError, pop.beginLocusOfRank, 0)
        self.assertEqual(pop.beginLocusOfRank(1), 0)
        self.assertEqual(pop.beginLocusOfRank(2), 1)
        self.assertEqual(pop.beginLocusOfRank(3), 10)
        self.assertRaises(exceptions.IndexError, pop.beginLocusOfRank, 4)
        # end chrom of rank
        self.assertRaises(exceptions.IndexError, pop.endLocusOfRank, 0)
        self.assertEqual(pop.endLocusOfRank(1), 1)
        self.assertEqual(pop.endLocusOfRank(2), 10)
        self.assertEqual(pop.endLocusOfRank(3), 16)
        self.assertRaises(exceptions.IndexError, pop.endLocusOfRank, 4)
        if mpiRank() == 0:
            self.assertRaises(exceptions.IndexError, pop.beginChrom)
            self.assertRaises(exceptions.IndexError, pop.endChrom)
            self.assertRaises(exceptions.IndexError, pop.beginLocus)
            self.assertRaises(exceptions.IndexError, pop.endLocus)
        elif mpiRank() > 0:
            self.assertEqual(pop.beginChrom(), pop.beginChromOfRank(mpiRank()))
            self.assertEqual(pop.endChrom(), pop.endChromOfRank(mpiRank()))
            self.assertEqual(pop.beginLocus(), pop.beginLocusOfRank(mpiRank()))
            self.assertEqual(pop.endLocus(), pop.endLocusOfRank(mpiRank()))
        #

    def testChromMapFromSetOptions(self):
        'Testing chromMap set by setOptions'
        # the map is 2, 1, 3
        # use system setOptions 
        pop = population(10, loci=[1,2,3,4,5,1])
        # chromosomes
        # 0 | 1 2 | 3 4 5 | 6 7 8 9 | 10 11 12 13 14 | 15
        # rank
        # chrom: 0 1 | 2 | 3 4 5
        # locus: 0 1 2 | 3 4 5 | 6 7 8 9 10 11 12 13 14 15
        self.assertEqual(pop.chromMap(), (2, 1, 3))
        self.assertEqual(pop.rankOfChrom(0), 1)
        self.assertEqual(pop.rankOfChrom(1), 1)
        self.assertEqual(pop.rankOfChrom(2), 2)
        self.assertEqual(pop.rankOfChrom(3), 3)
        self.assertEqual(pop.rankOfChrom(4), 3)
        self.assertEqual(pop.rankOfChrom(5), 3)
        self.assertRaises(exceptions.IndexError, pop.rankOfChrom, 6)
        # rank of locus
        self.assertEqual(pop.rankOfLocus(0), 1)
        self.assertEqual(pop.rankOfLocus(1), 1)
        self.assertEqual(pop.rankOfLocus(3), 2)
        self.assertEqual(pop.rankOfLocus(7), 3)
        self.assertEqual(pop.rankOfLocus(9), 3)
        self.assertEqual(pop.rankOfLocus(10), 3)
        self.assertEqual(pop.rankOfLocus(14), 3)
        self.assertEqual(pop.rankOfLocus(15), 3)
        self.assertRaises(exceptions.IndexError, pop.rankOfLocus, 16)
        # begin chrom of rank
        self.assertRaises(exceptions.IndexError, pop.beginChromOfRank, 0)
        self.assertEqual(pop.beginChromOfRank(1), 0)
        self.assertEqual(pop.beginChromOfRank(2), 2)
        self.assertEqual(pop.beginChromOfRank(3), 3)
        self.assertRaises(exceptions.IndexError, pop.beginChromOfRank, 4)
        # end chrom of rank
        self.assertRaises(exceptions.IndexError, pop.endChromOfRank, 0)
        self.assertEqual(pop.endChromOfRank(1), 2)
        self.assertEqual(pop.endChromOfRank(2), 3)
        self.assertEqual(pop.endChromOfRank(3), 6)
        self.assertRaises(exceptions.IndexError, pop.endChromOfRank, 4)
        # begin locus of rank
        self.assertRaises(exceptions.IndexError, pop.beginLocusOfRank, 0)
        self.assertEqual(pop.beginLocusOfRank(1), 0)
        self.assertEqual(pop.beginLocusOfRank(2), 3)
        self.assertEqual(pop.beginLocusOfRank(3), 6)
        self.assertRaises(exceptions.IndexError, pop.beginLocusOfRank, 4)
        # end chrom of rank
        self.assertRaises(exceptions.IndexError, pop.endLocusOfRank, 0)
        self.assertEqual(pop.endLocusOfRank(1), 3)
        self.assertEqual(pop.endLocusOfRank(2), 6)
        self.assertEqual(pop.endLocusOfRank(3), 16)
        self.assertRaises(exceptions.IndexError, pop.endLocusOfRank, 4)
        if mpiRank() == 0:
            self.assertRaises(exceptions.IndexError, pop.beginChrom)
            self.assertRaises(exceptions.IndexError, pop.endChrom)
            self.assertRaises(exceptions.IndexError, pop.beginLocus)
            self.assertRaises(exceptions.IndexError, pop.endLocus)
        elif mpiRank() > 0:
            self.assertEqual(pop.beginChrom(), pop.beginChromOfRank(mpiRank()))
            self.assertEqual(pop.endChrom(), pop.endChromOfRank(mpiRank()))
            self.assertEqual(pop.beginLocus(), pop.beginLocusOfRank(mpiRank()))
            self.assertEqual(pop.endLocus(), pop.endLocusOfRank(mpiRank()))



if __name__ == '__main__':
    unittest.main()
