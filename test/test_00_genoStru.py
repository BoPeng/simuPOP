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

class TestCarray(unittest.TestCase):
    
    def testGenotypeCarray(self):
        'Testing allele carray type returned by genotype'
        pop = population(size=2, loci=[2,1])
        InitByValue(pop, [1,2,3])
        arr = pop.genotype(True)
        arr[:] = [0,1,2]*4
        # can print
        # expression
        if AlleleType() != 'binary':
            self.assertEqual( str(arr), "[0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]")
        else:
            self.assertEqual( str(arr), "[0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]")
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
        arr[:] = [0,1,2]*4
        if AlleleType() != 'binary':
            self.assertEqual( arr.tolist(), [0,1,2]*4)
        else:
            self.assertEqual( arr.tolist(), [0,1,1]*4)
            self.assertNotEqual( arr.tolist(), [0,1,2]*4)
        # convert to list
        if AlleleType() != 'binary':
            self.assertEqual( arr, [0,1,2]*4)
        else:
            self.assertEqual( arr, [0,1,1]*4)
        # slice
        arr[:] = [0,1,2]*4
        arr1 = arr[:3]
        if AlleleType() != 'binary':
            self.assertEqual( arr1, [0, 1, 2])
        else:
            self.assertEqual( arr1, [0, 1, 1])
        arr1 = arr[3:5]
        if AlleleType() != 'binary':
            self.assertEqual( arr1, [0, 1])
        else:
            self.assertEqual( arr1, [0, 1])
        # assign slice
        arr1[:] = 5
        # IMPORTANT NOTE that arr will also be affected
        if AlleleType() != 'binary':
            self.assertEqual( arr1, [5, 5] )
            self.assertEqual( arr, [0,1,2,5,5,2,0,1,2,0,1,2])
        else:
            self.assertEqual( arr1, [1, 1] )
            self.assertEqual( arr, [0,1,1,1,1,1,0,1,1,0,1,1])
        # assign vector
        arr1[:] = [0,0]
        self.assertEqual( arr1, [0, 0] )
        if AlleleType() != 'binary':
            self.assertEqual( arr, [0,1,2,0,0,2,0,1,2,0,1,2])
        else:
            self.assertEqual( arr, [0,1,1,0,0,1,0,1,1,0,1,1])
        # assign vector of different length
        try:
            arr1[:] = [50, 60, 70]
        except exceptions.ValueError:
            pass
        try:
            arr1[:] = arr[1:2]
        except exceptions.ValueError:
            pass
        # assign from another part
        arr[:6] = arr[6:12]
        if AlleleType() != 'binary':
            self.assertEqual( arr, [0,1,2,0,1,2,0,1,2,0,1,2])
        else:
            self.assertEqual( arr, [0,1,1,0,1,1,0,1,1,0,1,1])

    def testBlah():
        ''
        pass

if __name__ == '__main__':
    unittest.main()
