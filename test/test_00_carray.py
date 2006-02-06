#!/usr/bin/env python
#
# Purpose:
#
# This is a unittest file for population object
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

  def testFloatCarray(self):
    'Testing float carray type returned by arrLociPos'
    pop = population(loci=[3,4], lociPos=[1,2,3,4,5,6,7])
    arr = pop.arrLociPos()
    # can print
    # print arr, ignore
    # expression
    self.assertEqual(str(arr), "[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]")
    # count
    self.assertEqual(arr.count(2), 1)
    arr[1] = 1
    self.assertEqual(arr.count(1.0), 2)
    self.assertEqual(arr.count(2), 0)
    # index
    self.assertRaises(exceptions.ValueError,
      arr.index, 2)
    self.assertEqual(arr.index(1), 0)
    self.assertEqual(arr.index(6), 5)
    # can read write
    arr[0] = 0.5
    arr[3] = 3.5
    self.assertEqual(arr[0], 0.5)
    self.assertEqual(arr[3], 3.5)
    # convert to list
    self.assertEqual( arr.tolist(), [0.5, 1.0, 3.0, 3.5, 5.0, 6.0, 7.0])
    # direct comparison
    self.assertEqual( arr == [0.5, 1.0, 3.0, 3.5, 5.0, 6.0, 7.0], True)
    # convert to list
    self.assertEqual( list(arr), [0.5, 1.0, 3.0, 3.5, 5.0, 6.0, 7.0])
    # slice
    arr[:] = [1,2,3,4,5,6,7]
    arr1 = arr[:3]
    self.assertEqual( arr1, [1, 2, 3])
    arr1 = arr[3:5]
    self.assertEqual( arr1, [4, 5])
    # assign slice
    arr1[:] = 10
    self.assertEqual( arr1, [10, 10] )
    # IMPORTANT NOTE that arr will also be affected
    self.assertEqual( arr, [1,2,3,10,10,6,7])
    # assign vector
    arr1[:] = [30,40]
    self.assertEqual( arr1, [30, 40] )
    self.assertEqual( arr, [1,2,3,30,40,6,7])
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
    arr[1:3] = arr[3:5]
    self.assertEqual( arr, [1, 30,40, 30,40,6,7])
    
  def testGenotypeCarray(self):
    'Testing allele carray type returned by arrGenotype'
    pop = population(size=2, loci=[1,2])
    arr = pop.arrGenotype()
    arr[:] = [0,1,2]*4
    # can print
    # print arr
    # expression
    if alleleType() != 'binary':
      self.assertEqual( str(arr), "[0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]")
    else:
      self.assertEqual( str(arr), "[0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]")
    # count
    if alleleType() != 'binary':
      self.assertEqual(arr[5], 2)
      self.assertEqual(arr.count(2), 4)
    else:
      self.assertEqual(arr[5], 1)
      self.assertEqual(arr.count(2), 0)
      self.assertEqual(arr.count(1), 8)
    arr[0] = 1
    if alleleType() != 'binary':
      self.assertEqual(arr.count(1), 5)
    else:
      self.assertEqual(arr.count(1), 9)
    # index
    self.assertRaises(exceptions.ValueError,
      arr.index, 10)
    self.assertEqual(arr.index(1), 0)
    # can read write
    arr[3] = 3
    if alleleType() != 'binary':
      self.assertEqual(arr[3], 3)
    else:
      self.assertEqual(arr[3], 1)
    # convert to list
    arr[:] = [0,1,2]*4
    if alleleType() != 'binary':
      self.assertEqual( arr.tolist(), [0,1,2]*4)
    else:
      self.assertEqual( arr.tolist(), [0,1,1]*4)
      self.assertNotEqual( arr.tolist(), [0,1,2]*4)
    # convert to list
    if alleleType() != 'binary':
      self.assertEqual( arr, [0,1,2]*4)
    else:
      self.assertEqual( arr, [0,1,1]*4)
    # slice
    arr[:] = [0,1,2]*4
    arr1 = arr[:3]
    if alleleType() != 'binary':
      self.assertEqual( arr1, [0, 1, 2])
    else:
      self.assertEqual( arr1, [0, 1, 1])
    arr1 = arr[3:5]
    if alleleType() != 'binary':
      self.assertEqual( arr1, [0, 1])
    else:
      self.assertEqual( arr1, [0, 1])
    # assign slice
    arr1[:] = 5
    # IMPORTANT NOTE that arr will also be affected
    if alleleType() != 'binary':
      self.assertEqual( arr1, [5, 5] )
      self.assertEqual( arr, [0,1,2,5,5,2,0,1,2,0,1,2])
    else:
      self.assertEqual( arr1, [1, 1] )
      self.assertEqual( arr, [0,1,1,1,1,1,0,1,1,0,1,1])
    # assign vector
    arr1[:] = [0,0]
    self.assertEqual( arr1, [0, 0] )
    if alleleType() != 'binary':
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
    if alleleType() != 'binary':
      self.assertEqual( arr, [0,1,2,0,1,2,0,1,2,0,1,2])
    else:
      self.assertEqual( arr, [0,1,1,0,1,1,0,1,1,0,1,1])

if __name__ == '__main__':
  unittest.main()
