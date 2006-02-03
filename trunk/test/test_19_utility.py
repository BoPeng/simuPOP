#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of random number selector
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

class TestRNG(unittest.TestCase):

  def testSetRNG(self):
    'Testing all RNG types'
    for rg in listAllRNG():
      setRNG(rg)

  def testDefaultRNG(self):
    'Testing default RNG'
    rg = rng()
    self.assertEqual(rg.name(), 'mt19937')

  def testBinomial(self):
    'Testing binomial distribution'
    rg = rng()
    for n in range(1,10):
      rg.randBinomial(10, .7)
  
  def testUniform01(self):
    'Testing uniform distribution generator'
    rg = rng()
    for n in range(1,10):
      rg.randUniform01()
  
  def testLogging(self):
    'Testing logging output to a file'
    setLogOutput("session.log")
    # output stuff
    setLogOutput()
    os.remove('session.log')
  
if __name__ == '__main__':
  unittest.main()
