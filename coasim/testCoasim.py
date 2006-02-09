#!/usr/bin/env python
#
# Purpose:
#   Test the python binding of coaSim (http://www.birc.dk/Software/CoaSim/)
#
# Author:
#   Bo Peng
#
# $LastChangedRevision: 154 $
# $LastChangedDate: 2006-02-06 01:05:52 -0600 (Mon, 06 Feb 2006) $
# 

from coaSim import *
import unittest, os, sys, exceptions

class TestCoaSim(unittest.TestCase):
  
  def testSNPMarker(self):
    'Testing SNP markers'
    a = SNPMarker(0.11, 0.01, 0.90)
    self.assertEqual(a.position(), 0.11)
    self.assertEqual(a.type(), 'snp')
    # set position
    a.position(0.5)
    self.assertEqual(a.position(), 0.5)
    # wrong position, (strange, exception does not work right now)
    # a.position(2)
    self.assertEqual(a.default_value(), 0)
    self.assertEqual(a.low_freq(), 0.01)
    self.assertEqual(a.high_freq(), 0.90)
    # a.add_value(2)

  def testMicroMarker(self):
    'Testing microsatellite markers'
    a = MicroSatelliteMarker(0.11, theta=1e-6, K=50)
    self.assertEqual(a.position(), 0.11)
    self.assertEqual(a.type(), 'ms')
    # set position
    a.position(0.5)
    self.assertEqual(a.position(), 0.5)
    # wrong position, (strange, exception does not work right now)
    # a.position(2)
    self.assertEqual(a.default_value(), 0)

  def testTraitMarker(self):
    'Testing trait markers'
    a = TraitMarker(0.11, 0.01, 0.90)
    self.assertEqual(a.position(), 0.11)
    self.assertEqual(a.type(), 'trait')
    # set position
    a.position(0.5)
    self.assertEqual(a.position(), 0.5)
    # wrong position, (strange, exception does not work right now)
    # a.position(2)
    self.assertEqual(a.default_value(), 0)
    self.assertEqual(a.low_freq(), 0.01)
    self.assertEqual(a.high_freq(), 0.90)
    # a.add_value(2)

  def testConfiguration(self):
    'Testing configure object'
    no_leaves = 0
    conf = Configuration(
      popSizes=[1,2,3], 
      markers = [SNPMarker(0.1, 0.1, 0.9), SNPMarker(0.1, 0.1, 0.9)],
      epochs = [], 
      rho=0.001, Q=0.001, gamma=0.001, growth=0.001)



if __name__ == '__main__':
  unittest.main()
