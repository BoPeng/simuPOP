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
  
  def assertMarkerProperties(self, m, pos, name, default, run_first):
    'Testing behavior to all marker types'
    self.assertEqual(m.position(), pos)
    m.position(0.5)
    self.assertEqual(m.position(), 0.5)
    self.assertEqual(m.type(), name)
    self.assertEqual(m.default_value(), default)
    # wrong position
    self.assertRaises(exceptions.ValueError,
      m.position, 2)
    # run_first (??)
    self.assertEqual(m.run_first(), run_first)
    
  def testSNPMarker(self):
    'Testing SNP markers'
    a = SNPMarker(0.11, 0.01, 0.90)
    # 
    self.assertMarkerProperties(a, pos=0.11, 
      name='snp', default=0, run_first=False)
    #
    self.assertEqual(a.low_freq(), 0.01)
    self.assertEqual(a.high_freq(), 0.90)
    # can not add_value
    self.assertRaises(exceptions.ValueError,
      a.add_value, 2)

  def testMicroSatelliteMarker(self):
    'Testing microsatellite markers'
    a = MicroSatelliteMarker(0.11, theta=1e-6, K=50)
    #
    self.assertMarkerProperties(a, pos=0.11, 
      name='ms', default=0, run_first=False)
    #
    self.assertEqual( a.theta(), 1e-6)
    self.assertEqual( a.K(), 50)
    # add value should be allowed FIXME: why add_value?

  def testTraitMarker(self):
    'Testing trait markers'
    a = TraitMarker(0.11, low_freq=0.01, high_freq=0.90)
    #
    self.assertMarkerProperties(a, pos=0.11, 
      name='trait', default=0, run_first=True)
    #
    self.assertEqual(a.low_freq(), 0.01)
    self.assertEqual(a.high_freq(), 0.90)
    # can not add value 
    self.assertRaises(exceptions.ValueError,
      a.add_value, 2)

  def testDescender(self):
    'Testing descender'
    #dec = Descender(conf)
    #dec.evolve()

  def assertEpochProperties(self, epo, st, et):
    'Assert properties for all epochs'
    self.assertEqual(epo.start_time(), st)
    self.assertEqual(epo.end_time(), et)
    print epo.earliest_event()
    
  def testMigrationEvent(self):
    'Testing migration events (epochs.hh, configuration.hh)'
    a = Migration(source=0, destination=1,
      migration_rate=0.001, start_time=0,
      end_time= 0.5)
    #
    self.assertEpochProperties(a, st=0, et=0.5)
    self.assertEqual(a.source(), 0)
    self.assertEqual(a.destination(), 1)
    self.assertEqual(a.migration_rate(), 0.001)
  
  def testBottleNeckEpoch(self):
    'Testing BottleNeckEpoch (epochs.hh)'
    a = BottleNeckEpoch(population=0,  
      scale_fraction=1.2, 
      start_time=0, end_time= 0.5)
    #
    self.assertEpochProperties(a, st=0, et=0.5)
    self.assertEqual( a.scale_fraction(), 1.2)
    
  def testGrowthEpoch(self):
    'Testing GrowthEpoch (epochs.hh)'
    a = GrowthEpoch(population=0,  
      beta=1.2, start_time=0, end_time= 0.5)
    #
    self.assertEpochProperties(a, st=0, et=0.5)
    self.assertEqual( a.beta(), 1.2)    
   
  def testScheduler(self):
    'Testing scheduler (configuration.hh)'
    #a = Scheduler()
    #a.add_event(e)
    #a.remove_event(e)
    #print next_event
    
  def testConfiguration(self):
    'Testing configure object (configuration.hh)'
    conf = Configuration(
      popSizes=[1000,2000,3000], 
      markers = [
        SNPMarker(0.1, 0.1, 0.9), 
        SNPMarker(0.2, 0.1, 0.9)
        ],
      epochs = [ GrowthEpoch(0, 1.2, 0, 0.5)], 
      rho=0.1, Q=0.01, gamma=0.001, growth=0.0001)
    # 
    self.assertEqual( conf.pop_sizes(), (1,2,3) )
    # FIXME: not sure why no_leaves is 6
    self.assertEqual( conf.no_leaves(), 6)
    self.assertEqual( conf.no_markers(), 2)
    self.assertEqual( conf.position(0), 0.1)
    self.assertEqual( conf.position(1), 0.2)
    self.assertRaises( exceptions.IndexError,
      conf.position, 2)
    self.assertEqual( conf.rho(), 0.1)
    self.assertEqual( conf.Q(), 0.01)
    self.assertEqual( conf.gamma(), 0.001)
    self.assertEqual( conf.growth(), 0.0001)
    #
    self.assertEqual( conf.marker(0).type(), 'snp')
    self.assertRaises( exceptions.IndexError,
      conf.marker, 2)
    # FIXME: what are first_marker and plain_marker
    assert not conf.is_first_marker(0)
    assert not conf.is_first_marker(1)
    assert conf.is_plain_marker(0)
    assert conf.is_plain_marker(1)
    self.assertRaises( exceptions.IndexError,
      conf.first_marker, 0)
    self.assertEqual(conf.plain_marker(0).position(), 0.1)
    self.assertEqual(conf.plain_marker(1).position(), 0.2)

  def testDistributionFunctions(self):
    'Testing distribution functions'
    # Should we simply trust these functions?
    expdev(2.)
    expdev(2, .2)
    expdist(.2, 3.)
    uniform()
    uniform(1,2)
    # what is this guy?
    uniform(1,2,3)  
    random_sign()
    irand(10)
    # can not test two_int_rand


  def testPopulation(self):
    'Testing population'
    # need ARG...

if __name__ == '__main__':
  unittest.main()
