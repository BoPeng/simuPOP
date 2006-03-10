#!/usr/bin/env python
#
# Purpose:
#   Testing quantitative trait.
# 
# Author:
#   Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 146 $
# $LastChangedDate: 2006-02-03 00:18:04 -0600 (Fri, 03 Feb 2006) $
# 

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions, math, random

class TestQuanTrait(unittest.TestCase):
  
  def setUp(self):
    self.pop = population(subPop=[5000], 
      ploidy=2, loci = [1])
    InitByValue(self.pop, 
      value = [[0,0],[0,1],[1,1]],
      indRange = [[0,1250], [1250,3750],[3750,5000]])
 
  def stdev(self, x):
    'calculate standard devisaion'
    mean = sum(x)/len(x)
    return math.sqrt( sum([(i-mean)**2 for i in x])/(len(x)-1) )
    
  def testMapQuanTrait(self):
    'Testing map qtrait'
    MapQuanTrait(self.pop, locus=0, 
      qtrait={'0-0':0, '0-1':1, '1-1':2})
    self.assertEqual(self.pop.dvars().qtrait, [0]*1250+[1]*2500+[2]*1250)
    # with variance
    MapQuanTrait(self.pop, locus=0, sigma=0.5,
      qtrait={'0-0':1, '0-1':1, '1-1':1})
    assert abs( sum(self.pop.dvars().qtrait)/5000 - 1) < 0.1
    # syandard deviation
    assert abs( self.stdev( self.pop.dvars().qtrait ) - 0.5) < 0.1
    
  def testMaQuanTrait(self):
    'Testing multi-allele qtrait'
    MaQuanTrait(self.pop, locus=0, wildtype=0, sigma=0,
      qtrait=[0, 1, 1])
    self.assertEqual(self.pop.dvars().qtrait, [0]*1250+[1]*2500+[1]*1250)
    # with variance
    MaQuanTrait(self.pop, locus=0, sigma=0.5, wildtype=0,
      qtrait=[0, .3, .5])
    assert abs( sum(self.pop.dvars().qtrait[:1250])/1250 - 0) < 0.1
    assert abs( sum(self.pop.dvars().qtrait[1250:3750])/2500 - 0.3) < 0.1
    assert abs( sum(self.pop.dvars().qtrait[3750:5000])/1250 - 0.5) < 0.1
    # syandard deviation
    assert abs( self.stdev( self.pop.dvars().qtrait[:1250] ) - 0.5) < 0.1
    assert abs( self.stdev( self.pop.dvars().qtrait[1250:3750] ) - 0.5) < 0.1
    assert abs( self.stdev( self.pop.dvars().qtrait[3750:5000] ) - 0.5) < 0.1
    # different sigma  
    MaQuanTrait(self.pop, locus=0, sigma=[0.5, 0.8, 1], wildtype=0,
      qtrait=[0, .3, .5])
    assert abs( sum(self.pop.dvars().qtrait[:1250])/1250 - 0) < 0.1
    assert abs( sum(self.pop.dvars().qtrait[1250:3750])/2500 - 0.3) < 0.1
    assert abs( sum(self.pop.dvars().qtrait[3750:5000])/1250 - 0.5) < 0.15
    # syandard deviation
    assert abs( self.stdev( self.pop.dvars().qtrait[:1250] ) - 0.5) < 0.1
    assert abs( self.stdev( self.pop.dvars().qtrait[1250:3750] ) - 0.8) < 0.1
    assert abs( self.stdev( self.pop.dvars().qtrait[3750:5000] ) - 1.0) < 0.1
    
  def testMultiLocusMaQuanTrait(self):
    'Test the multi-locus version of maQTrait'
    pop = population(1000, loci=[3,5])
    InitByFreq(pop, [.3, .7])
    #
    MaQuanTrait(pop, loci=[3,5], wildtype=0, sigma=0.5,
      qtrait=[0, .3, .5, 0.3, 0.6, 0.8, 0.1, 1, 0.8])

    
  def testMlQuanTrait(self):
    'Testing multi-locus qtrait'
    pop = population(1000, loci=[3,5])
    InitByFreq(pop, [.3, .7])
    #
    MlQuanTrait(pop, [
      maQuanTrait(locus=0,  wildtype=0, sigma=0.5,
        qtrait=[0, .3, .5]),
      mapQuanTrait(locus=1, sigma=1, 
        qtrait={'0-0':0, '0-1':1, '1-1':1})
      ],
      mode=PEN_Additive
    )
    #
    MlQuanTrait(pop, [
      maQuanTrait(locus=2,  wildtype=0, sigma=0.5, 
        qtrait=[0, .3, .5]),
      mapQuanTrait(locus=4, 
        qtrait={'0-0':0, '0-1':1, '1-1':1})
      ],
      mode=PEN_Multiplicative
    )

  def testPyQuanTrait(self):
    def qt(geno):
      if geno == [0, 0]:
        return random.normalvariate(0, 0.5)
      elif geno == [0, 1]:
        return random.normalvariate(0.5, 1)
      elif geno == [1, 0]:
        return random.normalvariate(0.5, 1)
      else:
        return random.normalvariate(1, 2)
    PyQuanTrait(self.pop, loci=[0], func=qt)
    assert abs( sum(self.pop.dvars().qtrait[:1250])/1250 - 0) < 0.1
    assert abs( sum(self.pop.dvars().qtrait[1250:3750])/2500 - 0.5) < 0.1
    assert abs( sum(self.pop.dvars().qtrait[3750:5000])/1250 - 1) < 0.3
    # syandard deviation
    assert abs( self.stdev( self.pop.dvars().qtrait[:1250] ) - 0.5) < 0.1
    assert abs( self.stdev( self.pop.dvars().qtrait[1250:3750] ) - 1) < 0.3
    assert abs( self.stdev( self.pop.dvars().qtrait[3750:5000] ) - 2) < 0.3
    #
    # multi-locus
    pop = population(1000, loci=[3,5])
    InitByFreq(pop, [.3, .7])
    def qt1(geno):
      assert len(geno) == 4
      return random.normalvariate(0, 0.5*sum(geno) )
    PyQuanTrait(pop, loci=[2,6], func=qt1)

if __name__ == '__main__':
  unittest.main()
