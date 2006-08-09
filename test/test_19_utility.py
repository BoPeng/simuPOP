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
    for rg in ListAllRNG():
      SetRNG(rg)

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

	def testBernulliTrials(self):
		'Testing bernullitrials'
		rg = rng()
		p = [0.00001, 0.001, 0.5, 0.99]
		N = 1000000
		bt = BernulliTrials(rg, p, N)
		bt.doTrial()
		for i in range(len(p)):
			prop = bt.succProp(i)
			# binomial, mean p, variance = p(1-p)/n
			std = math.sqrt(p[i]*(1.-p[i])/N)
			assert  prop > p - 4*std and prop < p - 4*std

	def testSeed(self):
		'repeated set rng() without seed, and see if the seed repeat'
		seed = []
		for i in 1000:
			SetRNG(ListAllRNG()[0])
			sd = rng.seed()
			assert sd not in seed
			seed.append(sd)
		# test set seed
		for i in 1000:
			SetRNG(ListAllRNG()[0], i)
			self.assertEqual(rng().seed(), i)
			
  
if __name__ == '__main__':
  unittest.main()
