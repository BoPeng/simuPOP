#!/usr/bin/env python
#
#  This is a unittest file for taggers
#
#  Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
# 

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

class  TestTagger(unittest.TestCase):
  
  def testParentsTagger(self):
    'Testing parents tagger.'
    simu = simulator( 
      population(size=20, ploidy=2, loci=[2,4], subPop=[5,15]),
      randomMating(numOffspring=2))
    simu.step([parentsTagger()])
    pop = simu.population(0)
    # check if all siblings have the same parents
    for sp in range(pop.numSubPop()):
      for i in range(pop.subPopSize(sp)/2):
        self.assertEqual( pop.individual(i*2,sp).tag(),
          pop.individual(i*2+1,sp).tag() )
      # note that the last one may be left alone

  def testInheritTagger(self):
    'Testing inherit tagger (pass info from parents to offspring'
    # this operator pass tag from one or both parents to offspring
    # the game is not:
    # who is the offspring of one parent?
    pop = population(size=20, ploidy=2, loci=[2,4], subPop=[5,15])
    pop.individual(0).setTag(1)
    pop.individual(5).setTag(2)
    simu = simulator( pop, randomMating())
    # other mode include TAG_Maternal, TAG_Both
    simu.step([inheritTagger(mode=TAG_Paternal)])
    # we only know subpopulation 0 can not have tag 2
    # we only know subpopulation 1 can not have tag 1
    for i in range(pop.subPopSize(0)):
      self.assertNotEqual( pop.individual(i,0).tag(), 2 )
    for i in range(pop.subPopSize(1)):
      self.assertNotEqual( pop.individual(i,1).tag(), 1 )
    # from this test, we can see that genetic drift 
    # can easily remove a signal (tag) from population.
    
    
if __name__ == '__main__':
  unittest.main()   
