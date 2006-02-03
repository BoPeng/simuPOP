#!/usr/bin/env python
#
# This is a unittest file for migration operators
#
# Bo Peng (bpeng@rice.edu)
# 
# Last Modified: Sep, 2005
# 

from simuPOP import *
from simuRPy import *
import unittest, time, exceptions

class TestMigration(unittest.TestCase):
  
  def testMigrateByCounts(self):
    pop = population(subPop=[2000,4000,4000], loci=[2])
    Migrate(pop, mode=MigrByCounts, 
      rate = [ [0, 50, 50],
               [0, 0, 0],
               [0, 0, 0] ])
    assert pop.subPopSizes() == (1900, 4050, 4050)
    Migrate(pop, mode=MigrByCounts, 
      rate = [ [0, 50, 50],
               [50, 0, 0],
               [50, 0, 0] ])
    assert pop.subPopSizes() == (1900, 4050, 4050)
    # rate should be a matrix
    self.assertRaises(exceptions.ValueError, 
      Migrate, pop, [ [0, 50],
               [50, 0, 0],
               [50, 0, 0] ], MigrByCounts)
    
  def testMigrateByProportion(self):
    pop = population(subPop=[2000,4000,4000], loci=[2])
    # now if we want to inject a mutation whenever fixation happens
    Migrate(pop, mode=MigrByProportion, 
      rate = [ [0, .05, .05],
               [0.025, 0, 0],
               [0.025, 0, 0] ])
    assert pop.subPopSizes() == (2000, 4000, 4000)
    Migrate(pop, mode=MigrByProportion, 
      rate = [ [0, .25, .25],
               [0.25, 0, 0],
               [0, 0.25, 0] ])
    assert pop.subPopSizes() == (2000, 4500, 3500)
    
  def testMigrateByProbability(self):
    pop = population(subPop=[2000,4000,4000], loci=[2])
    # now if we want to inject a mutation whenever fixation happens
    Migrate(pop, mode=MigrByProbability, 
      rate = [ [0, .05, .05],
               [0.025, 0, 0],
               [0.025, 0, 0] ])
    # print pop.subPopSizes()
    assert abs(pop.subPopSize(0) - 2000) < 100
    assert abs(pop.subPopSize(1) - 4000) < 100
    assert abs(pop.subPopSize(2) - 4000) < 100
    Migrate(pop, mode=MigrByProbability, 
      rate = [ [0, .25, .25],
               [0.25, 0, 0],
               [0, 0.25, 0] ])
    # print pop.subPopSizes()
    assert abs(pop.subPopSize(0) - 2000) < 100
    assert abs(pop.subPopSize(1) - 4500) < 100
    assert abs(pop.subPopSize(2) - 3500) < 100

  def testMigrateFromTo(self):
    pop = population(subPop=[2000,4000,4000], loci=[2])
    # now if we want to inject a mutation whenever fixation happens
    Migrate(pop, mode=MigrByProbability, 
      fromSubPop = [0], toSubPop = [1,2], 
      rate = [.05, .05] )
    # print pop.subPopSizes()
    assert abs(pop.subPopSize(0) - 1800) < 50
    assert abs(pop.subPopSize(1) - 4100) < 50
    assert abs(pop.subPopSize(2) - 4100) < 50
    # other parameter form can be used as well
    pop = population(subPop=[2000,4000,4000], loci=[2])
    Migrate(pop, mode=MigrByProbability, 
      fromSubPop = 0, toSubPop = [1,2], 
      rate = [[.05, .05]] )
    assert abs(pop.subPopSize(0) - 1800) < 50
    assert abs(pop.subPopSize(1) - 4100) < 50
    assert abs(pop.subPopSize(2) - 4100) < 50
  
  def testMigrConstAlleleFreq(self):
    ' test if allele frequency changes '
    pop = population(subPop=[2000,4000,4000], loci=[2])
    InitByFreq(pop, [.2, .8])
    Stat(pop, alleleFreq=[0])
    af = pop.dvars().alleleFreq[0][1]  # ~.2
    # migrate and check if allele frequency changes
    Migrate(pop, mode=MigrByProbability, 
      fromSubPop = [0], toSubPop = [1,2], 
      rate = [.05, .05] )
    Stat(pop, alleleFreq=[0])
    assert pop.dvars().alleleFreq[0][1] == af
    Migrate(pop, mode=MigrByProbability, 
      rate = [ [0, .05, .05],
               [0.025, 0, 0],
               [0.025, 0, 0] ])
    assert pop.dvars().alleleFreq[0][1] == af
    Migrate(pop, mode=MigrByProbability, 
      rate = [ [0, .25, .25],
               [0.25, 0, 0],
               [0, 0.25, 0] ])
    assert pop.dvars().alleleFreq[0][1] == af
       
  def testPyMigrator(self):
    pop = population(subPop=[2,4,4], loci=[2,6])
    InitByFreq(pop, [.2,.4,.4])
    ind1 = pop.individual(0).arrGenotype()
    PyMigrate(pop, subPopID=carray('i',[2,2,0,0,2,2,1,1,1,1]))
    assert ind1 == pop.individual(6).arrGenotype()
    PyMigrate(pop, subPopID=carray('i',[0,0,2,2,2,2,1,1,1,1]))
    assert ind1 == pop.individual(2).arrGenotype()
    self.assertRaises(exceptions.ValueError, 
      PyMigrate, pop, carray('i',[0,0,2,2,2,2,1,1]))
    
if __name__ == '__main__':
  unittest.main()
   


