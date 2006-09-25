#!/usr/bin/env python
#
# testing operator behaviors for simupoop.
# 
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys

# record active generations in pop.dvars().hist
def genRecorder(pop):
  try:
    pop.dvars().hist.append(pop.gen())
  except:
    pop.dvars().hist = [pop.gen()]
  return True
  
def opRecorder(*args, **kwargs):
  return pyOperator(func=genRecorder, *args, **kwargs)
  
class TestOperator(unittest.TestCase):

  def testMemberFunctions(self):
    'Testing common operator member functions'
    d = dumper()
    self.assertEqual(d.canApplyPreMating(), False)
    self.assertEqual(d.canApplyDuringMating(), False)
    self.assertEqual(d.canApplyPostMating(), True)
    # apply to any replicate
    self.assertEqual(d.applicableReplicate(), REP_ALL)
    # apply to any group
    self.assertEqual(d.applicableGroup(), GRP_ALL)
    # set grp
    d.setApplicableGroup(GRP_ALL)
    self.assertEqual(d.applicableGroup(), GRP_ALL)
    #
    d.setApplicableReplicate(REP_LAST)
    self.assertEqual(d.applicableReplicate(), REP_LAST)
  
  def testActiveGen(self):
    'Testing active generation specifications'
    def getActiveGens(endGen=20, *args, **kwargs):
      d = opRecorder(*args, **kwargs)
      simu = simulator(population(), noMating())
      simu.evolve(ops=[d], end=endGen)
      return simu.population(0).dvars().hist
    self.assertEqual(getActiveGens(begin=2, end=10), 
      range(2,11))
    self.assertEqual(getActiveGens(begin=2, end=10, step=2), 
      range(2,11,2))
    self.assertEqual(getActiveGens(begin=2, step=2), 
      range(2,22,2))
    self.assertEqual(getActiveGens(step=2), range(0,22,2))
    self.assertEqual(getActiveGens(), range(0,21))
    self.assertEqual(getActiveGens(at=[2,5,9]), [2,5,9])
    self.assertEqual(getActiveGens(at=[2,5,-1]), [2,5,20])
    self.assertEqual(getActiveGens(begin=-10), range(11,21))
    # 20=-1, 16=-5
    self.assertEqual(getActiveGens(begin=-10, end=-5), range(11,17))
    # 
    self.assertEqual(getActiveGens(begin=-10, step=2, end=-5), range(11,17,2))
    self.assertRaises( exceptions.ValueError,
      getActiveGens, begin=-10, step=-3, end=-5 )
    
  def testGroup(self):
    'Testing group related functions'
    simu = simulator(population(), noMating(), rep=3)
    simu.setGroup([1,1,2])
    simu.evolve(
      ops = [opRecorder(grp=1)], 
      end=10
    )
    self.assertEqual(simu.population(0).dvars().hist, range(11))
    self.assertEqual(simu.population(1).dvars().hist, range(11))
    try:
      simu.population(2).dvars().hist
    except exceptions.AttributeError:
      pass

  def testReplicate(self):
    'Testing replicate related functions'
    simu = simulator(population(), noMating(), rep=3)
    simu.evolve(
      ops = [opRecorder(rep=REP_LAST)], 
      end=10
    )
    try:
      simu.population(0).dvars().hist
    except exceptions.AttributeError:
      pass
    try:
      simu.population(1).dvars().hist
    except exceptions.AttributeError:
      pass
    self.assertEqual(simu.population(2).dvars().hist, range(11))
    
  def assertFileContent(self, file, text):
    f = open(file)
    t = f.read()
    f.close()
    self.assertEqual(t, text)
    
  def testOutput(self):
    'Testing output specifications'
    simu = simulator( population(), 
        noMating(), rep=5)
    simu.evolve([
      output("a", output=">a.txt"),
      ], end=10)
    # although everyone have written to this file,
    # only the last one will be kept
    self.assertFileContent("a.txt", 'a')
    os.remove('a.txt')
    #
    # you can ignore >
    simu.setGen(0)
    simu.evolve([
      output("a", output="a.txt"),
      ], end=10)
    # although everyone have written to this file,
    # only the last one will be kept
    self.assertFileContent("a.txt", 'a')
    os.remove('a.txt')
    #
    # >>
    simu.setGen(0)
    simu.evolve([
      output("a", output=">>a.txt"),
      ], end=10)
    # a is appended 5 rep * 11 generations
    self.assertFileContent("a.txt", 'a'*55)
    os.remove('a.txt')
    #
    # rep = ...
    simu.setGen(0)
    simu.evolve([
      output("a", output=">>a.txt", rep=REP_LAST),
      ], end=10)
    # a is appended 5 rep * 11 generations
    self.assertFileContent("a.txt", 'a'*11)
    # if we use >>>, append to the end
    simu.setGen(0)
    simu.setGroup([0,0,1,1,1])
    simu.evolve([
      output("b", output=">>>a.txt", grp=1),
      ], end=10)
    # a is appended 5 rep * 11 generations
    self.assertFileContent("a.txt", 'a'*11+'b'*33)
    os.remove('a.txt')
    #
    # now, we can use eval instead of output
    simu.setGen(0)
    simu.setGroup([0,0,1,1,1])
    simu.evolve([
      pyEval("gen", output=">>a.txt", grp=1),
      ], end=10)
    # a is appended 5 rep * 11 generations
    self.assertFileContent("a.txt", 
      ''.join( [ str(x)*3 for x in range(11)] ))
    os.remove('a.txt')

  def testOutputExpr(self):
    'Testing the usage of output expression'
    simu = simulator( population(), 
      noMating(), rep=5)
    # each replicate
    simu.evolve([
      output("a", outputExpr="'rep%d.txt'%rep"),
      ], end=10)
    # although everyone have written to this file,
    # only the last one will be kept
    for i in range(5):
      self.assertFileContent("rep%d.txt"%i, 'a')
      os.remove('rep%d.txt'%i)
    #
    # you can ignore >
    simu.setGen(0)
    simu.evolve([
      output("a", outputExpr="'>rep%d.txt'%rep"),
      ], end=10)
    # although everyone have written to this file,
    # only the last one will be kept
    for i in range(5):
      self.assertFileContent("rep%d.txt"%i, 'a')
      os.remove('rep%d.txt'%i)
    #
    # >>
    simu.setGen(0)
    simu.evolve([
      output("a", outputExpr="'>>rep%d.txt'%rep"),
      ], end=10)
    # a is appended 1 rep * 11 generations
    for i in range(5):
      self.assertFileContent("rep%d.txt"%i, 'a'*11)
      os.remove('rep%d.txt'%i)
    # each generation?
    simu.setGen(0)
    simu.evolve([
      output("a", outputExpr="'>>gen%d.txt'%gen"),
      ], end=10)
    # a is appended 1 rep * 11 generations
    for i in range(11):
      self.assertFileContent("gen%d.txt"%i, 'a'*5)
      os.remove('gen%d.txt'%i)
   
if __name__ == '__main__':
  unittest.main()
