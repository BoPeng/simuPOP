#!/usr/bin/env python
#
# testing operator behaviors for simupoop.
# 
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision$
# $LastChangedDate$
#

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
    def getActiveGens(endGen=20, *args, **kwargs):
      d = opRecorder(*args, **kwargs)
      simu = simulator(population(), noMating())
      simu.evolve(ops=[d], end=endGen)
      return simu.population(0).dvars().hist
    self.assertEqual( getActiveGens(begin=2, end=10), 
      range(2,11))
    self.assertEqual( getActiveGens(begin=2, end=10, step=2), 
      range(2,11,2))
    self.assertEqual( getActiveGens(begin=2, step=2), 
      range(2,22,2))
    self.assertEqual( getActiveGens(step=2), range(0,22,2))
    self.assertEqual( getActiveGens(), range(0,21))
    self.assertEqual( getActiveGens(at=[2,5,9]), [2,5,9])
    self.assertEqual( getActiveGens(at=[2,5,-1]), [2,5,20])
    self.assertEqual( getActiveGens(begin=-10), range(11,21))
    # 20=-1, 16=-5
    self.assertEqual( getActiveGens(begin=-10, end=-5), range(11,17))
    # 
    self.assertEqual( getActiveGens(begin=-10, step=2, end=-5), range(11,17,2))
    self.assertRaises( exceptions.ValueError,
      getActiveGens, begin=-10, step=-3, end=-5 )
    
  def testGroup(self):
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
    
  def testOutput(self):
    simu = simulator(
        population(size=10, ploidy=2, loci=[2, 3 ]), 
        randomMating(), rep=5)
    # this will not change generation number
    simu.apply([ initByFreq([.2, .3, .5]),
                 dumper(alleleOnly=1)])
    # specify output is kind of complicated (because we sometimes need to 
    # let ops output to the same file.)
    #
    # here is the general rule: (details please refer to online documents.)
    #
    #   ">file.txt"  output to this file, once done, close.
    #                if mutliple ops output to the same file
    #                only the last one will succeed.
    #   ">>file.txt" output to this file, append, close when evolve stop.
    #   ">>file%rep.txt" different filename for each replicates.
    #                %rep will be substituted by replicate number.
    #   ">>file%gen.txt" different filename for different generation.
    #                %grp will be substituted by generation number
    #   ">>%{sim}%rep.txt variables (set by setStringVar etc) can be used.
    #   "|pipeName"   A pipe can be used. I.e., everything will be
    #                written to a memory file. This is useful only
    #                when the memory file will be read by a
    #                DataProvider.
    #
    # for example:
    #
    simu.setGen(0)
    # note the use of sep
    # count allele at each locus
    # header format: replicate_loci_allele
    simu.evolve([ stat(alleleFreq=[2],
                                output=">", sep=",",rep=1),
                  output("\n", rep=REP_LAST)
               ], end=10)
    
    #  if you would like to supress the last comma, use a
    # \b backspace character. However, this is useful
    # only when the output is standard output
    # (If you write to a file, ',' and \b will be written)
    simu.setGen(0)
    simu.evolve([ stat(alleleFreq=[0],output=">>", sep=',',rep=2),
                  output("\b\n", rep=REP_LAST)
               ], end=10)
    
    # output to another file
    simu.setGen(0)
    simu.evolve([ stat(alleleFreq=[0],output="a.txt"),
                  output("\n", rep=REP_LAST, output="a.txt")
               ], end=10)
    
    # Nothing!! (this usage is totally python)
    print open("a.txt").read()
    
    # you should use >> so that output will
    # append instead of replace output from stat
    # In the previous case, only the last output write to
    # a.txt, i.e., only a \n was written.
    #
    simu.setGen(0)
    simu.evolve([ stat(alleleFreq=[0],
                  output=">>a.txt", rep=2),
                  output("\n", rep=REP_LAST, output=">>a.txt")
               ], end=10)
    
    
    print open("a.txt").read()
    
    # note that once you specify "append mode", you do not have
    # to do this each time.
    simu.setGen(0)
    simu.evolve([ stat(alleleFreq=[0],output=">>a.txt", sep=','),
                  output("\b\n", rep=REP_LAST, output="a.txt")
               ], end=10)
    
    print open("a.txt").read()
    
    # replicate specific output
    simu.setGen(0)
    simu.evolve([ stat(alleleFreq=[0],
       outputExpr=r"'>>a'+str(rep)+'.txt'",
                 sep=','),
       output("\b\n", outputExpr=r"'>>a'+str(rep)+'.txt'")
      ], end=10)
    
    print open("a1.txt").read()
    print open("a2.txt").read()
    
    # generation specific
    simu.setGen(0)
    simu.evolve([ stat(alleleFreq=[0],
                                outputExpr=r"'>>b'+str(rep)+'.txt'" ,sep=',')
               ], end=10)
    
    # note that no header is output.
    #print open("b1.txt").read()
    #print open("b2.txt").read()
    
    # you can also specify "no output" 
    simu.setGen(0)
    simu.evolve([ stat(alleleFreq=[0],output="")], end=10)
    
    # you can actually combine replication, generation, group specific
    # outputs....
    # to avoid generting too many files, no example is given here.
    
    # remove all generated files
    
    import os
    os.remove("a.txt")
    for n in range(5):
        os.remove( 'a%d.txt' % n)
    #for n in range(0,11):
    #    os.remove( 'b%d.txt' % n)
   
if __name__ == '__main__':
  unittest.main()
