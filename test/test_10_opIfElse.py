#!/usr/bin/env python
#
# This is a unittest file for ifElse operator
#
# Bo Peng (bpeng@rice.edu)
# 
# Last Modified: Sep, 2005
# 

from simuPOP import *
from simuRPy import *
import unittest, time

class TestIfElseOperator(unittest.TestCase):
   
  def setUp(self):
    ' set up a simulator'
    self.simu = simulator(population(100, loci=[2]), 
      randomMating())
  
  def testIfElseOperator(self):
    # now if we want to inject a mutation whenever fixation happens
    self.simu.evolve(
      preOps = [ initByFreq([.5,.5]) ],
      ops = [
        # count number of allels at this locus
        stat(alleleFreq=[0]),
        # inject 50% of allele 2 if this allele get lost
        ifElse('alleleFreq[0][1]==1.', kamMutator(rate=.5, maxAllele=2, atLoci=[0]) ),
        # the other way around?
        ifElse('alleleFreq[0][1]==0.', kamMutator(rate=.5, maxAllele=2, atLoci=[0]) ),
        # print allele frequencies, we should see patterns of reaching zero
        # jump to .5
        # pyEval(r'"%.2f\n" % alleleFreq[0][1]')
        varPlotter('alleleFreq[0][1]', update=100, ylim=[0,1],
          title='allele freq adjust to .5 when one allele get lost')
      ],
      end=1000
    )
    # hold the image for 10 seconds, 
    time.sleep(10)
    
if __name__ == '__main__':
  unittest.main()
   


