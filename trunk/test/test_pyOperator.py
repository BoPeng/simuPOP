#!/usr/bin/env python

############################################################################
#    Copyright (C) 2004 by Bo Peng                                         
#    bpeng@rice.edu                                                        
#                                                                          
#    $LastChangedDate$          
#    $Rev$                       
#                                                                          
#    This program is free software; you can redistribute it and/or modify  
#    it under the terms of the GNU General Public License as published by  
#    the Free Software Foundation; either version 2 of the License, or     
#    (at your option) any later version.                                   
#                                                                          
#    This program is distributed in the hope that it will be useful,       
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
#    GNU General Public License for more details.                          
#                                                                          
#    You should have received a copy of the GNU General Public License     
#    along with this program; if not, write to the                         
#    Free Software Foundation, Inc.,                                       
#    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
############################################################################
#
# testing puOperator
# 
# load module.
from simuPOP import *
from simuRPy import *
import unittest, time
    
class TestPyOperator(unittest.TestCase):
  def setUp(self):
    self.pop = population(size=10000, ploidy=2, 
      loci=[2, 3])

  # define a function
  def myFunc(self, pop):
    Stat(pop, alleleFreq=[0])
    assert abs(pop.dvars().alleleFreq[0][1] - 0.2) < 0.05
    return True
    
  def testSimpleFunc(self):
    InitByFreq(self.pop, [.2, .3, .5])
    simu = simulator(self.pop, randomMating())
    simu.evolve( ops=[pyOperator(self.myFunc)], 
      end=2)
    
  def testCopyCone(self):
    op = pyOperator(self.myFunc)
    op1 = op
    op2 = op.clone()
    InitByFreq(self.pop, [.2, .3, .5])
    # all copied version are working fine.
    op.apply(self.pop)
    op1.apply(self.pop)
    op2.apply(self.pop)
    
  def myFuncWithParam(self, pop, para):
    ' para is (allele, freq) pair '
    Stat(pop, alleleFreq=[0])
    assert abs(pop.dvars().alleleFreq[0][para[0]] - para[1]) < 0.05
    return True

  def testFuncWithParam(self):
    InitByFreq(self.pop, [.2, .3, .5])
    simu = simulator(self.pop, randomMating())
    simu.evolve( ops=[
      pyOperator(func=self.myFuncWithParam, param=(1,.2)),
      pyOperator(func=self.myFuncWithParam, param=(2,.3)),
      pyOperator(func=self.myFuncWithParam, param=(3,.5))
      ], 
      end=2)

  def myFuncAsTerminator(self, pop):
    # print pop.gen()
    if pop.gen() == 3:
      return False
    else:
      return True

  def testTerminator(self):
    simu = simulator(self.pop, randomMating())
    simu.evolve( ops=[ pyOperator(self.myFuncAsTerminator) ],
      end = 10 )
    assert simu.gen() == 4
    
  def dynaMutator(self, pop, param):
    ''' this mutator mutate common loci with low mutation rate
    and rare loci with high mutation rate, as an attempt to
    bring allele frequency of these loci at an equal level.'''
    # unpack parameter
    (cutoff, mu1, mu2) = param;
    Stat(pop, alleleFreq=range( pop.totNumLoci() ) )
    for i in range( pop.totNumLoci() ):
      # 1-freq of wild type = total disease allele frequency
      if 1-pop.dvars().alleleFreq[i][1] < cutoff:
        KamMutate(pop, maxAllele=2, rate=mu1, atLoci=[i])
      else:
        KamMutate(pop, maxAllele=2, rate=mu2, atLoci=[i])
    return True

  def testDynaMutator(self):
    simu = simulator(self.pop, randomMating())
    simu.evolve(
      preOps = [ 
        initByFreq( [.6, .4], atLoci=[0,2,4]),
        initByFreq( [.8, .2], atLoci=[1,3]) ],
      ops = [ 
        pyOperator( func=self.dynaMutator, param=(.5, .1, 0) ),
        stat(alleleFreq=range(5)),
        varPlotter('(alleleFreq[0][2],alleleFreq[1][2])', 
          title='allele frequency at loci 0 and 1',
          update=5, varDim=2),
        ],
      end = 30
    )        
    print "The R window will be closed after five seconds..."
    time.sleep(5)
    
if __name__ == '__main__':
  unittest.main()
   


