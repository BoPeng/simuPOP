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
# Purpose:
#  testing of interfaces of mutators  of simuPOP
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# June, 2004.
#
# Important:  set path to simuPOP.py
#
#  if you do not set PYTHONPATH, you can uncomment
#  the following two lines and specify path to
#  simuPOP.py here.
#

from simuPOP import *
from simuUtil import *

import unittest

class TestMutator(unittest.TestCase):
  def testMutator(self):
    simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
      randomMating(), rep=5)
    simu.apply( [ initByFreq([.2,.8])])
    simu.evolve([ kamMutator(rate=0.1)], end=200)

  def testPyMutator(self):
    pop = population(size=10, loci=[2])
    # cutom mutator
    def mut(x):
      return 8
    m = pyMutator(rate=1, func=mut)
    m.apply(pop)
    assert pop.individual(0).allele(0) == 8, \
      "PyMutator failed"
    
  def testMutationCount(self):
    N = 10000
    r = 0.001
    G = 100
    pop = population(size=N, ploidy=2, loci=[5])
    simu = simulator(pop, randomMating())
    mut = kamMutator(rate=r, maxAllele=10)
    simu.evolve( [mut], end=G)
    print mut.mutationCounts()
    assert abs( mut.mutationCount(0) - 2*N*r*G) < 200, \
      "Number of mutation event is not as expected."
    
  def testPointMutator(self):    
    # test point mutator
    pop = population(size=10, ploidy=2, loci=[5])
    InitByValue(pop, value=[[1]*5, [2]*5], proportions=[.3,.7])
    PointMutate(pop, inds=[1,2,3], toAllele=3, atLoci=[1,2])
    assert pop.individual(1).allele(1,0) == 3
    assert pop.individual(1).allele(1,1) != 3
    PointMutate(pop, inds=[1,2,3], atPloidy=[1], 
      toAllele=4, atLoci=[1,2])
    assert pop.individual(1).allele(2,1) == 4
    assert pop.individual(1).allele(2,0) != 4

if __name__ == '__main__':
  unittest.main()
