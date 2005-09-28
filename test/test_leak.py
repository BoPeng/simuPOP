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

# test memory leak in various cases

# use command:
# 
#  valgrind --tool=memcheck --leak-check=yes python test_leak.py  
# 
# interactively, loading simuPOP itself will leak 1049 in 6 blocks
# this seems to be a python/swig problem.

# using script, everything is fine. 
import sys
import types


def rc(o):
  # one for parameter of rc, one for getrefcount(o)
  return sys.getrefcount(o)-3 

from simuPOP import *
from simuUtil import *

## 
pop = population(subPop=[50,10,10], ploidy=2,loci=[2,4])
simu = simulator(pop, randomMating(), rep=3)

# any operator?
op = dumper()
op1 = pyEval("gen")
op2 = pyExec("s=1")
 
op3 = initByFreq([.2,.3,.5], indRange=[[1,5],[6,10]])
op4 = stat(alleleFreq=range(0,pop.totNumLoci()))
op5 = stat(LD=[[1,2],[3,4]])
op6 = stat(Fst=[1,2])
# no leak till now
simu.apply([op3,op4])

def simulate():
  pop = population(subPop=[50,10,10], ploidy=2,loci=[2,4])
  InitByFreq(pop, [.2,.3,.5])  
  simu = simulator(pop, randomMating(), rep=3)
  
  # any operator?
  op = dumper()
  op1 = pyEval("gen")
  op2 = pyExec("s=1")
  
  op3 = initByFreq([.2,.3,.5], indRange=[[1,5],[6,10]])
  op4 = stat(alleleFreq=range(0,pop.totNumLoci()))
  op5 = stat(LD=[[1,2],[3,4]])
  op6 = stat(Fst=[1,2])
  
  # no leak till now
  simu.evolve([op6],end=10)
  return simu

##
for i in range(0,1):
  simu = simulate()

#for i in range(1,10):
#  simu = simulator(pop, randomMating(), rep=3)
#  simu.apply([ initByFreq([.2,.3,.5], indRanges=[[1,5],[6,10]]),
#    stat(LD=[[1,2]])])

 
#a = simu.vars(0)
#del simu
#print sys.getrefcount(a)

pop=population(size=1000, loci=[2,3])
InitByFreq(pop, [.2,.5,.3])
Stat(pop, alleleFreq=[2,4])
listVars(pop.vars())
pop.checkRefCount()
