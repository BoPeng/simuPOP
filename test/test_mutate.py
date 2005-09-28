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
simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
                  randomMating(), rep=5)

simu.apply( [ initByFreq([.2,.8])])
d = dumper(alleleOnly=1, rep=5)
simu.apply([d])

simu.evolve([ kamMutator(rate=0.1)], end=200)


# have a look at a single mutator
# help(kamMutator)
m = kamMutator(rate=0.5, maxAllele=9)
print m.maxAllele()
m.setMaxAllele(5)
print m.rate()
# 
m.setRate([0.1,.02],[1,2])
print m.rate()
m.setRate(0.2)
print m.rate()

simu.setGen(0)
simu.evolve([m, dumper(step=5,rep=5)], end=12)

# it is easier to see mutation if no mating is invoolved.

simu = simulator(population(size=10, ploidy=2, loci=[2, 3]),
                  randomMating(), rep=3)
simu.evolve([kamMutator(rate=0.1), dumper(rep=3)], end=5)

pop = population(size=10)
d=dumper()
d.apply(pop)
# cutom mutator
def mut(x):
  return 8

m = pyMutator(rate=1, func=mut)
m.apply(pop)
d.apply(pop)


from simuUtil import *
# test point mutator
pop = population(size=10, ploidy=2, loci=[5])
InitByValue(pop, value=[[1]*5, [2]*5], proportions=[.3,.7])
Dump(pop)
PointMutate(pop, inds=[1,2,3], toAllele=3, atLoci=[1,2])
Dump(pop)


#
#
# test P113 of HC 3nd
#
