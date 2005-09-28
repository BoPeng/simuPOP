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

#!/usr/bin/env python
#
# Purpose:
#   selection.
# 
# Author:
#   Bo Peng (bpeng@rice.edu)
#
# load module.
from simuPOP import *

#
# 1. basic selection model
#
#    specify relative fitness: w11, w12/w21, w22
#
#

# random mating with subPopulations, multiple offspring
setArrayVar("popSize", [1000,1000])
simu = simulator(
    population(size=2000, ploidy=2, loci=[1], subPop=[1000,1000]),
    randomMating(numOffsprings=4, newSubPopSizeExpr="%popSize=%popSize*1.02"))

# initial allele frequency
count =  stat(
      alleleFreq=[0[0,1,2]],
      genotypes=[[0,1,1,1,2,2,2]]) 
sel = basicSelector(locus=0, fitness={101:0.6, 102:1, 202:.8},
                    begin=100)
visual = matlabPlotter(
    varDataSource("%gen%genotypeFreq0_sp0[102]%genotypeFereq0_sp1[102]"),
    win=300, rep=-1,
    update=5)
visual.addPostPlotCmd("ylim([0,1.3])")

simu.setGen(0)
simu.evolve([
    count,
    sel,
    visual,
    ],
    preOps [initByFreq(alleleFreq=[.2])],
    end=10)

#
# test pyMating
#

def mateFunc(pop):
  # keep population size
  idx=[0]*(pop.popSize()*2)
  i = 0
  for sp in range(pop.numSubPop()):
    if pop.subPopSize(sp) >= 2:
      for ind in range(pop.subPopSize(sp)):
        idx[i] = pop.subPopBegin(sp)
        idx[i+1] = pop.subPopBegin(sp)+1
        i += 2
    elif pop.subPopSize(sp) == 1:
        idx[i] = pop.subPopBegin(sp)
        idx[i+1] = idx[i]
        i += 2
  print idx
  return idx

pop = population(subPop=[4,6], loci=[1,2])
mateFunc(pop)


simu = simulator(pop, pyMating(mateFunc))
simu.evolve([
    parentsTagger(),
    ],
    preOps = [initByFreq(alleleFreq=[.2, .8])],
    end=10)

Dump(simu.population(0))
