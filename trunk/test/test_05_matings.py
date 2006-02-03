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
## from simuPOP import *
## 
## #
## # 1. map selection model
## #
## #    specify relative fitness: w11, w12/w21, w22
## #
## #
## 
## # random mating with subPopulations, multiple offspring
## setArrayVar("popSize", [1000,1000])
## simu = simulator(
##     population(size=2000, ploidy=2, loci=[1], subPop=[1000,1000]),
##     randomMating(numOffsprings=4, newSubPopSizeExpr="%popSize=%popSize*1.02"))
## 
## # initial allele frequency
## count =  stat(
##       alleleFreq=[0[0,1,2]],
##       genotypes=[[0,1,1,1,2,2,2]]) 
## sel = mapSelector(locus=0, fitness={101:0.6, 102:1, 202:.8},
##                     begin=100)
## visual = matlabPlotter(
##     varDataSource("%gen%genotypeFreq0_sp0[102]%genotypeFereq0_sp1[102]"),
##     win=300, rep=-1,
##     update=5)
## visual.addPostPlotCmd("ylim([0,1.3])")
## 
## simu.setGen(0)
## simu.evolve([
##     count,
##     sel,
##     visual,
##     ],
##     preOps [initByFreq(alleleFreq=[.2])],
##     end=10)
## 
## #
## # test pyMating
## #
## 
## def mateFunc(pop):
##   # keep population size
##   idx=[0]*(pop.popSize()*2)
##   i = 0
##   for sp in range(pop.numSubPop()):
##     if pop.subPopSize(sp) >= 2:
##       for ind in range(pop.subPopSize(sp)):
##         idx[i] = pop.subPopBegin(sp)
##         idx[i+1] = pop.subPopBegin(sp)+1
##         i += 2
##     elif pop.subPopSize(sp) == 1:
##         idx[i] = pop.subPopBegin(sp)
##         idx[i+1] = idx[i]
##         i += 2
##   print idx
##   return idx
## 
## pop = population(subPop=[4,6], loci=[1,2])
## mateFunc(pop)
## 
## 
## simu = simulator(pop, pyMating(mateFunc))
## simu.evolve([
##     parentsTagger(),
##     ],
##     preOps = [initByFreq(alleleFreq=[.2, .8])],
##     end=10)
## 
## Dump(simu.population(0))
 

# set new size
def newSize(gen,sz):
  os = range(0,len(sz))
  for i in range(0, len(os)):
    os[i] = sz[i] * 1.2
  return os

simu = simulator(population(subPop=[10,20,30]),
  randomMating(newSubPopSizeFunc=newSize))

simu.evolve(
  ops = [
    stat( popSize = 1),
    pyEval('subPopSize'),
    endl()
    ],
  end = 10
  )

# ancestral population

simu = simulator(population(subPop=[10,20,30], ancestralDepth=3),
  randomMating(newSubPopSizeFunc=newSize))

simu.evolve(
  ops = [
    stat( popSize = 1),
    pyEval('subPopSize'),
    endl()
    ],
  end = 10
  )

Dump( simu.population(0), ancestralPops=True, infoOnly=True)


# test numOffsprings
simu = simulator(population(10, loci=[2]), randomMating(numOffsprings=2))
simu.step(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()])

# numOffspringsFunc
def nos(gen):
  return gen+1

simu = simulator(population(10, loci=[2]), randomMating(numOffspringsFunc=nos))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# MATE_NumOffspringsEachFamily
import random
def nos(gen):
  return random.randrange(1,4)

simu = simulator(population(10, loci=[2]), randomMating(numOffspringsFunc=nos))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffspringsFunc=nos, mode=MATE_NumOffspringsEachFamily))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# MATE_GeometricDistribution
simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffsprings=.3, mode=MATE_GeometricDistribution))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# MATE_BinomialDistribution
simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffsprings=.3, maxNumOffsprings=5, mode=MATE_BinomialDistribution))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# MATE_PoissonDistribution
simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffsprings=.3, mode=MATE_PoissonDistribution))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# func with distribution
def p(gen):
  return min(gen/10., .9)

turnOnDebug(DBG_MATING)
simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffspringsFunc=p,
    maxNumOffsprings=5,
    mode=MATE_BinomialDistribution))
simu.evolve(preOps=[initByFreq([.4,.6])],
  ops=[parentsTagger()], end=10)


