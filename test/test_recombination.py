#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of mutators  of simuPOP
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#

from simuPOP import *
from simuUtil import *

#
# 
#
# verify table 4 of H&C 3nd edition P49
N = 10000
r = 0.1

genoDad = [[1,1],[1,1],[1,1],[1,1],[1,2],[1,2],[1,2],[2,1],[2,1],[2,2]]
genoMom = [[1,1],[1,2],[2,1],[2,2],[1,2],[2,1],[2,2],[2,1],[2,2],[2,2]]

#turnOnDebug(DBG_RECOMBINATOR)
for i in range(0, len(genoDad)):
  pop = population(size=N, loci=[2])
  InitByValue(pop, value=[genoDad[i], genoMom[i]], proportions=[0.5,0.5])
  simu = simulator(pop, randomMating())
  simu.step(
      [ recombinator(rate=r),
        stat(haploFreq=[[0,1]]) ], steps=2)
  # listVars(simu.dvars(0).haploFreq)
  
# table is correct.
# 1/2*(1-r)=0.5*0.9=0.45


#
#
# verify Dn=(1-r)^n D0
#
i=1
pop = population(size=N, loci=[5,7,8], sexChrom=True)
InitByValue(pop, value=[genoDad[i], genoMom[i]], proportions=[0.5,0.5])
simu = simulator(pop, randomMating())
simu.evolve(
    [ recombinator(rate=r, maleRate=r/2., maleAfterLoci=[2,6,7]),
      stat(haploFreq=[[0,1]]) ],
    end=100)
