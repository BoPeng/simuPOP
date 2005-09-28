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


# testing statistics calculation

from simuUtil import *
import sys

simu = simulator(population(subPop=[500,100,1000], ploidy=2,
   loci = [1]), randomMating()       )
 

init = initByValue(
  value = [[1,1],[1,2],[2,2],[1,1],[1,2],[2,2],[1,2],[1,2],[2,2]],
  indRange = [[0,125], [125,375],[375,500],[500,550],
    [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])

 
count = stat(
  popSize = 1,
  numOfMale = 1,
  numOfAffected = 1,
  alleleFreq = [0],
  heteroFreq = [0],
  homoFreq = [0],
  expHetero = [0],
  genoFreq = [0],
  Fst = [ 0 ]
  )

pene = basicPenetrance(locus=0, stage=PreMating,
  penetrance={'1-1':0,'1-2':1,'2-2':0} )


simu.apply([init, pene])

SaveFstat(simu.population(0), "p1.dat", maxAllele=2)
# Fst is compaared with result from Fstat.
#
#For locus : loc_0_
#Allele  Capf	 Theta	Smallf	 Relat	Relatc	 Sig_a	 Sig_b	 Sig_w
#    1	-0.103	 0.102	-0.229	 0.228
#    2	-0.103	 0.102	-0.229	 0.228
#  All	-0.103	 0.102	-0.229	 0.228	 0.373	 0.051	-0.102	 0.550
#


simu.apply([count])

listVars(simu.dvars(0) )


# test number of alleles
pop = population(subPop=[500,1000,1000], ploidy=2, loci = [20])
InitByFreq(pop, [.2]*5)
Stat(pop, LD=[[1,2],[2,3,4,5]])
listVars(pop.dvars())


tests = [
  ['','s=simu.dvars(0)','',''],
  ['','listVars(s)','',''],
  ['len(s.subPop)','',3,'']
]


# perform test
for i in range(0, len(tests)):
  print "Testing ", i
  testExpr( expr=tests[i][0], stmts=tests[i][1],\
    res = tests[i][2], excpt = tests[i][3], d=vars() )

# test haplotype frequency
pop = population(subPop=[5,10], ploidy=2, loci = [10])
InitByValue(pop, value=[[1]*10,[2]*10,[3]*10], proportions=[.2,.3,.5])
Stat(pop, haploFreq=[[0,1],[0,5]], homoFreq=[0,3], heteroFreq=[0,3])
listVars(pop.dvars())


## # check LD
## p = validPops[0].dvars()
## Stat(validPops[0], LD=[[0,1]])
## print p.AvgLD_prime[0][1]

## DP = 0.
## for i in range(1,3):
##   for j in range(1,6):
##     PAB = p.haploFreq['0-1'].setdefault("%d-%d"%(i,j),0)
##     PA = p.alleleFreq[0][i]
##     PB = p.alleleFreq[1][j]
##     DAB = PAB-PA*PB
##     if DAB < 0:
##       DMAX = min(PA*PB, (1-PA)*(1-PB))
##     else:
##       DMAX = min(PA*(1-PB), PB*(1-PA))
##     print DAB, DMAX, DAB/DMAX
##     DP += PA*PB*abs(DAB/DMAX)
## print DP

pop = population(subPop=[5,10], ploidy=2, loci = [20]*20)
InitByFreq(pop, [.2,.3,.5])
Stat(pop, LD=[[4,5],[100,120]])
listVars(pop.dvars(), subPop=False)
