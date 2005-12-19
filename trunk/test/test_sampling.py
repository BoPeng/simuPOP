## #!/usr/bin/env python
## 
## ############################################################################
## #    Copyright (C) 2004 by Bo Peng                                         
## #    bpeng@rice.edu                                                        
## #                                                                          
## #    $LastChangedDate$          
## #    $Rev$                       
## #                                                                          
## #    This program is free software; you can redistribute it and/or modify  
## #    it under the terms of the GNU General Public License as published by  
## #    the Free Software Foundation; either version 2 of the License, or     
## #    (at your option) any later version.                                   
## #                                                                          
## #    This program is distributed in the hope that it will be useful,       
## #    but WITHOUT ANY WARRANTY; without even the implied warranty of        
## #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
## #    GNU General Public License for more details.                          
## #                                                                          
## #    You should have received a copy of the GNU General Public License     
## #    along with this program; if not, write to the                         
## #    Free Software Foundation, Inc.,                                       
## #    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
## ############################################################################
## 
## # test sampling issues
## #
## 
## from simuUtil import *
## 
## simu = simulator(
##     population(subPop=[100,200], ploidy=2, loci=[5,10],
##       ancestralDepth=1, maxAllele=9),
##     randomMating(numOffsprings=2))
## simu.evolve([
##     stat( alleleFreq=[0,1], genoFreq=[0,1]),
##     migrator(rate=[[0.1,0.1],[0.1,0.1]]),
##     basicPenetrance(locus=0,
##       penetrance={'1-1':0,'1-2':.7,'2-2':1}),
##     parentsTagger(),
##     ],
##  preOps=[  initByFreq(alleleFreq=[.2,.8], atLoci=[0]),
##     initByFreq(alleleFreq=[.2]*5, atLoci=range(1, simu.totNumLoci()))   ],
##  end=4
## )
## pop = simu.getPopulation(0)
## Dump(pop, ancestralPops=True)
## 
## # random sampling
## s = RandomSample(pop, 10)
## Dump(s[0])
## 
## s = RandomSample(pop,[2,8])
## Dump(s[0])
## 
## # case control sampling.
## s = CaseControlSample(pop, 10,10)
## Dump(s[0])
## listVars(s[0].dvars())
## 
## s =  CaseControlSample(pop, [5,5],[5,5])
## Dump(s[0])
## listVars(s[0].dvars())
## 
## # find sibpairs
## s=AffectedSibpairSample(pop)
## Dump(s[0])
## listVars(pop.dvars())
## 
## # from each subpop
## s = AffectedSibpairSample(pop,10)
## Dump(s[0], ancestralPops=True)
## listVars(s[0].dvars())
## s = AffectedSibpairSample(pop,[5,5])
## Dump(s[0])
## listVars(s[0].dvars())
## # 
## Dump(s[0], ancestralPops=True)
## 
## # save to linkage format
## # SaveLinkage(s[0], output='s0')
## 
## SaveRandFam(s[0], output='s0')
