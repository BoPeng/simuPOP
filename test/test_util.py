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
## # test utilities
## #
## 
## from simuPOP import *
## from simuUtil import *
## 
## simu = simulator(population(subPop=[100,200], loci=[4,3]),
##   randomMating(), rep=2)
## 
## count = stat(alleles=range(0,5),
##                 genotypes=[(1,1,2),(2,1,3)],
##                 haplotypes=[(1,2,1,2),(1,2,2,3)])
## 
## ps = popStat([StatNumOfMale, StatNumOfFemale, StatPopSize])
## 
## simu.apply([
##   initByFreq([.3,.5,.2]),
##   ps, count, hetero(2) ]
## )
