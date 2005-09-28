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

# testing fstat read/write/Fst calculation

from simuUtil import *
simu = simulator(population(subPop=[500,100,1000], ploidy=2,
   loci = [1], maxAllele=2), randomMating(), rep=2)

simu.apply( [
  initByValue(
    value = [[1,1],[1,2],[2,2],[1,1],[1,2],[2,2],[1,2],[1,2],[2,2]],
    indRange = [[0,125], [125,375],[375,500],[500,550],
    [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
 ] )

# now, write to a file
SaveFstat(simu.population(0), "p0.dat")
print open("p0.dat").read()
pop = LoadFstat("p0.dat")
dumper().apply(pop)

#
count = stat(
  popSize = 1,
  numOfMale = 1,
  alleleFreq = [0],
  heteroFreq = [0],
  genoFreq = [0],
  Fst = [ 0]
  )

simu.apply([count])
listVars(simu.dvars(0))
