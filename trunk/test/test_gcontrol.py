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

# testing for Genomic Control

from simuUtil import *
from simuRPy import *
pop = population(subPop=[3000,3000], ploidy=2,
   loci = [5,20,6], maxAllele=2)
simu = simulator(pop, randomMating(), rep=2)

# SNP markers
init = initByFreq([.3,.7])

count = stat(
  genoFreq = [ x for x in range(0,pop.totNumLoci())],
  )

# turnOnDebug(DBG_STATOR)
simu.apply([init, count,
  GC(case=0, control=1, loci=range(0, pop.totNumLoci())  )])
GControl(case=0,control=1,rep=0, loci=range(0,pop.totNumLoci()))

saveInFstatFormat(simu.population(0), "a1.dat")

# the results totally agree with what fstat got
simu.apply([ saveFstat( outputExpr=r"'ran%d.dat' % rep")])

stat(1).listVars(2)

