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

# This script test the mutation, drift equilibrium
# for a biallelic case
#

# infinite allele model
import simuOpt
simuOpt.setOptions(longAllele=True)
from simuPOP import *
from simuUtil import *
from simuRPy import *


N, numRep, mu = 5000, 4, 1e-4
simu = simulator(
  population(size=N, ploidy=2, loci=[1]),
  randomMating(),
  rep=numRep)
simu.evolve(
  preOps = [ initByFreq([.1]*10) ],
  ops = [
    #kamMutator(rate=mu, maxAllele=1000),
    smmMutator(rate=mu, maxAllele=1000),
    stat(homoFreq=[0]),
    varPlotter('homoFreq[0]', numRep=numRep, byRep=0, ylim=[0,1],
      update=50, win=2000),
    pyEval("gen, homoFreq[0]", rep=REP_LAST, step=50),
    endl(rep=REP_LAST, step=50),
    pause(stopOnKeyStroke=True),
    ]
  )
