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


# get options
from simuOpt import simuOptions

if simuOptions['Optimized'] == False and simuOptions['LongAllele'] == False:
  from simuPOP_std import *
elif simuOptions['Optimized'] == True and simuOptions['LongAllele'] == False:
  from simuPOP_op import * 
elif simuOptions['Optimized'] == False and simuOptions['LongAllele'] == True:
  from simuPOP_la import *
else:
  from simuPOP_laop import *

if not simuOptions['Quiet']:
  print "simuPOP : Copyright (c) 2004-2005 Bo Peng"
  # compile date, compiler etc are macros that are replaced during compile time.
  print ("Version %s (Revision %d, %s) for Python %s" % (simuVer(), simuRev(), compileDate(),
    compilePyVersion() ))
  print compileCompiler()
  print "Random Number Generator is set to", rng().name()
  print "Maximum allele number per locus is %d." % MaxAllele
  if optimized():
    print "You are running in optimized mode at maximum speed."
  else:
    print "You are running in standard mode with strict boundary check etc."
  print "For more information, please visit http://simupop.sourceforge.net,"
  print "or email simupop-list@lists.sourceforge.net (subscription required)."

