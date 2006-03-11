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

if simuOptions['Optimized'] == False and simuOptions['AlleleType'] == 'short':
  from simuPOP_std import *
elif simuOptions['Optimized'] == True and simuOptions['AlleleType'] == 'short':
  from simuPOP_op import * 
elif simuOptions['Optimized'] == False and simuOptions['AlleleType'] == 'long':
  from simuPOP_la import *
elif simuOptions['Optimized'] == True and simuOptions['AlleleType'] == 'long':
  from simuPOP_laop import *
elif simuOptions['Optimized'] == False and simuOptions['AlleleType'] == 'binary':
  from simuPOP_ba import * 
elif simuOptions['Optimized'] == True and simuOptions['AlleleType'] == 'binary':
  from simuPOP_baop import *
else:
  print "Warning: options unrecognized (AlleleType=%s, Optimized=%d). Use standard library. " \
    % (simuOptions['AlleleType'], simuOptions['Optimized'])
  from simuPOP_std import *

if simuOptions['Debug'] != []:
  for g in simuOptions['Debug']:
    TurnOnDebugWithName(g)

if not simuOptions['Quiet']:
  print "simuPOP : Copyright (c) 2004-2005 Bo Peng"
  # compile date, compiler etc are macros that are replaced during compile time.
  print ("Version %s (Revision %d, %s) for Python %s" % (simuVer(), simuRev(), compileDate(),
    compilePyVersion() ))
  print compileCompiler()
  print "Random Number Generator is set to", rng().name()
  # MaxAllele + 1 since 0 is one of the allelic states
  print "This is the %s allele version with %d maximum allelic states." % (alleleType(), MaxAllele+1)
  if optimized():
    print "You are running in optimized mode at maximum speed."
  else:
    print "You are running in standard mode with strict boundary check etc."
  print "For more information, please visit http://simupop.sourceforge.net,"
  print "or email simupop-list@lists.sourceforge.net (subscription required)."

