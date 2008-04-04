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
import os, sys

if simuOptions['Optimized']:
    if simuOptions['AlleleType'] == 'short':
        from simuPOP_op import *
    elif simuOptions['AlleleType'] == 'long':
        from simuPOP_laop import *
    elif simuOptions['AlleleType'] == 'binary':
        from simuPOP_baop import * 
    else:
        from simuPOP_op import *
else:
    if simuOptions['AlleleType'] == 'short':
        from simuPOP_std import *
    elif simuOptions['AlleleType'] == 'long':
        from simuPOP_la import *
    elif simuOptions['AlleleType'] == 'binary':
        from simuPOP_ba import * 
    else:
        from simuPOP_std import *


if not simuOptions['Quiet']:
    print "simuPOP : Copyright (c) 2004-2008 Bo Peng"
    # compile date, compiler etc are macros that are replaced during compile time.
    if simuVer() == '9.9.9':
        # this is the subversion version of simuPOP
        print ("Developmental Version (%s) for Python %s" % (ModuleDate(), ModulePyVersion() ))
    else:
        # this is the released version
        print ("Version %s (Revision %d, %s) for Python %s" % (simuVer(), simuRev(), ModuleDate(),
            ModulePyVersion() ))
    print ModuleCompiler()
    print "Random Number Generator is set to %s with random seed 0x%08x" % (rng().name(), rng().seed())
    # MaxAllele + 1 since 0 is one of the allelic states
    if Optimized():
        print "This is the optimized %s allele version with %d maximum allelic states." % (AlleleType(), MaxAllele()+1)
    else:
        print "This is the standard %s allele version with %d maximum allelic states." % (AlleleType(), MaxAllele()+1)
    print "For more information, please visit http://simupop.sourceforge.net,"
    print "or email simupop-list@lists.sourceforge.net (subscription required)."

    if simuOptions['Debug'] != []:
        for g in simuOptions['Debug']:
            if g not in ['', None]:
                print "Turn on debug '%s'" % g
                TurnOnDebug(g)

