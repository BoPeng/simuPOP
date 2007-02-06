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
    if simuOptions['MPI']:
        if simuOptions['AlleleType'] == 'short':
            from simuPOP_opmpi import *
        elif simuOptions['AlleleType'] == 'long':
            from simuPOP_laopmpi import *
        elif simuOptions['AlleleType'] == 'binary':
            from simuPOP_baopmpi import * 
        else:
            from simuPOP_opmpi import *
    else:
        if simuOptions['AlleleType'] == 'short':
            from simuPOP_op import *
        elif simuOptions['AlleleType'] == 'long':
            from simuPOP_laop import *
        elif simuOptions['AlleleType'] == 'binary':
            from simuPOP_baop import * 
        else:
            from simuPOP_op import *
else:
    if simuOptions['MPI']:
        if simuOptions['AlleleType'] == 'short':
            from simuPOP_mpi import *
        elif simuOptions['AlleleType'] == 'long':
            from simuPOP_lampi import *
        elif simuOptions['AlleleType'] == 'binary':
            from simuPOP_bampi import * 
        else:
            from simuPOP_mpi import *
    else:
        if simuOptions['AlleleType'] == 'short':
            from simuPOP_std import *
        elif simuOptions['AlleleType'] == 'long':
            from simuPOP_la import *
        elif simuOptions['AlleleType'] == 'binary':
            from simuPOP_ba import * 
        else:
            from simuPOP_std import *

if simuOptions['Debug'] != []:
    for g in simuOptions['Debug']:
        print "Turn on debug ", g
        TurnOnDebugWithName(g)

# seed rng() if necessay
if not os.path.isfile('/etc/urandom') and not os.path.isfile('/etc/random'):
    import time, random, sys
    rng().setSeed(int(time.time() + random.randint(0, rng().maxSeed())) % rng().maxSeed())

# IMPORTANT:
#  the MPI modules work in a master slave mode in that
#  the master intepret the user input and execute the script, it calls
#  the slave when certain operation needs to be done for slaves.
#  the slave accepts commands from master and execute the command.
#  it does not execute the script itself.
#  
#  Therefore, after the module is loaded, the slaves go to a 'feed me'
#  mode

if mpi() and mpiRank() > 0:
    slaveExecutionLoop()
    sys.exit(0)

if not simuOptions['Quiet']:
    if mpi():
        print "simuPOP/MPI : Copyright (c) 2004-2006 Bo Peng"
    else:
        print "simuPOP : Copyright (c) 2004-2006 Bo Peng"
    # compile date, compiler etc are macros that are replaced during compile time.
    print ("Version %s (Revision %d, %s) for Python %s" % (simuVer(), simuRev(), compileDate(),
        compilePyVersion() ))
    print compileCompiler()
    print "Random Number Generator is set to %s with random seed 0x%08x" % (rng().name(), rng().seed())
    # MaxAllele + 1 since 0 is one of the allelic states
    if optimized():
        print "This is the optimied %s allele version with %d maximum allelic states." % (alleleType(), MaxAllele+1)
    else:
        print "This is the standard %s allele version with %d maximum allelic states." % (alleleType(), MaxAllele+1)
    print "For more information, please visit http://simupop.sourceforge.net,"
    print "or email simupop-list@lists.sourceforge.net (subscription required)."

