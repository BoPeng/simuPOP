#!/usr/bin/env python
#
# Demonstrate the Hardy-Weinberg equilibrium 
#
# Author: Yaji Xu (Yaji.Xu@uth.tmc.edu)
#
# $LastChangedDate: 2007-03-02 14:05:06 -0600 (Fri, 02 Mar 2007) $
# $Rev: 824 $ 

"""
This program demonstrate the Hardy-weinberg equilibrium when the  
allele frequencies in females and males are different.
"""

import simuOpt, os, sys, types, time

options = [
    {'arg':'h', 
     'longarg':'help', 
     'default':False, 
     'description':'Print this usage message.', 
     'jump':-1
    },
    {'arg':'s:', 
     'longarg':'size=', 
     'default':100000, 
     'label':'Population Size', 
     'allowedTypes':[types.IntType, types.LongType],
     'validate':simuOpt.valueGT(0),
     'description':'Population size'
    },
    {'arg':'e:', 
     'longarg':'endGen=', 
     'default':50,
     'allowedTypes':[types.IntType, types.LongType],
     'label':'Ending Generation', 
     'description':'Length of evolution',
     'validate':simuOpt.valueGT(0)
    },
    {'arg':'m:', 
     'longarg':'malleleFreq=', 
     'default':0.5,
     'allowedTypes':[types.FloatType, types.LongType],
     'label':'Male Allele Frequency', 
     'description':'Allele Frequency in males',
     'validate':simuOpt.valueBetween(0, 1)
    },
    {'arg':'f:', 
     'longarg':'falleleFreq=', 
     'default':0.5,
     'allowedTypes':[types.FloatType, types.LongType],
     'label':'Female Allele Frequency', 
     'description':'Allele Frequency in females',
     'validate':simuOpt.valueBetween(0, 1)
    },
    {'arg':'s:',
     'longarg':'saveFigure=',
     'label':'Save figure to filename',
     'default':'',
     'allowedTypes':[types.StringType],
     'description':'file the last figure to this filenameXX.eps .'
    },
    {'longarg':'saveConfig=', 
     'default':'', 
     'allowedTypes':[types.StringType],
     'description':'Save current paremeter set to specified file.'
    },
    {'arg':'v', 
     'longarg':'verbose', 
     'default':False, 
     'description':'Verbose mode.'},
    ]


from simuPOP import *

try:
    from simuRPy import *
except:
    print "simuRPy import failed. Please check your rpy installation."
    print "HWE values will not be plotted"
    useRPy = False
else:
    useRPy = True

# get all parameters
allParam = simuOpt.getParam(options, __doc__)

if len(allParam) > 0:    # successfully get the params
    (help, popSize, endGen, malleleFreq, falleleFreq, saveFigure, saveConfig, verbose) = allParam
else:
    sys.exit(0)

if saveConfig != '':
    simuOpt.saveConfig(options, saveConfig, allParam)

if help:
    print simuOpt.usage(options, __doc__)
    sys.exit(1)
    
# print out info if in verbose mode
if verbose:
    print "Pop size: ", popSize
    print "End gen: ", endGen
    print "Save figure to: ", saveFigure
    
# diploid population, one chromosome with 1 loci
# random mating with sex
simu = simulator(
    population(size=popSize, ploidy=2, loci=[1]),
    randomMating(),
    )

#if useRPy:
#    plotter = varPlotter(methodplot, win=endGen, 
#        ylim = [0,1], xlab="generation", saveAs=saveFigure, update=endGen,
#        ylab=method, title="Decay of Linkage Disequilibrium r=%.3f" )
#else:
#    plotter = noneOp()
   
# simulation
simu.evolve(
    preOps = [initByFreq( maleFreq=1, indRange=[0,5000], alleleFreq=[[malleleFreq, 1-malleleFreq]] ),
              initByFreq( maleFreq=0, indRange=[5000,10000], alleleFreq=[[0.5, 0.5]] )
    ],
    ops = [
        stat( alleleFreq=[0], genoFreq=[0] ),
        pyEval(r"'%.3f\t%s\t%s\t%s\n' % (alleleFreq[0][0], genoFreq[0][0][0], genoFreq[0][0][1], genoFreq[0][1][1])"),
    ],
    end=endGen
)

# wait five seconds before exit
if useRPy:
    print "Figure will be closed after five seconds."
    time.sleep(5)
