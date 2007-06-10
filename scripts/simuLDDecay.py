#!/usr/bin/env python
#
# Demonstrate the decay of linkage disequilibrium 
#
# Author: Bo Peng (bpeng@rice.edu)
#
# $LastChangedDate$
# $Rev$ 

"""
This program demonstrate the decay of linkage disequilibrium due 
to recombination.
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
     'default':1000, 
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
    {'arg':'r:', 
     'longarg':'recRate=', 
     'default':0.01,
     'label':'Recombination Rate', 
     'allowedTypes':[types.FloatType],
     'description':'Recombination rate',
     'validate':simuOpt.valueBetween(0.,1.),     
    },
    {'arg':'n:', 
     'longarg':'replicate=', 
     'default':5, 
     'label':'Number of Replicate',
     'allowedTypes':[types.IntType, types.LongType],
     'description':'Number of replicates',
     'validate':simuOpt.valueGT(0)
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
    {'longarg':'method=', 
     'default':'D',
     'label':'Choose method',
     'description':'Choose method to compute linkage disequilibrium.',
     'chooseOneOf':['D', "D'", 'R2'],
     'validate': simuOpt.valueOneOf(['D', "D'", 'R2']),
    },
    {'arg':'v', 
     'longarg':'verbose', 
     'default':False, 
     'description':'Verbose mode.'
    },
]


from simuPOP import *

try:
    from simuRPy import *
except:
    print "simuRPy import failed. Please check your rpy installation."
    print "LD values will not be plotted"
    useRPy = False
else:
    useRPy = True

# get all parameters
allParam = simuOpt.getParam(options, __doc__)

if len(allParam) > 0:    # successfully get the params
    (help, popSize, endGen, recRate, numRep, saveFigure, 
    saveConfig, method, verbose) = allParam
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
    print "Recombination rate: ", recRate
    print "Number of replicates: ", numRep
    print "Save figure to: ", saveFigure
    
# diploid population, one chromosome with 2 loci
# random mating with sex
simu = simulator(
    population(size=popSize, ploidy=2, loci=[2]),
    randomMating(),
    rep = numRep)

# get method value used to plot and evolve
if method=="D'":
    methodplot = "LD_prime[0][1]"
    upperlim = 1
    methodeval = r"'%.4f\t' % LD_prime[0][1]"
elif method=='R2':
    methodplot = "R2[0][1]"
    upperlim = 1
    methodeval = r"'%.4f\t' % R2[0][1]"
else:
    methodplot = "LD[0][1]"
    upperlim = 0.25
    methodeval = r"'%.4f\t' % LD[0][1]"


if useRPy:
    plotter = varPlotter(methodplot, numRep=numRep, win=endGen, 
        ylim = [0,upperlim], xlab="generation", saveAs=saveFigure, update=endGen,
        ylab=method, title="Decay of Linkage Disequilibrium r=%.3f" % recRate)
else:
    plotter = noneOp()
    
# everyone will have the same genotype: 01/10
simu.evolve(
    preOps = [initByValue([0,1,1,0])],
    ops = [
        recombinator( rate = recRate),
        stat( alleleFreq=[0], LD=[0,1] ),
        pyEval(methodeval),
        pyEval(r"'\n'", rep=REP_LAST),
        plotter
    ],
    end=endGen
)

# wait five seconds before exit
if useRPy:
    print "Figure will be closed after five seconds."
    time.sleep(5)
