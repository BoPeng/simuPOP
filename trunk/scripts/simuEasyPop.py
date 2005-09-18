#!/usr/bin/env python
#
# repeat part of the functions of easyPOP. Hopefully easier to use
#
# Author: Bo Peng (bpeng@rice.edu)
#
# Date:   Apr. 2005
# Last Updated: Apr. 2005
# 

# 
# OPTIONS:
#

"""
This program 
"""


import simuOpt, os, sys, types

print "Not finished yet."

sys.exit(0)

options = [
  {'arg':'h', 'longarg':'help', 'default':False, 
   'description':'Print this usage message.', 'jump':-1},
  {'arg':'s:', 'longarg':'size=', 'default':1000, 
   'configName':'popSize', 'allowedTypes':[types.IntType, types.LongType],
   'prompt':'Population size (1000): ',
   'description':'Population size'},
  {'arg':'e:', 'longarg':'endGen=', 'default':50,
   'allowedTypes':[types.IntType, types.LongType],
   'configName':'endGen', 'prompt':'Length of evolution (50): ',
   'description':'Length of evolution'},
  {'arg':'r:', 'longarg':'recRate=', 'default':0.01,
   'configName':'recombinationRate', 'allowedTypes':[types.FloatType],
   'prompt': 'Recombination rate (0.01): ', 
   'description':'Recombination rate'},
  {'arg':'n:', 'longarg':'replicate=', 'default':5, 
   'configName':'numOfReplicate',
   'allowedTypes':[types.IntType, types.LongType],
   'prompt':'Number of replicates (5): ',
   'description':'Number of replicates'},
  {'longarg':'saveConfig=', 'default':None, 'allowedTypes':[types.StringType],
   'description':'Save current paremeter set to specified file.'},
  {'arg':'v', 'longarg':'verbose', 'default':False, 
   'description':'Verbose mode.'},
  ]


from simuPOP import *
from simuUtil import *
from simuRPy import *

# -h will be handled before anyone else
(help, popSize, endGen, recRate, numRep, saveConfig, verbose) = \
  simuOpt.getParam(options)

if saveConfig != None:
  simuOpt.saveConfig(options, saveConfig, 
    (help, popSize, endGen, recRate, numRep, saveConfig, verbose))

if help:
  print simuOpt.usage(options, __doc__)
  sys.exit(1)
  
# print out info if in verbose mode
if verbose:
  print "Pop size: ", popSize
  print "End gen: ", endGen
  print "Recombination rate: ", recRate
  print "Number of replicates: ", numRep
  
