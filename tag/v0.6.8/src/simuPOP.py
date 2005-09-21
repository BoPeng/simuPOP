#!/usr/bin/env python

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
  showSimuPopInfo()
