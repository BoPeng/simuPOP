#!/usr/bin/env python
#
# Testing demographic models
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 149 $
# $LastChangedDate: 2006-02-03 15:51:04 -0600 (Fri, 03 Feb 2006) $
#

import unittest, os, sys
from simuOpt import setOptions
setOptions(quiet=True)
new_argv = []
for arg in sys.argv:
    if arg in ['short', 'long', 'binary', 'mutant', 'lineage']:
        setOptions(alleleType = arg)
    elif arg.startswith('-j'):
        setOptions(numThreads = int(arg[2:]))
    else:
        new_argv.append(arg) 

sys.argv=new_argv
from simuPOP import *
from time import sleep

from simuPOP.demography import *

class TestPlotter(unittest.TestCase):
    def testInstantChangeModel(self):
        demo = InstantChangeModel(T=100, N0=200)
        recordprintDemographicModel(demo)

    
if __name__ == '__main__':
    unittest.main()
