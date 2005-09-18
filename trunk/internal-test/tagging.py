#!/usr/bin/env python
# Important:  set path to simuPOP.py
#
#  if you do not set PYTHONPATH, you can uncomment
#  the following two lines and specify path to
#  simuPOP.py here.
#
#import sys
#sys.path.append('../lib/debug')
# 
from simuPOP import *
simu = simulator( population(size=20, ploidy=2, loci=[2,4], subPop=[5,15]),
         randomMating())
simu.apply( [initByFreq([.2,.3,.5])])
#
# randomMating with tag of parents
#simu.setGen(0)
simu.evolve([parentsTagger(), dumper(alleleOnly=1)], end=3)


