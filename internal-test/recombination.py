#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of mutators  of simuPOP
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# June, 2004.
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
simu = simulator( population(size=10, ploidy=2, loci=[20, 30]),
                  randomMating(), rep=5)

simu.apply( [ initByFreq([.2,.4,.4])])
d = dumper(alleleOnly=1, rep=5)
simu.apply([d])

#  I need one example to show exactly what happens when recombination happens.
# maybe a combination of verboselevel, noMating ?

# recombinator
# help(unifRecombinator)
r = unifRecombinator( rate = 0.1 )
n = noRecombination()

simu.setGen(0)
simu.evolve([r, d], end=20)
simu.setGen(0)
simu.evolve([n, dumper(step=100)], end=100)
