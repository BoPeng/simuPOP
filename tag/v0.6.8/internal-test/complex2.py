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
pop = population(size=100, ploidy=2, loci=[3,4],subPop=[20,30,50])
d = dumper(alleleOnly=1)
initByFreq([.2,.8]).apply(pop)
#
d.apply(pop)
m = migrator([[0,.2,.1],[.2,0,.1],[.1,.2,0]])
simu = simulator(pop, randomMating(), rep=5)
simu.setGen(0)

simu.evolve([output("%gen\t", output=">>a.txt",rep=1),
             noRecombination(),
             subPopStat([1],output=">>a.txt" ),
             output("\n", output=">>a.txt", rep=5),
             m ] , end=20)

import os
os.remove("a.txt")

