#!/usr/bin/env python
#
# Purpose:
#  testing of simuPOP, complex examples
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
pop =  population(size=100, ploidy=2, loci=[10,20,20],
                  subPop=[20, 30, 50])

initByFreq([.2,.8]).apply(pop)

# migration
m = migrator([[0,.2,.1],[.2,0,.1],[.1,.2,0]], mode=MigrByProbability)
# recombination
r = unifRecombinator(rate=.1)
# mutation
mu = kamMutator(rate=(0.01))

simu = simulator(pop, randomMating(),rep=5)
#
simu.setGen(0)

# this seems to be running OK.
# note that subPopStat by default does not output
simu.evolve( [# m,r,mu,
              output("%gen\t", rep=1), 
              subPopStat([1,2],output=">",rep=1),
              output("\n",rep=1)],
            end=20)

# check variables
listVars(1)
