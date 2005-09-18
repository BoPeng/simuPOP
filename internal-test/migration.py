#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of mutators  of simuPOP
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# June, 2004.
#
# Important:  set path to simuPOP.py
#
#  if you do not set PYTHONPATH, you can uncomment
#  the following two lines and specify path to
#  simuPOP.py here.
#

from simuPOP import *
pop = population(size=200, ploidy=2, loci=[3,4],subPop=[50,50,100])
d = dumper(alleleOnly=1)
initByFreq([.2, .8]).apply(pop)
d.apply(pop)

#migrator by probability
m = migrator([[0,.2,.1],[.2,0,.1],[.1,.2,0]], mode=MigrByProbability)
# mode1 will show accumulative rates
print m.rates()
simu = simulator(pop, randomMating(),rep=5)
simu.setGen(0)
simu.evolve([
           output("%gen\t", rep=1),
           m, subPopStat([1],output=">"), output("\t"),
           output("\n",rep=5)
    ],end=5)


# second mode of migration: by proportion
m = migrator([[0,.2,.1],[.2,0,.1],[.1,.2,0]], mode=MigrByProportion)
print m.rates()
simu = simulator(pop, randomMating(), rep=5)
simu.setGen(0)
simu.evolve([output("%gen\t",rep=1),
             subPopStat([1]), output("\n", rep=5), m ] , end=5)

#
m = migrator([[0, 20, 20], [18, 0, 20], [20, 20, 0]], mode=MigrByCounts)
print m.rates()

del simu


# create new subPop
#
#
simu = simulator(pop, randomMating(), rep=5)
m = migrator(rates=[[ .3]], fromSubPop=[0], toSubPop=[4])
simu.setGen(0)
simu.evolve([output("%gen\t",rep=1),
             subPopStat([0], output=">"), output("\t"),
             output("\n", rep=5), m ] , end=2)

# 
