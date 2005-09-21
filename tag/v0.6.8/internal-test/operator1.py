#!/usr/bin/env python
#
# Purpose:
#   testing operator behaviors for simupoop.
# 
# Author:
#   Bo Peng (bpeng@rice.edu)
#
# Note:
#   1. %gen -- generation
#      %rep -- replicate
#      %grp -- group
#
# Important:  set path to simuPOP.py
#
#  if you do not set PYTHONPATH, you can uncomment
#  the following two lines and specify path to
#  simuPOP.py here.
#
#import sys
#sys.path.append('../lib/debug')
# 
# load module.
from simuPOP import *
# 
pop = population(size=10, ploidy=2, loci=[2, 3 ], subPop=[])
initByFreq([.2, .3, .5], maleFreq=0.8).apply(pop)

# here is an operator we have seen
d = dumper()
d.apply(pop)

# what are these ops?
# baseOperator (Operator in C++) defines the common interface
# for all ops
#help(baseOperator)
#help(dumper)
print d.canApplyPreMating()
print d.canApplyDuringMating()
# so dumper is a post mating operator
print d.canApplyPostMating()
# apply to any replicate
print d.applicableReplicate()
# apply to any group of replicates
print d.applicableGroup()
# always on... if current is 1 (we do not know current

# operator can take its own options (this option is unique to dumper)
d1 = dumper(alleleOnly=1)
# note that the output might be different in optimized mode
d1.apply(pop)
# most arguments can be inquired and reset. The name convention is
# argumentName() and setArgumentName()
print d1.alleleOnly()
d1.setAlleleOnly(0)
d1.apply(pop)
#
#
# Opertors can be used individually (especially
# initializer and outputer but most of them
# are used along with a simulator.
#
# this simulator controls two replicates of pop
# and will use randomMating() to generate next
# generation.
simu = simulator(pop, randomMating(), rep=2)
#
# have a look at properties of simu.
#help(simulator)
# current generation?
print simu.gen()
#
# initialy, groups are set to replicate index
# it can be set by setGroups().
# this usage will be shown in operator2.py
print simu.group()
simu.setGroup([2,3])
print simu.group()
print _vars
#
# get a reference to the first replicate
p1 = simu.population(1)
# have a look
d.apply(p1)
# NOTE that p1 is now the same as pop
# but will be different after simulator
# evolves.


# simu can apply an array of pre, post mating ops
# without evolution.
# in this case, every replicate will be dumped.
simu.apply([d])
# generation is not changed
print simu.gen()
#
# you can evolve step by step, each time print all population
simu.step([d])
simu.step([d])
# or, you can only look at the second replicate (details later)
simu.step([dumper(rep=2)])
# evolve using random mating and NO additional ops
simu.evolve([ ], end=25)
print simu.gen()
