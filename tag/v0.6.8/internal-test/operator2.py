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
#   2.
#
# load module.
from simuPOP import *
#
# easier ways to set up a simulator and initialize it
simu = simulator(
    population(size=10, ploidy=2, loci=[2, 3 ]), 
    binomialSelection(), rep=5)
#
# initialize populations and display them.
#
# this will not change generation number
#
simu.apply([ initByFreq([.2, .3, .5]),
             dumper(alleleOnly=1, rep=0)])

# Now, I am going to demonstrate how to let ops apply
# on certain generations.
# to print out generation number, we need to use a 'special'
# operator output, which print out generation, replicate
# etc given a printf-like string.
simu.evolve([
    eval("gen", begin=5, end=20, step=3, headers=["gen"]),
    eval(r"'\n'", begin=5, end=20, step=3, rep=REP_LAST)
     ], end=25)

simu.setGen(0)

simu.evolve([ eval(r"str(gen)+'\n'", begin=5, end=20, step=3) ], end=25)

# the above works fine but it print out generation for each
# replicate. (Remember, the ops will be applied to
# each and every populations in a simulator.)
#
# We can use rep parameter to apply output
# only to one replicate.
# But we need first to reset generation counter
simu.setGen(0)
# only output for replicate 5.
endl= eval(r"'\n'", begin=5, end=20, step=3, rep=REP_LAST)
simu.evolve([ eval("str(gen)+'\t'", begin=5, end=20, step=3), endl ], end=25)

# now you can see output is active only at 5,8,11,etc generation.
# but what if I want to output from the last five generations
simu.setGen(0)
simu.evolve([ eval(r"str(gen)+'\n'", begin=-5, end=-1, step=2, rep=1) ], end=25)
#
# active at certain generations --- use of 'at' parameter
# now I also output current replicate
simu.setGen(0)
simu.evolve([eval(r"str(gen)+'\t'+str(rep)+'\n'", at=[1,-1,-2], rep=1) ], end=25)
# NOTE that these parameters are common to *all* ops

# now, what if you want to output something for each replicate
# but an additional return for the last one?
# use two output.
# rep=-1 means last replicate
simu.setGen(0)
simu.evolve([ eval(r"str(gen)+':'+str(rep)+'\t'"),
              eval(r"'\n'",rep=REP_LAST) ], end=25)

# what if you want to split replicate into two groups and
# only output first group?
simu.setGroup([1,1,2,2,2])
