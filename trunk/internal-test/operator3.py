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
## Important:  set path to simuPOP.py
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
# easier ways to set up a simulator and initialize it
simu = simulator(
    population(size=10, ploidy=2, loci=[2, 3 ]), 
    randomMating(), rep=5)
# this will not change generation number
simu.apply([ initByFreq([.2, .3, .5]),
             dumper(alleleOnly=1)])

 
# specify output is kind of complicated (because we sometimes need to 
# let ops output to the same file.)
#
# here is the general rule: (details please refer to online documents.)
#
#   ">file.txt"  output to this file, once done, close.
#                if mutliple ops output to the same file
#                only the last one will succeed.
#   ">>file.txt" output to this file, append, close when evolve stop.
#   ">>file%rep.txt" different filename for each replicates.
#                %rep will be substituted by replicate number.
#   ">>file%gen.txt" different filename for different generation.
#                %grp will be substituted by generation number
#   ">>%{sim}%rep.txt variables (set by setStringVar etc) can be used.
#   "|pipeName"   A pipe can be used. I.e., everything will be
#                written to a memory file. This is useful only
#                when the memory file will be read by a
#                DataProvider.
#
# for example:
#
simu.setGen(0)
# note the use of sep
# count allele at each locus
# header format: replicate_loci_allele
simu.evolve([ alleleCounter(alleles=[[0]],
                            output=">", sep=",",rep=1),
              eval(r'"\n"', rep=REP_LAST)
           ], end=10)

#  if you would like to supress the last comma, use a
# \b backspace character. However, this is useful
# only when the output is standard output
# (If you write to a file, ',' and \b will be written)
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]],output=">>", sep=',',rep=2),
              eval(r"'\n'", rep=-1)
           ], end=10)

# output to another file
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]],output="a.txt"),
              eval(r'"\n"', rep=-1, output="a.txt")
           ], end=10)

# Nothing!! (this usage is totally python)
print open("a.txt").read()

# you should use >> so that output will
# append instead of replace output from alleleCounter
# In the previous case, only the last output write to
# a.txt, i.e., only a \n was written.
#
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]],
              output=">>a.txt", rep=2),
              eval(r'"\n"', rep=REP_LAST, output=">>a.txt")
           ], end=10)


print open("a.txt").read()

# note that once you specify "append mode", you do not have
# to do this each time.
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]],output=">>a.txt", sep=','),
              eval(r"'\b\n'", rep=REP_LAST, output="a.txt")
           ], end=10)

print open("a.txt").read()

# replicate specific output
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]],
                            outputExpr="'>>a' + str(rep) + '.txt'", sep=','),
              eval(r"'\b\n'", outputExpr="'a'+str(rep)+'.txt'")
           ], end=10)

print open("a0.txt").read()
print open("a2.txt").read()

# generation specific
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]],
        outputExpr="'>>b'+str(gen) + '.txt'",sep=',')
           ], end=10)

# note that no header is output.
#print open("b1.txt").read()
#print open("b2.txt").read()

# you can also specify "no output" 
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]],output="")], end=10)

# you can actually combine replication, generation, group specific
# outputs....
# to avoid generting too many files, no example is given here.

# remove all generated files

import os
os.remove("a.txt")
for n in range(1,6):
    os.remove( 'a%d.txt' % n)
#for n in range(0,11):
#    os.remove( 'b%d.txt' % n)



# advanced usage: shared variables
#
# ops can save/load sharedvariables to a gloval
# shared variable area. This is one way (another one is
# through readable files) to share information through
# ops.
#
# Details please refer to misc.py
#
#
#
# area=0 is the public, default area.
#
name = "simulation1"

# we can use variable names
# now, in operator,
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]],
    outputExpr=r"'>>'+name+'.txt'", sep=','),
  eval(r"'\b\n'", rep=REP_LAST,
       outputExpr="name+'.txt'")
    ], end=10)

# we can not use this name in python though.
print open("simulation1.txt").read()
os.remove("simulation1.txt")

