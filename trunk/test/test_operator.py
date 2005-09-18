#!/usr/bin/env python
#
# Purpose:
#   testing operator behaviors for simupoop.
# 
# Author:
#   Bo Peng (bpeng@rice.edu)
#
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
print simu.grp()
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
simu.apply([ initByFreq([.2, .3, .5]), dumper(alleleOnly=1, rep=1)])

# Now, I am going to demonstrate how to let ops apply
# on certain generations.
# to print out generation number, we need to use a 'special'
# operator output, which print out generation, replicate
# etc given a printf-like string.
simu.evolve([ pyEval(r"str(gen)+'\n'", begin=5, end=20, step=3) ],
            end=25)

simu.setGen(0)
simu.evolve([ pyEval(r"str(gen)+'\n'",end=2) ],
          end=25)

simu.setGen(0)
simu.evolve([ pyEval(r"str(gen)+'\n'",end=20, step=5) ],
            end=25)
# the above works fine but it print out generation for each
# replicate. (Remember, the ops will be applied to
# each and every populations in a simulator.)
#
# We can use rep parameter to apply output
# only to one replicate.
# But we need first to reset generation counter
simu.setGen(0)
# only output for replicate 5.
simu.evolve([pyEval(r"str(gen)+'\n'" , begin=5, end=20, step=3, rep=1) ], end=25)

# now you can see output is active only at 5,8,11,etc generation.
# but what if I want to output from the last five generations
simu.setGen(0)
simu.evolve([pyEval(r"str(gen)+'\n'" , begin=-5, end=-1, step=2, rep=1) ], end=25)
#
# active at certain generations --- use of 'at' parameter
# now I also output current replicate
simu.setGen(0)
simu.evolve([pyEval(r"str(gen)+'\t'+str(rep)+'\n'", at=[1,-1,-2], rep=1) ], end=25)
# NOTE that these parameters are common to *all* ops

# now, what if you want to output something for each replicate
# but an additional return for the last one?
# use two output.
# rep=-1 means last replicate
simu.setGen(0)
simu.evolve([ pyEval(r"str(gen)+'\t'+str(rep)"), output("\n",rep=REP_LAST)
              ], end=25)

# what if you want to split replicate into two groups and
# only output first group?
simu.setGroup([1,1,2,2,2])
simu.setGen(0)
simu.evolve([ pyEval(r"'%d:%d:%d\t'%(gen,rep,grp)", grp=2),
              output("\n",rep=REP_LAST) ], end=25)
#
# this feature is useful when you want to initialize
# or evolve two sets of populations.
#
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
simu.evolve([ stat(alleleFreq=[2],
                            output=">", sep=",",rep=1),
              output("\n", rep=REP_LAST)
           ], end=10)

#  if you would like to supress the last comma, use a
# \b backspace character. However, this is useful
# only when the output is standard output
# (If you write to a file, ',' and \b will be written)
simu.setGen(0)
simu.evolve([ stat(alleleFreq=[0],output=">>", sep=',',rep=2),
              output("\b\n", rep=REP_LAST)
           ], end=10)

# output to another file
simu.setGen(0)
simu.evolve([ stat(alleleFreq=[0],output="a.txt"),
              output("\n", rep=REP_LAST, output="a.txt")
           ], end=10)

# Nothing!! (this usage is totally python)
print open("a.txt").read()

# you should use >> so that output will
# append instead of replace output from stat
# In the previous case, only the last output write to
# a.txt, i.e., only a \n was written.
#
simu.setGen(0)
simu.evolve([ stat(alleleFreq=[0],
              output=">>a.txt", rep=2),
              output("\n", rep=REP_LAST, output=">>a.txt")
           ], end=10)


print open("a.txt").read()

# note that once you specify "append mode", you do not have
# to do this each time.
simu.setGen(0)
simu.evolve([ stat(alleleFreq=[0],output=">>a.txt", sep=','),
              output("\b\n", rep=REP_LAST, output="a.txt")
           ], end=10)

print open("a.txt").read()

# replicate specific output
simu.setGen(0)
simu.evolve([ stat(alleleFreq=[0],
   outputExpr=r"'>>a'+str(rep)+'.txt'",
             sep=','),
   output("\b\n", outputExpr=r"'>>a'+str(rep)+'.txt'")
  ], end=10)

print open("a1.txt").read()
print open("a2.txt").read()

# generation specific
simu.setGen(0)
simu.evolve([ stat(alleleFreq=[0],
                            outputExpr=r"'>>b'+str(rep)+'.txt'" ,sep=',')
           ], end=10)

# note that no header is output.
#print open("b1.txt").read()
#print open("b2.txt").read()

# you can also specify "no output" 
simu.setGen(0)
simu.evolve([ stat(alleleFreq=[0],output="")], end=10)

# you can actually combine replication, generation, group specific
# outputs....
# to avoid generting too many files, no example is given here.

# remove all generated files

import os
os.remove("a.txt")
for n in range(5):
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
