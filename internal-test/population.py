#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of population 
#  of simuPOP
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# June, 2004.
#
# Note:
#  - individuals can not be created independently since
#    their genotypic info is managed by population
#  - all indexes start from zero. This is because 
#    these functions can be called not only by python,
#    but also by internal implementation.
#  - these functions may look useful but in reality,
#    you will have very few chances to use them directly.
#    Almost everything should be handled by ops.
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
# load module
from simuPOP import *

# create a diploid population with 10 individuals
#
# size: population size (overall size)
# subPop: subPopulation sizes, the size of last subPop can
#         be ignored.
# ploidy: ploidy. 2 means diploid
# loci:   number of loci on each chromosome.
#         [5,7] means two chromosomes with 5, 7 loci respectively
#
turnOnDebug(DBG_POPULATION)

a = population(size=10, ploidy=2, loci=[5, 7], subPop=[2, 8])
# create several ops.
d = dumper()
#
# initialize each locus by allele 1,2,3 with probability .2,.3,.5
initByFreq([.2, .3, .5 ]).apply(a)
#
# dump the population.
# tag     male/female affected/unaffected chroms copy 1 | copy 2 |
# ( 0, 0) MU XXXXXXX|XXXXXXXXXX|
# 
d.apply(a)

# other parameter sets
a = population(  ploidy=2, loci=[5, 7],
               subPop=[2000, 8000])
# other parameter sets
a = population(  ploidy=2, loci=[5, 7],
                 lociDist=[range(5),range(7)],
                 alleleNames=[r"_","A","C","T","G"],
               subPop=[2000, 8000])

#
# print out basic information
print a.popSize()
print a.ploidy()
print a.ploidyName()
print a.numChrom()
#print a.numLoci(5)  # should triger exception.
print a.numLoci(0)
print a.numLoci(1)
# the 11th locus, actually on the second chromosome
print a.locusDist(10)
# to convert between absolute index and relative indices
# use the following functions
print a.absLocusIndex(1,5)
#
# you are supposed to see (1,5) if this is not the case
# report a bug
print a.chromLocusPair(10)
#
print a.locusDist( a.absLocusIndex(1,2) )
#
#
# print subPopulation structure
print a.numSubPop()
print a.totNumLoci()
print a.subPopSize(1)
# again, conversion between absolute and relative indices
print a.absIndIndex(1,1)
#
# you are supposed to see (1,2) if this is not the case
# report a bug
print a.subPopIndPair(4)
#
print a.subPopSize(1)
# dump the population. There will be more info in debug mode than optimized mode.
d.apply(a)


# if you have SNP or DNA sequence as markers
# you can specify allele names
# the first one should always be name for "unknown allele"
a = population(size=10, ploidy=2, loci=[5, 7], subPop=[2, 8],
                alleleNames=['_','A','C','T','G'])
# create several ops.
d = dumper()
initByFreq([.2, .3, .2, .3 ]).apply(a)
d.apply(a)

# you can also specify separator for the output
d.setSeparator(" ")
d.apply(a)



# you can save a population to a file in various
# format and load them.
#

# save it in text format
a.savePopulation("a.txt")
a.savePopulation("a.xml", format="xml")
a.savePopulation("a.bin", format="bin")

# load it in another population
b = loadPopulation("a.txt")
b = loadPopulation("a.xml", format="xml")
b = loadPopulation("a.bin", format="bin")

# populations in a simulator can be saved together
# mating type need to be re-specified during reloading.
#
simu = simulator(a, randomMating(), rep=10)
simu.saveSimulator("s.txt")
simu.saveSimulator("s.xml", format="xml")
simu.saveSimulator("s.bin", format="bin")

# load them
simu1 = loadSimulator("s.txt", randomMating())
simu1 = loadSimulator("s.xml", randomMating(), format="xml")
simu1 = loadSimulator("s.bin", randomMating(), format="bin")

import os
os.remove("a.txt")
os.remove("a.xml")
os.remove("a.bin")
os.remove("s.txt")
os.remove("s.xml")
os.remove("s.bin")

