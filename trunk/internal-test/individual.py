#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of individuals
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
#    these interfaces can be used not only by python,
#    but also by internal implementation which uses
#    zero index arrays. Fortunately, these interfaces
#    are not supposed to be used directly a lot.
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

# create a population, intialize it.
pop = population(10, ploidy=2, loci=[5, 7], subPop=[2,8])
# create several ops.
d = dumper()
# sexFreq is by default
initByFreq(alleleFreq=[.2, .3, .5 ], maleFreq=0.5).apply(pop)

d.apply(pop)
#
# get individual info
ind = pop.individual(9)
# genotypic information can be retrieved from
# individuals the same way as from population.
print ind.ploidy()
print ind.ploidyName()
print ind.numChrom()
print ind.numLoci(0)
print ind.numLoci(1)
print ind.locusDist(10)
print ind.totNumLoci()

# for a list of functions of individual, run
# help(individual)
#
# individualized information
print ind.affected()
print ind.affectedChar()
ind.setAffected(1)
print ind.affectedChar()
# compare with dumper output, see the meaning of M/F,U/A .
d.apply(pop)
print ind.sex()
print ind.sexChar()
#
#
# How to use info of individuals?
#
# The following handles migration by hand.
# This is *not* the way to work with simuPOP, but a
# demonstration of features.
#
# Operators will handle all things like this.
#
pop.setIndInfoWithSubPopID()
# info is not shown
d.apply(pop)
# if we turn on debug 
turnOnDebug(DBG_INDIVIDUAL)
d.apply(pop)
print pop.individual(5).info()
pop.individual(0).setInfo(1)
pop.individual(4).setInfo(0)
pop.individual(2).setInfo(4)
d.apply(pop)
pop.rearrangeByIndInfo()
d.apply(pop)
pop.setSubPopByIndInfo()
# adjust for shallow copying 
# this is done automatically if within subPop
# allele iterators are requested.
# some ops requires all alleles of a subPopulation
# within its allelic boundary ...
pop.adjustGenoPosition()
d.apply(pop)
#
