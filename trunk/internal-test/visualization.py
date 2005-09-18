#!/usr/bin/env python
#
# Purpose:
#  testing of simuPOP, complex examples
#
# Author
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
#import sys
#sys.path.append('../lib/debug')
# 
from simuPOP import *
pop = population(size=200, ploidy=2, loci=[3,4],subPop=[50,50,100])
simu = simulator(pop, randomMating(), rep=5)

d = dumper(alleleOnly=1)
simu.apply( [ initByFreq([.2,.8])])
migr = migrator([[0,.2,.1],[.2,0,.1],[.1,.2,0]])

#turnOnDebug(DBG_ALL)

# first let us see what the data look like
mem=">"
simu.setGen(0)
simu.evolve([
    output("%gen\t",output=mem,rep=1),
    subPopStat([0], output=mem, sep="\t", rep=3 ),
    output("\n", output=mem, rep=5),
    migr
    ] , end=20)

# then go to a file?
mem=">>File.txt"
simu.setGen(0)
simu.evolve([
    output("%gen\t",output=mem,rep=1),
    subPopStat([0], output=mem, sep="\t" ),
    output("\n", output=mem, rep=5),
    migr
    ] , end=20)

print open("File.txt").read()

# what if we want to add to File.txt
mem=">>>File.txt"
simu.setGen(0)
simu.evolve([
    output("%gen\t",output=mem,rep=1),
    subPopStat([0], output=mem, sep="\t" ),
    output("\n", output=mem, rep=5),
    migr
    ] , end=20)

#
# pick up from this file?
# Note that if the file is in >>> mode, the datasource
# might be in trouble since there may be two outputs
# (this is usually not a problem since DataSource ONLY
# provide data from the last line.


print open("File.txt").read()

import os
os.remove("File.txt")

# we can also use a pipe
mem="|mem"
simu.setGen(0)
simu.evolve([
    output("%gen\t",output=mem,rep=1),
    subPopStat([0], output=mem, sep="\t" ),
    output("\n", output=mem, rep=5),
    migr
    ] , end=20)


# Now, let us do visualization:
#
#first, have a look at data
simu.setGen(0)
simu.evolve([
    output("%gen\t",rep=1),
             subPopStat([0], output=">",rep=1 ),
             output("\n", rep=5),
             migr ] , end=20)

# data source?
mem="|memFile"
simu.setGen(0)
simu.evolve([
    output("%gen\t",output=mem,rep=1),
    subPopStat([0], output=mem,rep=1 ),
    output("\n", output=mem, rep=5),
    migr,
    outputDataSource(fileDataSource(source=mem),rep=5),
    output("\n",rep=5)
  ] , end=20)

#
# visualize this datasource
#
# only be called at fifth replicate
visual = matlabPlotter(fileDataSource(mem), win=50, rep=5, update=5)
simu.setGen(0)
simu.evolve([
    output("%gen\t",output=mem,rep=1),
    subPopStat([0], output=mem,rep=1 ),
    output("\n", output=mem, rep=5),
    migr,
    visual
  ] , end=20)

# you can do what ever you want through
# most likely, you would like to save the image.
#visual.evalString("data_1=[]")
#visual.evalString("clf")

#
#
# the second method is easier
# we make use the fact that subPopStat write shared
# array subPopSize
#
# let us see
listVars(1)
#
# subPopsize is there for each replicate.
#
# first view datasource
simu.setGen(0)
simu.evolve([
  subPopStat([0], output="" ),
  migr,
  outputDataSource(
    # datasource is for one replicate only
    varDataSource(preSource="%gen",source="%subPopSize",rep=1), 
    rep=5), # only call varDataSource for one replicate
  output("\n",rep=5)
  ] , end=20)


# now visualize

simu.setGen(0)
simu.evolve([
  subPopStat([0], output="" ),
  migr,
  matlabPlotter(
    varDataSource(preSource="%gen",source="%subPopSize",rep=1), 
  win=50, rep=5)
  ] , end=20)
