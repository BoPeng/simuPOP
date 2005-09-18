#!/usr/bin/env python
#
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
# step 1: load simuPOP
from simuPOP import *
simu = simulator( population(size=10,ploidy=1,loci=[4,3]),
    binomialSelection())

# check fixation
# try one step. 
simu.step( [ alleleCounter([[0]],output=">"), output("\n"),
             terminateIf("any( %alleleNum0 == 0 )")],
             preop = [initByFreq([.2,.4, .4])] )
#
# go! using fixation checker so DO NOT specify end.
# to avoid problem, set end=1000 is a good practice.
simu.setGen(0)
simu.evolve( ops = [
    output("%gen,"),
    alleleCounter(alleles=[[0]],output=">", sep=','),
    output("\n"),
    terminateIf("any( %alleleNum0 == 0 )",output=">>fix.txt")]
    , preop = [initByFreq([.2,.4, .4])],
    end=1000 )

# fix.txt record the generation and fixed alleles of the replicates.
# in this case, there is only one replicate.
print open("fix.txt").read()

# If we do not want to track allele numbers,
simu.setGen(0)
simu.evolve( ops = [
    alleleCounter(alleles=[[0]]),
    terminateIf("any( %alleleNum0 == 0 )",output=">")]
    , preop = [initByFreq([.2,.4, .4])],
    end=1000 )

#
#
# another try, differnt skills
#
# if we want all population start from the same
p = population(size=10,ploidy=1,loci=[4,3])
initByFreq([.2,.4,.4]).apply( p)

simu = simulator( p, binomialSelection(), rep=5)

simu.setGen(0)
# after all replicates reaches fixation, evolve will
# stop.
# fixation checker is classified as terminator for this reason.
simu.evolve( [
    output("%gen,", rep=5),
    alleleCounter(alleles=[[0]], output=">>", sep=',', rep=5),
    output("\b\n",rep=5),
    terminateIf("any( %alleleNum0 == 0 )", output="fix.txt")],
    end=1000 )

print open("fix.txt").read()


# note that for bigger simulation, you can check alleles
# once a while by using
# alleleCounter(...., step=10)

import os
os.remove("fix.txt")

