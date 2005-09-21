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
simu.step( [ stat([[0,1,2,3]],output=">"), output("\n"),
             terminateIf(r"min( alleleNum0) == 0")],
             preOps = [initByFreq([.2,.4, .4])] )
#
# go! using fixation checker so DO NOT specify end.
# to avoid problem, set end=1000 is a good practice.
simu.setGen(0)
simu.evolve( ops = [
    pyEval(r"str(gen)+','"),
    stat(alleleFreq=[0[0,1,2,3]],output=">", sep=','),
    output("\n"),
    terminateIf("min( alleleNum0) == 0",output=">>fix.txt")]
    , preOps = [initByFreq([.2,.4, .4])],
    end=1000 )

# fix.txt record the generation and fixed alleles of the replicates.
# in this case, there is only one replicate.
print open("fix.txt").read()

# If we do not want to track allele numbers,
simu.setGen(0)
simu.evolve( ops = [
    stat(alleleFreq=[0[0,1,2,3]]),
    terminateIf(r"min(alleleNum0) == 0",output=">") ]
    , preOps= [initByFreq([.2,.4, .4])],
    end=1000 )

#
#
# another try, differnt skills
#
# if we want all population start from the same
p = population(size=10,ploidy=1,loci=[4,3])
initByFreq([.2,.4,.4]).apply( p)

simu = simulator( p, binomialSelection(), rep=5)


