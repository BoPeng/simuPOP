#!/usr/bin/env python
#
# Demonstrate the Hardy-Weinberg equilibrium
#
# Author: Yaji Xu (Yaji.Xu@uth.tmc.edu)
#
# $LastChangedDate: 2007-03-02 14:05:06 -0600 (Fri, 02 Mar 2007) $
# $Rev: 824 $

"""
This program demonstrate the Hardy-weinberg equilibrium when the
allele frequencies in females and males are different.
"""

import os, sys, types, time

import simuOpt
from simuPOP import *

options = [
    {
     'name':'size',
     'default':100000,
     'label':'Population Size',
     'type':[int, long],
     'validator':simuOpt.valueGT(0),
     'description':'''population size. HWE assumes infinite population size
         so large population size improves approximity to theoretical estimates.'''
    },
    {
     'name':'endGen',
     'default':5,
     'type':[int, long],
     'label':'Ending Generation',
     'description':'Length of evolution',
     'validator':simuOpt.valueGT(0)
    },
    {
     'name':'malleleFreq',
     'default':0.4,
     'type':[float, long],
     'label':'Male Allele Frequency',
     'description':'Initial allele frequency in males,',
     'validator':simuOpt.valueBetween(0, 1)
    },
    {
     'name':'falleleFreq',
     'default':0.7,
     'type':[float, long],
     'label':'Female Allele Frequency',
     'description':'Initial allele frequency in females.',
     'validator':simuOpt.valueBetween(0, 1)
    },
]



# get all parameters
par = simuOpt.Params(options, __doc__)
if not par.getParam():
    sys.exit(1)

(popSize, endGen, malleleFreq, falleleFreq) = par.asList()

# diploid population, one chromosome with 1 loci
# random mating with sex
pop = Population(size=popSize, ploidy=2, loci=[1])
pop.setVirtualSplitter(SexSplitter())

# simulation
print "p\tP00 (p^2)\tP01 (2p(1-p))\tP11 ((1-p)^2)"
pop.evolve(
    initOps = [
        InitSex(),
        InitGenotype(freq=[malleleFreq, 1-malleleFreq], subPops=[(0, 0)]),
        InitGenotype(freq=[falleleFreq, 1-falleleFreq], subPops=[(0, 1)])
    ],
    matingScheme = RandomMating(),
    postOps = [
        Stat(alleleFreq=[0], genoFreq=[0]),
        PyEval(r"'%.3f\t%.3f (%.3f)\t%.3f (%.3f)\t%.3f (%.3f)\n' % (alleleFreq[0][0], "\
            "genoFreq[0][(0,0)], alleleFreq[0][0]*alleleFreq[0][0], "\
            "genoFreq[0][(0,1)], 2*alleleFreq[0][0]*(1-alleleFreq[0][0]), "\
            "genoFreq[0][(1,1)], (1-alleleFreq[0][0])*(1-alleleFreq[0][0]) )"),
    ],
    gen=endGen
)


