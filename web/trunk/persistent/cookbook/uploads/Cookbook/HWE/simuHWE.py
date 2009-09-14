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

import simuOpt, os, sys, types, time

options = [
    {'arg':'s:',
     'longarg':'size=',
     'default':100000,
     'label':'Population Size',
     'allowedTypes':[types.IntType, types.LongType],
     'validate':simuOpt.valueGT(0),
     'description':'''Population size. HWE assumes infinite population size
         so large population size improves approximity to theoretical estimates.'''
    },
    {'arg':'e:',
     'longarg':'endGen=',
     'default':5,
     'allowedTypes':[types.IntType, types.LongType],
     'label':'Ending Generation',
     'description':'Length of evolution',
     'validate':simuOpt.valueGT(0)
    },
    {'arg':'m:',
     'longarg':'malleleFreq=',
     'default':0.4,
     'allowedTypes':[types.FloatType, types.LongType],
     'label':'Male Allele Frequency',
     'description':'Initial allele frequency in males,',
     'validate':simuOpt.valueBetween(0, 1)
    },
    {'arg':'f:',
     'longarg':'falleleFreq=',
     'default':0.7,
     'allowedTypes':[types.FloatType, types.LongType],
     'label':'Female Allele Frequency',
     'description':'Initial allele frequency in females.',
     'validate':simuOpt.valueBetween(0, 1)
    },
]


from simuPOP import *

# get all parameters
par = simuOpt.simuOpt(options, __doc__)
if not par.getParam():
    sys.exit(1)

(popSize, endGen, malleleFreq, falleleFreq) = par.asList()

# diploid population, one chromosome with 1 loci
# random mating with sex
pop = population(size=popSize, ploidy=2, loci=[1])
pop.setVirtualSplitter(sexSplitter())
simu = simulator(pop, randomMating())

# simulation
print "p\tP00 (p^2)\tP01 (2p(1-p))\tP11 ((1-p)^2)"
simu.evolve(
    preOps = [
        initSex(),
        initByFreq(alleleFreq=[malleleFreq, 1-malleleFreq], subPops=[(0, 0)], initSex=False),
        initByFreq(alleleFreq=[falleleFreq, 1-falleleFreq], subPops=[(0, 1)], initSex=False)
    ],
    ops = [
        stat(alleleFreq=[0], genoFreq=[0]),
        pyEval(r"'%.3f\t%.3f (%.3f)\t%.3f (%.3f)\t%.3f (%.3f)\n' % (alleleFreq[0][0], "\
            "genoFreq[0][(0,0)], alleleFreq[0][0]*alleleFreq[0][0], "\
            "genoFreq[0][(0,1)], 2*alleleFreq[0][0]*(1-alleleFreq[0][0]), "\
            "genoFreq[0][(1,1)], (1-alleleFreq[0][0])*(1-alleleFreq[0][0]) )"),
    ],
    gen=endGen
)

