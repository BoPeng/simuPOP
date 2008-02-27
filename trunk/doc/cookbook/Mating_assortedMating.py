#!/usr/bin/env python
# 
# Test the evolution of LD, with non-random mating schemes.
#

from simuPOP import *
#from rpy import *
import random
    
# Assuming some individuals to have a preference to mate with another 
# individual with similar phenotype
size = 2000
# weight of random mating
w = 10
#
pop = population(size, loci=[1])
pop.setVirtualSplitter(genotypeSplitter(locus=0, alleles=[[0, 0], [0, 1, 1, 1], [0, 1], [1, 1]]), 0)

simu = simulator(pop, heteroMating([
    randomMating(weight = w), # whole population random mating
    randomMating(virtualSubPop=0, weight = 1),
    randomMating(virtualSubPop=1, weight = 1)
    ]))
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        stat(popSize=True, alleleFreq=[0]),
        pyEval(r"'Individuals with genotype 00: %d, 01: %d, 11: %d, allele frequency of allele 0 %.3f\n' " + \
            "% (virtualPopSize[0][0], virtualPopSize[0][2], virtualPopSize[0][3], alleleFreq[0][0])")
    ],
    end = 2000
)

        
