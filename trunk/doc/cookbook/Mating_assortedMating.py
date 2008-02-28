#!/usr/bin/env python
# 
# File: Mating_assortedMating.py
# Author: Bo Peng (bpeng@mdanderson.org)
# 
# Purpose:
#   This scripts demonstrate how to implement assorted mating, namely mating
#   with preference to individuals with similar phenotype.
#
#   The core of this script is a heteroMating mating scheme that use
#   1. general random mating among all individuals
#   2. random mating between individuals with homozygous wildtype allele (0, 0)
#   3. random mating between individuals having at least one mutant (0, 1) or
#       (1,1).
#   
#   A parameter w determines the proportion of offspring produced by these
#   three mating schemes (w, 1, 1).
#
#   During the evolution, number of individuals having genotype (0, 0), (0, 1)
#   and (1, 1) are printed, along with the frequency of allele 0.
#   
#   The simulation confirms that positive assorted mating would lead to the
#   loss of heterozygotes.
#   
# $Date$
# $Revision$
# $HeadURL$
#

import sys
from simuPOP import *
    
def simuAssortedMating(w, size, gen):
    '''
        w       weight of general random mating.
        size    population size
        gen     how many generation to run
    '''
    pop = population(size, loci=[1])
    # define four virtual subpopulations. Individuals in the first three virtual
    # subpopulation have genotype (0, 0), (0, 1) or (1, 0), and (1, 1) respectively,
    # and have at leat one mutant (allele 1) in the last virtual subpopulation.
    pop.setVirtualSplitter(genotypeSplitter(locus=0, 
        alleles=[[0, 0], [0, 1], [1, 1], [0, 1, 1, 1]]), 0)

    # positive weights w, 1, 1 determines the number of offspring produced by each
    # mating scheme.
    simu = simulator(pop, heteroMating([
        randomMating(weight = 0),            # whole population random mating
        randomMating(virtualSubPop=0, weight = 0), # homozygous wildtype
        randomMating(virtualSubPop=3, weight = 0)  # having at least one mutant
        ]))
    #
    simu.evolve(
        preOps = [
            initByFreq([0.5, 0.5]),
            pyExec('AaNum=[]'),  # initialize a list in population's local dictionary
            ],
        ops = [
            # calculate virtual population sizes, and allele frequency at locus 0.
            stat(popSize=True, alleleFreq=[0], stage=PreMating),
            # print size of virtual populations and allele frequency
            pyEval(r"'#inds with genotype AA %4d, Aa %4d, aa %4d, freq of A: %.1f\n' % "
                "(virtualPopSize[0][0], virtualPopSize[0][1], virtualPopSize[0][2],"
                "alleleFreq[0][0]*100)"),
            # append number of individuals with genotype Aa to list AaNum
            pyExec(r"AaNum.append(virtualPopSize[0][1])")
        ],
        end = gen - 1
    )
    return simu.dvars(0).AaNum


if __name__ == '__main__':
    # if there is no parameter
    if len(sys.argv) == 1:
        print simuAssortedMating(1, 2000, 200)
        sys.exit(0)
    # else, draw a figure like what is in Peng2008
    from rpy import *
    replicate = 10
    gen = 20
    popSize = 2000
    r.postscript('assorted.eps')
    r.plot(0, 0, xlim=[0, gen], ylim=[0, popSize], xlab='', ylab='', main='',
        type='n')
    for w in (0, 0.2, 0.4, 0.5, 1, 10):
        Avg_AaNum = [0]*gen
        for rep in range(replicate):
            AaNum = simuAssortedMating(w, popSize, gen)
            Avg_AaNum = [Avg_AaNum[x] + AaNum[x] for x in range(len(AaNum))]
        Avg_AaNum = [x * 1.0 / replicate for x in Avg_AaNum]
        print Avg_AaNum
        r.lines(range(gen), Avg_AaNum)
    r.dev_off()
