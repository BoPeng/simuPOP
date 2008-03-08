#!/usr/bin/env python

'''
File: Mating_assortativeMating.py
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates how to implement assortative mating, namely mating
  with preference to individuals with similar phenotype.

  The core of this script is a heteroMating mating scheme that use
  1. general random mating among all individuals
  2. random mating between individuals with homozygous wildtype allele (0, 0)
  3. random mating between individuals having at least one mutant (0, 1) or
      (1,1).
  
  A parameter w determines the proportion of offspring produced by the general
  random mating scheme. w = 1 means no assortative mating. The numbers of offspring
  produced by other two mating schemes are proportional to the size of the
  corresponding virtual subpopulations in the parental generation. For example,
  if the population size is 1000, w=0.5, and there are 200 homozygous wildtype
  individuals, the number of offspring produced by these three mating schemes
  will be 500, 100 and 400. The size of virtual subpopulations will change
  as a result of general random mating.

  During the evolution, number of individuals having genotype (0, 0), (0, 1)
  and (1, 1) are printed, along with the frequency of allele 0.
  
  The simulation confirms that positive assortative mating would lead to the
  loss of heterozygotes.
  
$Date$
$Revision$
$HeadURL$
'''

import sys
from simuPOP import *
    
def simuAssortativeMating(w, size, gen, vsp=[0, 4]):
    '''
        w       proportion of general random mating.
        size    population size
        gen     how many generation to run
        vsp     virtual subpopulations for assortative mating.
    '''
    pop = population(size, loci=[1])
    # define four virtual subpopulations. Individuals in the first three virtual
    # subpopulation have genotype (0, 0), (0, 1) or (1, 0), and (1, 1) respectively,
    # and have at leat one mutant (allele 1) in the last virtual subpopulation.
    pop.setVirtualSplitter(genotypeSplitter(locus=0, 
        alleles=[[0, 0], [0, 1], [1, 1], [0, 0, 0, 1], [0, 1, 1, 1]]), 0)

    # Negative weight means fixed size (weight * current subpopulation size).
    # In the case of no positive weight, zero weights means proportional to
    # parental (virtual) subpopulation size.
    simu = simulator(pop, heteroMating([randomMating(weight = -1*w)] + \
        [randomMating(virtualSubPop=x, weight = 0) for x in vsp]))
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
    simuAssortativeMating(0.1, 2000, 200)
