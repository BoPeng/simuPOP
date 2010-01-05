#!/usr/bin/env python

'''
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates how to implement assortative mating, namely mating
  with preference to individuals with similar phenotype.

  The core of this script is a HeteroMating mating scheme that use
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

Date: 2009-01-08
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
    pop = Population(size, loci=[1])
    # define four virtual subpopulations. individuals in the first three virtual
    # subpopulation.have genotype (0, 0), (0, 1) or (1, 0), and (1, 1) respectively,
    # and have at leat one mutant (allele 1) in the last virtual subpopulation.
    pop.setVirtualSplitter(GenotypeSplitter(loci=0,
        alleles=[[0, 0], [0, 1], [1, 1], [0, 0, 0, 1], [0, 1, 1, 1]]))

    pop.evolve(
        initOps = [
            InitSex(),
            InitGenotype(freq=[0.5, 0.5]),
            PyExec('AaNum=[]'),  # initialize a list in population's local dictionary
            ],
            # calculate virtual population sizes, and allele frequency at locus 0.
        preOps = Stat(popSize=True, alleleFreq=[0], subPops=[(0,0), (0,1), (0,2)]),
        # Negative weight means fixed size (weight * current subpopulation size).
        # In the case of no positive weight, zero weights means proportional to
        # parental (virtual) subpopulation size.
        matingScheme = HeteroMating([RandomMating(weight = -1*w),
            RandomMating(subPops=[(0, x) for x in vsp], weight = 0)]),
        postOps = [
            # print size of virtual populations and allele frequency
            PyEval(r"'#inds with genotype AA %4d, Aa %4d, aa %4d, freq of A: %.1f\n' % "
                "(subPopSize[0], subPopSize[1], subPopSize[2], alleleFreq[0][0]*100)"),
            # append number of individuals with genotype Aa to list AaNum
            PyExec(r"AaNum.append(subPopSize[1])")
        ],
        gen = gen
    )
    return pop.dvars().AaNum


if __name__ == '__main__':
    simuAssortativeMating(0.1, 2000, 200)
