#!/usr/bin/env python

'''
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates how to implement overlapping generation
  and age-structured population, using a non-random mating scheme.

  individuals in this population.has an information field 'age', which
  is manipulated as follows:

  1. individuals are initialized with random ages between 0 and maxAge.
  2. At each generation, Individual ages are increased by 1.
  3. individuals with age over maxAge are not involved in mating. (died)
  4. individuals not in mating ages are copied to the next generation.
  5. individuals in mating ages performs rangom mating and produce
     offspring with age 0.
  6. Rpeat from step 2 for desired generations.

Date: 2008-12-13
'''


from simuPOP import *

import random

def simuOverlappingGeneration(size, maxAge, minMatingAge, maxMatingAge, gen):
    '''
    size         population size.
    maxAge       maximum age. individuals with age > maxAge will die.
    minMatingAge minimal mating age.
    maxMatingAge maximal mating age.
    gen          generations to simulate
    '''
    pop = Population(size, loci=[2], infoFields=['age'])
    pop.setIndInfo([random.randint(0, maxAge) for x in range(size)], 'age')
    # define virtual subpopulations
    # age < minMatingAge
    # age >= minMatingAge and age < maxMatingAge + 0.1 (age <= maxMatingAge)
    # age >= maxMatingAge + 0.1 and age < maxAge + 0.1 (maxMatingAge < age <= maxAge)
    # age >= maxAge + 0.1 (age > maxAge)
    #
    # Note that we use a cutoff InfoSplitter here, it is also possible to
    # provide a list of values, each corresponding to a virtual subpopulation.
    pop.setVirtualSplitter(InfoSplitter('age',
        cutoff=[minMatingAge, maxMatingAge + 0.1, maxAge + 0.1]))
    #
    pop.evolve(
        initOps = [
            InitSex(),
            InitGenotype(freq=[0.5, 0.5])
        ],
            # increase age by 1
        preOps = InfoExec('age += 1'),
        matingScheme = HeteroMating(
            # age <= maxAge, copy to the next generation (weight=-1)
            [CloneMating(subPops=[(0, x) for x in (0, 1, 2)], weight=-1),
            # random mating for individuals in mating ages
            RandomMating(subPops=[(0, 1)])]),
        postOps = [
            # count the individuals in each virtual subpopulation
            Stat(popSize=True, subPops=[(0,0), (0,1), (0,2), (0,3)]),
            # print virtual subpopulation sizes (there is no individual with age > maxAge after mating)
            PyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in subPopSize]))")
        ],
        gen = gen
    )

if __name__ == '__main__':
    simuOverlappingGeneration(2000, 10, 4, 6, 100)
