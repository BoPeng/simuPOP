#!/usr/bin/env python
# 
# File: Mating_overlappingGeneration.py
# Author: Bo Peng (bpeng@mdanderson.org)
# 
# Purpose:
#   This scripts demonstrate how to implement overlapping generation
#   and age-structured population, using a non-random mating scheme.
#
#   Individuals in this population has an information field 'age', which
#   is manipulated as follows:
#   
#   1. Individuals are initialized with random ages between 0 and maxAge.
#   2. At each generation, individual ages are increased by 1.
#   3. Individuals with age over maxAge are not involved in mating. (died)
#   4. Individuals not in mating ages are copied to the next generation.
#   5. Individuals in mating ages performs rangom mating and produce
#      offspring with age 0.
#   6. Rpeat from step 2 for desired generations.
#   
# $Date$
# $Revision$
# $HeadURL$
#


from simuPOP import *

import random

def simuOverlappingGeneration(size, maxAge, minMatingAge, maxMatingAge, gen):
    '''
    size         population size.
    maxAge       maximum age. Individuals with age > maxAge will die.
    minMatingAge minimal mating age.
    maxMatingAge maximal mating age.
    gen          generations to simulate
    '''
    pop = population(size, loci=[2], infoFields=['age'])
    pop.setIndInfo([random.randint(0, maxAge) for x in range(size)], 'age')
    # define virtual subpopulations 
    # age < minMatingAge
    # age >= minMatingAge and age < maxMatingAge + 0.1 (age <= maxMatingAge)
    # age >= maxMatingAge + 0.1 and age < maxAge + 0.1 (maxMatingAge < age <= maxAge)
    # age >= maxAge + 0.1 (age > maxAge)
    #
    # Note that we use a cutoff infoSplitter here, it is also possible to 
    # provide a list of values, each corresponding to a virtual subpopulation.
    pop.setVirtualSplitter(infoSplitter('age',
        cutoff=[minMatingAge, maxMatingAge + 0.1, maxAge + 0.1]), 0)
    #
    simu = simulator(pop, heteroMating(
        # age <= maxAge, copy to the next generation (weight=-1)
        [cloneMating(virtualSubPop=x, weight=-1) for x in (0, 1, 2)] +
        # random mating for individuals in mating ages
        [randomMating(virtualSubPop=1)])
    )
    simu.evolve(
        preOps = [initByFreq([0.5, 0.5])],
        ops = [
            # increase age by 1
            infoExec('age += 1', stage=PreMating),
            # count the individuals in each virtual subpopulation
            stat(popSize=True),
            # print virtual subpopulation sizes (there is no individual with age > maxAge after mating)
            pyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in virtualPopSize[0]]))")
        ],
        end = gen - 1
    )

if __name__ == '__main__':
    simuOverlappingGeneration(2000, 10, 4, 6, 100)
