#!/usr/bin/env python

'''
File: Mating_consanguineousMating.py
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates how to implement consanguineous mating, namely mating
  with preference to relatives of individuals.

  The core of this script is a heteroMating mating scheme that use
  1. general random mating among all individuals
  2. consanguineous mating between first cousins
  I use two mating schemes because it is unlikely that all matings happens
  between relatives in a population.

  A parameter w determines the proportion of offspring produced by the general
  random mating scheme. w = 1 means no consanguineous mating. The numbers of offspring
  produced by other two mating schemes are proportional to the size of the
  corresponding virtual subpopulations in the parental generation.

  The most difficult part of such a mating scheme is to find relatives of each
  individual. Using a parentsTagger(), simuPOP would record the indexes of each
  individual's parents in information fields 'father_idx' and 'mother_idx'. However,
  going through the pedigree to locate relative can be tedious. Fortunately, simuPOP
  provides two functions that can hopefully reduce the troubles.

  Some minor points:
  1. A population should keep at least two parental generations to locate first cousin,
    and more generations to locate more distance relatives. (c.f. parameter ancGen
    of the population() function).
  2. If only one offspring is produced at each mating event, it would be extremely
    difficult to find full siblings, and cousins in a population. Use numOffspring = 2
    or more advanced offspring control schemes to increace relatives in a population.
  3. For the first several generation where there are not enough ancestral generations,
    no relatives could be found and consanguineousMating would locate another parent
    from the whole population.
  4. In this example, both mating schemes are applied to the whole population so there
    is no need to define virtual subpopulations.
  5. Number of sibling and offspring vary from individual to individual due
    to chances. It is therefore difficult to estimate how many information fields
    are needed to store all siblings and offsprings. Big number is likely to
    slow down your simulation, and small number may limit the number of relatives
    you can locate.

$Date: 2008-06-15 15:57:22 -0500 (Sun, 15 Jun 2008) $
$Revision: 1624 $
$HeadURL: https://simupop.svn.sourceforge.net/svnroot/simupop/trunk/doc/cookbook/Mating_consanguineousMating.py $
'''

import sys
from simuPOP import *

def findCousin(pop, fields):
    'Find cousins of each individual and put their indexes to cousinFields'
    # locate offsprings
    if pop.ancestralGens() < 1:
        return
    (parFields, sibFields, offFields, cousinFields) = fields
    # create a pedigree object from the population
    ped = pedigree(pop, infoFields=parFields + cousinFields)
    # add intermediate information fields
    ped.addInfoFields(sibFields + offFields)
    ped.locateRelatives(Offspring, offFields);
    ped.locateRelatives(FullSibling, sibFields);
    # Find parents -> siblings -> offspring
    # Another parameter pathSex can control the sex of each step.
    ped.traceRelatives(pathGen = [0, 1, 1, 0],
        pathFields = [parFields, sibFields, offFields],
        resultFields = cousinFields)
    # get indexes of cousins from the pedigree object
    pop.updateInfoFieldsFrom(cousinFields, ped)

def simuConsanguineousMating(w, size, gen, numFields=4):
    '''
        w       proportion of general random mating.
        size    population size
        gen     how many generation to run
    '''
    parFields = ['father_idx', 'mother_idx']
    sibFields = ['sibling%d' % x for x in range(numFields)]
    offFields = ['offspring%d' % x for x in range(numFields)]
    cousinFields = ['cousin%d' % x for x in range(numFields)]

    pop = population(size, loci=[1], ancGen=2, infoFields = parFields + cousinFields)

    simu = simulator(pop, heteroMating([
        randomMating(numOffspring=2, weight = w),
        consanguineousMating(infoFields = cousinFields, func=findCousin,
            param = [parFields, sibFields, offFields, cousinFields],
            numOffspring = 2, weight = 1 - w)
        ])
    )
    #
    simu.evolve(
        preOps = [ initByFreq([0.5, 0.5]) ],
        ops = [
            parentsTagger(),
            pyEval(r'"%d\n" % gen'),
        ],
        gen = gen
    )
    return True


if __name__ == '__main__':
    simuConsanguineousMating(0.5, 2000, 20, 4)
