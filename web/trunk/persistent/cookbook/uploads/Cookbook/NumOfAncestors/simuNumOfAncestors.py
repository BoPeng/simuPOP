#!/usr/bin/env python
#
# Count the number of ancestors to a DNA sequences
#
# Author: Haiyan Teng
#
# $LastChangedDate: 2008-10-15 13:10:09 -0500 (Wed, 15 Oct 2008) $

"""
This program counts the average number of ancestors who contribute
genotype to an offspring after a number of generations. The script
creates a population and initialize all alleles of an individual with
their indexes.

Ind 0: 000000000000000 000000000000000
Ind 1: 111111111111111 111111111111111
...
Ind x: xxxxxxxxxxxxxxx xxxxxxxxxxxxxxx

The population is evolved for a number of generations, following a
Mendelian inheritance model. Recombinations are applied so that offspring
can inherit genotype from both parental homologous chromosomes.

At the end of the evolution, an individual may have a genome that consists
of sequences inherited from multiple ancestors, which can be used to
count the number of ancestors who contribute genotype to this individual.
For example,

Ind 0: 111122222225555 000022222224444

has genotype from 5 ancestors (0, 1, 2, 4, 5). This script outputs
the average number of ancestors at each generation, for a given
number of replicate populations.

"""

import simuOpt, types
simuOpt.setOptions(alleleType='long')
from simuPOP import *

options = [
    {
     'name':'size',
     'default':1000,
     'label':'Population Size',
     'type':[int, long],
     'validator':simuOpt.valueGT(0),
     'description':'population size'
    },
     {
     'name':'lociNum',
     'default':100,
     'label':'Number of Loci',
     'type':[int, long],
     'validator':simuOpt.valueGT(0),
     'description':'Number of loci'
    },
    {
     'name':'gen',
     'default':20,
     'type':[int, long],
     'label':'Generations to evolve',
     'description':'Length of evolution',
     'validator':simuOpt.valueGT(0)
    },
    {
     'name':'recRate',
     'default':0.001,
     'label':'Recombination Rate',
     'type':[float],
     'description':'Recombination rate',
     'validator':simuOpt.valueBetween(0., 1.),
    },
    {
     'name':'rep',
     'default':5,
     'label':'Number of Replicate',
     'type':[int, long],
     'description':'Number of replicates',
     'validator':simuOpt.valueGT(0)
    }
]

def avgNumOfAncestors(pop):
    totAnc = 0
    for indi in pop.individuals():
        totAnc += len(set(indi.genotype()));
    print '%.2f\t' % (totAnc * 1.0 / pop.popSize()),
    return True

def simuNumOfAncestors(popSize, lociNum, gen, recRate, numRep):
    '''
    popSize: population size
    lociNum: number of loci
    gen:     generations to evolve
    recRate: recombination rate
    numRep:  number of replicates
    '''
    pop = Population(size=popSize, ploidy=2, loci=[lociNum])

    # initialize each individual with a different allele
    for idx,ind in enumerate(pop.individuals()):
        ind.setGenotype([idx])
    #
    simu = Simulator(pop, rep = numRep)
    simu.evolve(
        initOps = InitSex(),
        matingScheme = RandomMating(ops=Recombinator(rates=recRate)),
        postOps=[
            PyOperator(func=avgNumOfAncestors),
            PyOutput('\n', reps=-1)
        ],
        gen=gen
    )

if __name__ == '__main__':
    pars = simuOpt.Params(options, doc =
        'This script counts the average number of ancestors who contribute\n' +
        'their genotype to an offspring after a few generations.',
        details = __doc__)

    if not pars.getParam():
        sys.exit(0)

    simuNumOfAncestors(pars.size, pars.lociNum, pars.gen, pars.recRate, pars.rep)

