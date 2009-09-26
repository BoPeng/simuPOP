#
# Demonstrate the change of allele frequencies among a population due to
# mutation.
#

"""
This programs attempts to model the frequency of allele change due to mutation
equilibrium.
"""


import os, sys, types, time
from simuPOP import *

# start simuPOP program

options = [
    {
     'longarg':'PopSize=',
     'default':2000,
     'label':'Population Size',
     'allowedTypes':[types.IntType, types.LongType],
     'validate':params.valueGT(0),
     },
    {
     'longarg':'m=',
     'default':0.001,
     'label':'Mutation Rate',
     'allowedTypes':[types.FloatType],
     'validate':params.valueBetween(0., 1.),
     },
    {
     'longarg':'generations=',
     'default':500,
     'label':'Generations to evolve',
     'description':'Length of evolution',
     'allowedTypes':[types.IntType, types.LongType],
     'validate':params.valueGT(0)
     },
    {
     'longarg':'step=',
     'default':100,
     'label':'Steps to take per generation',
     'description':'Values displayed per generation',
     'allowedTypes':[types.IntType, types.LongType],
     'validate':params.valueGT(0)
     },
]

def simuMigration(PopSize, m, generations, step):

    pop = population(size=[PopSize], loci=[1, 1], infoFields='fitness')
# initialize population
#   set population size
#   two loci of interest are looked at
#   create field at where fitness values can be stored

    simu = simulator(pop, randomMating())
# simulate random mating within the population

    simu.evolve(
        preOps = initSex(),
    # before evolve function takes place
    #   initiate population with males and females
        ops = [
        # begin evolve function
            snpMutator(u=m),
        # mutation function with rate for "A -> a" set earlier
            maSelector(loci=0, fitness=[1, 1-10*m, 1-10*2*m]),
        # apply a purifying selection pressure to the first locus
            stat(alleleFreq=[0, 1], step=step),
        # stat function taken set amount of steps
            pyEval(
                r"'Generation: %.3f\n' % (gen)",
                step=step),
            pyEval(
                r"'frequency of allele at locus 1: %.3f\n' % (alleleFreq[0][1])",
                step=step),
            pyEval(
                r"'frequency of allele at locus 2: %.3f\n' % (alleleFreq[1][1])",
                step=step),
        # output freq of allele 1 at the first and second locus
        ],
        gen = generations
    # program run over set amount of generations
)

if __name__ == '__main__':
    # get all parameters
    pars = params.simuParam(options, __doc__)
    if not pars.getParam():
        sys.exit(0)

    simuMigration(pars.PopSize, pars.m, pars.generations, pars.step)
    # wait five seconds before exit
