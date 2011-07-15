#
# Demonstrate the change of allele frequencies among a population due to
# mutation.
#

"""
This programs attempts to model the frequency of allele change due to mutation
equilibrium.
"""


import os, sys, types, time
import simuOpt
from simuPOP import *

# start simuPOP program

options = [
    {
     'name':'PopSize',
     'default':2000,
     'label':'Population Size',
     'type':[int, long],
     'validator':simuOpt.valueGT(0),
     },
    {
     'name':'m',
     'default':0.001,
     'label':'Mutation Rate',
     'type':[float],
     'validator':simuOpt.valueBetween(0., 1.),
     },
    {
     'name':'generations',
     'default':500,
     'label':'Generations to evolve',
     'description':'Length of evolution',
     'type':[int, long],
     'validator':simuOpt.valueGT(0)
     },
    {
     'name':'step',
     'default':100,
     'label':'Steps to take per generation',
     'description':'Values displayed per generation',
     'type':[int, long],
     'validator':simuOpt.valueGT(0)
     },
]

def simuMigration(PopSize, m, generations, step):

    pop = Population(size=[PopSize], loci=[1, 1], infoFields='fitness')
# initialize Population
#   set population size
#   two loci of interest are looked at
#   create field at where fitness values can be stored

    pop.evolve(
        initOps = InitSex(),
    # before evolve function takes place
    #   initiate population with males and females
        preOps = [
        # begin evolve function
            SNPMutator(u=m),
        # mutation function with rate for "A -> a" set earlier
            MaSelector(loci=0, fitness=[1, 1-10*m, 1-10*2*m]),
        ],
        matingScheme = RandomMating(),
        postOps = [
        # apply a purifying selection pressure to the first locus
            Stat(alleleFreq=[0, 1], step=step),
        # Stat function taken set amount of steps
            PyEval(
                r"'Generation: %.3f\n' % (gen)",
                step=step),
            PyEval(
                r"'frequency of allele at locus 1: %.3f\n' % (alleleFreq[0][1])",
                step=step),
            PyEval(
                r"'frequency of allele at locus 2: %.3f\n' % (alleleFreq[1][1])",
                step=step),
        # output freq of allele 1 at the first and second locus
        ],
        gen = generations
    # program run over set amount of generations
)

if __name__ == '__main__':
    # get all parameters
    pars = simuOpt.Params(options, __doc__)
    if not pars.getParam():
        sys.exit(0)

    simuMigration(pars.PopSize, pars.m, pars.generations, pars.step)
    # wait five seconds before exit

