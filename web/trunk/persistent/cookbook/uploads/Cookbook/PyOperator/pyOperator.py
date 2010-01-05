#!/usr/bin/env python

"""

This is a small example of how to use PyOperator to do
custimized mating and condition checking (terminator) 
etc. This is *not* the efficient way to do that since
everything can be handled using a post-mating operator.

Evolution scenario:

A disease spread like this: when one of the chromosomes
of a diploid individual has this allele, it will spread
to the other copy.

    x x x - - - - -     ==>    x x x x x - - -
    x x x x x - - -     ==>    x x x x x - - -

Purpose of simulation:

    Given N, random mating, how many generations to 
    reach all X status?

Note that many simulations will fail because all x
chromosomes disappeared because of genetic drift.
    
"""

import simuOpt
from simuPOP import *
import types

options = [
    {'arg': 'N:',
     'longarg': 'N=',
     'default': 10000,
     'label': 'population size',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': simuOpt.valueGT(0)
    },
    {'arg': 'r:',
     'longarg': 'rep=',
     'default': 10,
     'label': 'Replicates',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': simuOpt.valueGT(0)
    }
]


# this will be used by a duringMating pyOpertor,
# however, for this particular problem, you can 
# use a postMating operator and handle all invidiauls
# one by one. The performance will be *much* better.
#
def offGen(pop, off, dad, mom):
    ''' how to pass allele? 
        We can of course do everything by ourself, but if we do not
        set isTransmitter of this PyOperator, we can let RandomMating()
        generate offspring genotype as usual, we just need to change
        1 2 to 2 2.
    '''
    geno = off.genotype()
    if geno[0] + geno[1] == 3: # in the case of 1 2 or 2 1
        geno[0] = 2
        geno[1] = 2
    return True

def allTwos(pop):
    ''' see if the genotype of the population is all 2'''
    geno = sum(pop.genotype())
    # terminate if all 1 (set flag to fail
    if geno == pop.popSize()*2:     # all 1
        pop.dvars().succ = False
		# return false to terminate simulation
        return False
    # terminate if all 2 (set flag to success)
    elif geno == pop.popSize()*4: # all 2
        pop.dvars().succ = True
		# return false to terminate simulation
        return False
    else:
		# return true to continue simulation
        return True
         
def simu(N):
    ' run the simulation! '
    pop = Population(N, loci=[1], ploidy=2)
    initSex(pop)
    initGenotype(pop, genotype=[1])
    # you can also use a PointerMutator ...
    pop.individual(0).setAllele(2,0)
    pop.individual(0).setAllele(2,1)
    pop.evolve(
        matingScheme = RandomMating(ops=[
            MendelianGenoTransmitter(),
            PyOperator(func=offGen)]),
        postOps = PyOperator(func=allTwos),
    )
    return (pop.dvars().succ, pop.dvars().gen)

if __name__ == '__main__':
    pars = simuOpt.Params(options,
        '''This program mimic the evolution of the infection process
    where X- chromosome will be turn to XX. We are concerned    about 
    the speed at which all Population becomes XX. ''', __doc__)
    
    if not pars.getParam():
        sys.exit(1)
    
    # get the parameters
    N = pars.N
    rep = pars.rep

    succCount = 0
    genRecord = []
    for r in range(rep):
        print "Replicate ", r,
        (succ, gen) = simu(N) 
        if succ:
            succCount += 1
            genRecord.append(gen)
            print " success at generation ", gen
        else:
            print " failed"

    print "population size: ", N
    print "Replicates: ", rep
    print "Successful counts: ", succCount
    if succCount > 0:
        print "Mean generation: ", sum(genRecord)*1.0/succCount
        print "Minimal generation: ", min(genRecord)
        print "Maxmimal generation: ", max(genRecord)
    
