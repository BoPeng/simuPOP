#!/usr/bin/env python
#
# Demonstrate changes of allele frequency due to genetic drift. 

"""
This program demonstrates changes of allele frequency on single locus due to genetic drift.
"""

import simuOpt, os, sys, time
from simuPOP import *

try:
    from simuPOP.plotter import VarPlotter
except:
    print "simuRPy import failed. Please check your rpy installation."
    print "Allele Frequencies will not be plotted"
    useRPy = False
else:
    useRPy = True

options = [
    {'name':'popSize',
     'default':100,
     'label':'Population Size',
     'type': 'integer',
     'validator': 'popSize > 0',
     },
    {'name':'p',
     'default':0.2,
     'type': 'number',
     'label':'Initial Allele Frequency',
     'validator': 'p > 0 and p < 1',
     },
    {'name':'generations',
     'default':100,
     'label':'Number of Generations',
     'type': 'integer',
     'validator': 'generations > 0',
     },
    {'name':'replications',
     'default':5,
     'label':'Number of Replicates',
     'type': 'integer',
     'validator': 'replications > 0',
     },
]


def simuGeneticDrift(popSize=100, p=0.2, generations=100, replications=5):
    '''Simulate the Genetic Drift as a result of random mating.'''
    # diploid population, one chromosome with 1 locus
    # random mating with sex
    pop = Population(size=popSize, loci=[1])
    simu=Simulator(pop, rep=replications)

    if useRPy:
        plotter = VarPlotter('alleleFreq[0][0]', ylim=[0, 1], ylab='allele frequency',
            update=generations-1, saveAs='geneticDrift.png')
    else:
        plotter = NoneOp()

    # if number of generation is smaller than 200, step is 10 generations,
    # if it's between 200 and 500, set step to be 20 generations,
    # otherwise, step = 50 generations.
    if generations <= 200:
        s = 10
    elif 200 < generations <= 500:
        s = 20
    else:
        s = 50
        
    simu.evolve(
        # everyone initially will have the same allele frequency
        initOps = [
            InitSex(),
            InitGenotype(freq=[p, 1-p])
        ],
        matingScheme = RandomMating(),
        postOps = [
            Stat(alleleFreq=[0]),
            PyEval(r'"Generation %d:\t" % gen', reps = 0, step = s),
	        PyEval(r"'%.3f\t' % alleleFreq[0][0]", step = s),
	        PyOutput('\n', reps=-1, step = s),
	        plotter,
        ],
        gen = generations
    )

if __name__ == '__main__':
    # get all parameters
    pars = simuOpt.Params(options, __doc__)
    # cancelled
    if not pars.getParam():
        sys.exit(0)
      
    simuGeneticDrift(pars.popSize, pars.p, pars.generations, pars.replications)

    # wait ten seconds before exit
    if useRPy:
        print "Figure will be closed after 5 seconds."
        time.sleep(5)
        
    

    
