#!/usr/bin/env python
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedDate: 2005-10-20 15:26:39 -0500 (Thu, 20 Oct 2005) $
# $Rev: 72 $
#
# Known bugs:
#     None
# 

"""

Introduction
=============

This program simulates the evolution of a chromosome with SNP markers, under 
the influence of mutation, migration, recombination and population size 
change. Starting from a small founder population, each simulation will go
through the following three stages:

    1. Burn-in the population with mutation and recombination
    2. Split and grow the population without migration
    3. Mix subpopulations at given migration level

Steps 2 and 3 are optional in the sense that you can manipulate the parameters
to simulate constant populations without population structure.

The program is written in Python using simuPOP modules. For more information,
please visit simuPOP website http://simupop.sourceforge.net .

"""

import simuOpt
simuOpt.setOptions(quiet=True, alleleType='binary')

from simuPOP import *
from simuPOP.utils import *
import os, sys, types, exceptions, os.path , math

try:
    from simuRPy import *
    hasRPy = True
except:
    print "RPy module not found. Can not plot histogram of allele frequency"
    hasRPy = False
#
# declare all options. getParam will use these information to get parameters
# from a tk/wxPython-based dialog, command line, config file or user input
#
# details about these fields is given in the simuPOP reference manual.
options = [
    {'longarg': 'numLoci=',
     'default': 200,
     'label': 'Number of SNP loci',
     'description': 'Number of SNP loci',
     'allowedTypes': [types.IntType, types.LongType],
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'lociPos=',
     'default': [],
     'label': 'Loci position',
     'description': '''Loci position on the chromosome. Assumed to be in kb
        the unit is not important.''',
     'allowedTypes': [types.ListType, types.TupleType],
     'validate':    simuOpt.valueListOf(types.FloatType)
    },
    {'longarg': 'initSize=',
     'default': 1000,
     'label': 'Initial population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial population size. This size will be maintained
        till the end of burnin stage''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'finalSize=',
     'default': 100000,
     'label': 'Final population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Final population size after population expansion.',
     'validate':    simuOpt.valueGT(0)
    }, 
    {'longarg': 'burnin=',
     'default': 1000,
     'label': 'Length of burn-in stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of generations of the burn in stage.',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'noMigrGen=',
     'default': 1500,
     'label': 'Length of split-and-grow stage',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Number of generations when migration is zero. This stage
                is used to build up population structure.''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'mixingGen=',
     'default': 50,
     'label': 'Length of mixing stage',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Number of generations when migration is present. This stage
                will mix individuals from subpopulations using an circular stepping stone
                migration model.''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'growth=',
     'default': 'exponential',
     'label': 'Population growth model',
     'description': '''How population is grown from initSize to finalSize.
                Choose between linear and exponential''',
     'chooseOneOf': ['exponential', 'linear'],
    },
    {'longarg': 'numSubPop=',
     'default': 5,
     'label': 'Number of subpopulations to split',
     'allowedTypes': [types.IntType],
     'description': 'Number of subpopulations to be split into after burnin stage.',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'migrModel=',
     'default': 'stepping stone',
     'label': 'Migration model',
     'allowedTypes': [types.StringType],
     'description': '''Migration model. Choose between stepping stone (circular),
                island and none. ''',
     'validate':    simuOpt.valueOneOf(['island', 'stepping stone', 'none']),
     'chooseOneOf': ['stepping stone', 'island', 'none']
    }, 
    {'longarg': 'migrRate=',
     'default': 0.,
     'label': 'Migration rate',
     'description': '''Migration rate during mixing stage. ''',
     'allowedTypes': [types.IntType, types.FloatType],
     'validate':    simuOpt.valueBetween(0,1)
    },
    {'longarg': 'mutaRate=',
     'default': 1e-4,
     'label': 'Mutation rate',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': '''Mutation rate''',
     'validate': simuOpt.valueBetween(0,1)
    },
    {'longarg': 'recRate=',
     'default': [1e-4],
     'label': 'Recombination rate',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Recombination rate between adjacent markers. 
        It can be a number or a list with length numLoci-1''',
     'validate':    simuOpt.valueListOf( simuOpt.valueBetween(0,1))
    },
    {'longarg': 'name=',
     'default': 'simu',
     'allowedTypes': [types.StringType],
     'label': 'Name of the simulation',
     'description': '''Name of the simulation, configuration and output will be saved to a directory
        with the same name''',
    },
]


def ExponentialExpansion(initSize, endSize, end, burnin=0, split=0, numSubPop=1, bottleneckGen=-1, bottleneckSize=0):
    '''
    Exponentially expand population size from intiSize to endSize
    after burnin, split the population at generation split.
    '''
    rate = (math.log(endSize)-math.log(initSize))/(end-burnin)
    def func(gen, oldSize=[]):
        if gen == bottleneckGen:
            if gen < split:
                return [bottleneckSize]
            else:
                return [bottleneckSize/numSubPop]*numSubPop
        # not bottleneck
        if gen <= burnin:
            tot = initSize
        else:
            tot = int(initSize*math.exp((gen-burnin)*rate))
        if gen < split:
            return [int(tot)]
        elif gen > end:
            return [int(endSize/numSubPop)]*numSubPop
        else:
            return [int(tot/numSubPop)]*numSubPop
    return func

def InstantExpansion(initSize, endSize, end, burnin=0, split=0, numSubPop=1, bottleneckGen=-1, bottleneckSize=0):
    '''
    Instaneously expand population size from intiSize to endSize
    after burnin, split the population at generation split.
    '''
    def func(gen, oldSize=[]):
        if gen == bottleneckGen:
            if gen < split:
                return [bottleneckSize]
            else:
                return [bottleneckSize/numSubPop]*numSubPop
        # not bottleneck
        if gen <= burnin:
            tot = initSize
        else:
            tot = endSize
        if gen < split:
            return [int(tot)]
        else:
            return [int(tot/numSubPop)]*numSubPop
    return func

def plotAlleleFreq(pop, param):
    'plot the histogram of allele frequency'
    if not hasRPy:
        return True
    # 
    freq = pop.dvars().alleleFreq
    freq0 = [min(freq[i][0], 1-freq[i][0]) for i in range(pop.totNumLoci())]
    from rpy import r
    r.postscript(os.path.join(param[0], '%s_%d.eps' % (param[0], pop.dvars().gen)))
    r.hist(freq0, nclass=50, xlim=[0,0.5], xlab='frequency', ylab='hist', 
        main='Histogram of allele frequencies at generation %d' % pop.dvars().gen)
    r.dev_off()
    return True


# simulate function, 
def simulate( numLoci, lociPos, initSize, finalSize, burnin, noMigrGen, mixingGen, 
        growth, numSubPop, migrModel, migrRate, mutaRate, recRate, name):
    ''' run the simulation, parameters are:
        numLoci:        number of SNP loci on the only chromosome
        lociPos:        loci position on the chromosome. Should be in 
            increasing order.
        initSize:     initial size
        finalSize:    ending population size
        burnin, noMigrGen, mixingGen: length of three stages
        growth, mumSubPop, migrModel: migration related parameters
        migrRate, muaRate, recRate:     rates
        name:   name of the simulation, and output directory and file
    '''
    # event generations
    split    = burnin 
    mixing = split + noMigrGen
    endGen = split + noMigrGen + mixingGen    
    # demographic model
    if growth == 'linear':
        popSizeFunc = LinearExpansion(initSize, finalSize, endGen,
            burnin, split, numSubPop)
    elif growth == 'exponential':
        popSizeFunc = ExponentialExpansion(initSize, finalSize, endGen,
            burnin, split, numSubPop)
    else:
        raise exceptions.ValueError("Growth model can be one of linear and exponential.")
    # a migrator, stepping stone or island
    if numSubPop > 1 and migrModel == 'island' and migrRate > 0:
        migrOp = migrator(migrIslandRates(migrRate, numSubPop),
            mode=MigrByProbability, begin=mixing) 
    elif numSubPop > 1 and migrModel == 'stepping stone' and migrRate > 0:
        migrOp = migrator(migrSteppingStoneRates(migrRate, numSubPop, circular=True),
            mode=MigrByProbability, begin=mixing) 
    else:
        migrOp = noneOp()
    # population
    # default loci position is 1,2,3,..., with no unit
    if lociPos == []:
        lociPos = range(1, numLoci+1)
    if len(lociPos) != numLoci:
        print "If loci position is given, it should have length numLoci"
        sys.exit(1)
    pop = population(size=popSizeFunc(0), ploidy=2,
        loci = [numLoci], lociPos = lociPos)
    # save name in var
    pop.dvars().name = name
    # simulator
    simu = simulator( pop, 
        randomMating(subPopSize=popSizeFunc, ops = recombinator(rates=recRate)),
        rep = 1)
    simu.evolve( 
        initOps = [
            initSex(),
            # initialize all loci with two haplotypes (111,222)
            initByValue(value=[[x]*numLoci for x in range(2)],
                proportions=[.5]*2)
            ],
        preOps = [
            # k-allele model for mutation of SNP
            kamMutator(rates=mutaRate, k=2),
            # split population after burnin, to each sized subpopulations
            splitSubPops(0, proportions=[1./numSubPop]*numSubPop, at=[split]),
            # migration
            migrOp
        ],
        postOps = [
            # report statistics
            stat(popSize=True, alleleFreq=range(pop.totNumLoci()), step=10),
            # report progress
            pyEval(r'"Generation %d, population size %d, allelefre=%.3g, %.3g\n" % (gen, popSize, alleleFreq[0][1], alleleFreq[1][1])', step=10),
            # plot histogram
            pyOperator(func=plotAlleleFreq, param=(name,), step=50),
            # save population every 1000 generations
            savePopulation(output="!'%s/%s_%d.txt'%(name,name,gen)",step=1000),
            # show elapsed time
            ticToc(at=[ split, mixing, endGen])
            ],
        gen=endGen
    )
    # save population 
    pop.save(os.path.join(name, '%s.pop' % name))


if __name__ == '__main__':
    par = simuOpt.simuParam(options, 
        '''This program simulates the evolution of a set SNP loci, subject 
     to the impact of mutation, migration, recombination and population size change. 
     Click 'help' for more information about the evolutionary scenario.''', __doc__)
    if not par.getParam():
        sys.exit(1)
    #
    if not os.path.isdir(par.name):
        os.makedir(par.name)
    #
    par.saveConfig(os.path.join(par.name, par.name + '.cfg'))
    simulate(*par.asList())
    
