#!/usr/bin/env python
#
# Purpose:    simulate the forming and detecton of recombination
#     hotspots.
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

Several samples will be drawn from the final population and be saved in a format
that can be analysis by LDhat (McVean, 2004). You can also modify the output
function and output to other recombination-estimation programs.

The program is written in Python using simuPOP modules. For more information,
please visit simuPOP website http://simupop.sourceforge.net .

Users can learn from this script how to write simuPOP populations into
another format.

"""

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
from simuUtil import *
import os, sys, types, exceptions, os.path, math

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
     'default': 500,
     'label': 'Initial population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial population size. This size will be maintained
                till the end of burnin stage''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'finalSize=',
     'default': 20000,
     'label': 'Final population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Final population size after population expansion.',
     'validate':    simuOpt.valueGT(0)
    }, 
    {'longarg': 'burnin=',
     'default': 200,
     'label': 'Length of burn-in stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of generations of the burn in stage.',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'noMigrGen=',
     'default': 150,
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
     'default': 1,
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
    {'longarg': 'sampleSize=',
     'default': 100,
     'label':    'Sample size',
     'allowedTypes':    [types.IntType, types.LongType],
     'description':    '''Size of the samples, that will mean N/4 affected sibpair families (of size 4),
                N/2 cases and controls etc. ''',
     'validate':    simuOpt.valueGT(1)
    },
    {'longarg': 'numSample=',
     'default': 2,
     'label':    'Sample number',
     'allowedTypes':    [types.IntType, types.LongType],
     'description':    '''Number of samples to draw for each penetrance function. ''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'sampleName=',
     'default': 'sample',
     'label':    'Sample name',
     'allowedTypes':    [types.StringType],
     'description':    '''Sample name. Sample index and suffix will be appended ''',
    },
    {'longarg': 'outputDir=',
     'default': '.',
     'allowedTypes': [types.StringType],
     'label': 'Output directory',
     'description': 'Directory into which the datasets will be saved. ',
     'validate':    simuOpt.valueValidDir()
    },
    {'longarg': 'configFile=',
     'default': sys.argv[0].split('.')[0]+'.cfg',
     'allowedTypes': [types.StringType, types.NoneType],
     'label': 'Save configuration',
     'description': 'Save current paremeter set to specified file.'
    },
]

def LinearExpansion(initSize, endSize, end, burnin=0, split=0, numSubPop=1, bottleneckGen=-1, bottleneckSize=0):
    '''
    Linearly expand population size from intiSize to endSize
    after burnin, split the population at generation split.
    '''
    inc = (endSize-initSize)/float(end-burnin)
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
            tot = initSize + inc*(gen-burnin)
        #
        if gen < split:
            return [int(tot)]
        elif gen > end:
            return [int(endSize/numSubPop)]*numSubPop
        else:
            return [int(tot/numSubPop)]*numSubPop
    return func

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


def getOptions(details = __doc__):
    ''' get options from options structure,
        if this module is imported, instead of ran directly,
        user can specify parameter in some other way.
    '''
    # get all parameters, __doc__ is used for help info
    allParam = simuOpt.getParam(options, 
        '''    This program simulates the evolution of a set SNP loci, subject 
     to the impact of mutation, migration, recombination and population size change. 
     Click 'help' for more information about the evolutionary scenario.''', details, nCol=2)
    # when user click cancel ...
    if len(allParam) == 0:
        sys.exit(1)
    # -h or --help
    if allParam[0]:    
        print simuOpt.usage(options, __doc__)
        sys.exit(0)
    # --saveConfig
    if allParam[-2] != None: # saveConfig
        simuOpt.saveConfig(options, allParam[-2], allParam)
    # --verbose or -v (these is no beautifying of [floats]
    if allParam[-1]:                 # verbose
        for p in range(len(options)):
            if options[p].has_key('label'):
                if type(allParam[p]) == types.StringType:
                    print options[p]['label'], ':\t"'+str(allParam[p])+'"'
                else:
                    print options[p]['label'], ':\t', str(allParam[p])
    # return the rest of the parameters
    return allParam[1:-2]


def SaveLDhat(pop, filename):
    ' save in .seq and .loc files that can be handled by LDhat directly '
    m = 100    # 100 characters each line
    seq = open(filename+'.seq', 'w')
    # haplotype phased data
    seq.write('%d %d %d\n' % (pop.popSize()*pop.ploidy(), pop.totNumLoci(), 1) )
    # output each haplotype
    for i in range(pop.popSize()):
        ind = pop.individual(i)
        for p in range(pop.ploidy()):
            seq.write('>genotype_%d_%d\n' % (i, p))
            gt = ind.genotype(p)
            for line in range(len(gt)/m):
                if m*(line +1) <= len(gt):
                    seq.write( ''.join([str(g) for g in gt[(m*line):(m*(line+1))]]) + '\n')
                else:
                    seq.write( ''.join([str(g) for g in gt[(m*line):]]) + '\n')
    seq.close()
    # write loc file
    loc = open(filename+'.loc', 'w')
    loc.write('%d %d %s\n' % (pop.totNumLoci(), pop.locusPos(pop.totNumLoci()-1), 'L') )
    for l in range(pop.totNumLoci()):
        loc.write('%f\n' % pop.locusPos(l)) 
    loc.close()
                
    
# simulate function, 
def simuHotSpot(numLoci, lociPos, initSize, finalSize, burnin, noMigrGen, mixingGen, 
        growth, numSubPop, migrModel, migrRate, mutaRate, recRate,
        sampleSize, numSample, sampleName, outputDir):
    ''' run the simulation, parameters are:
        numLoci:        number of SNP loci on the only chromosome
        lociPos:        loci position on the chromosome. Should be in 
            increasing order.
        initSize:     initial size
        finalSize:    ending population size
        burnin, noMigrGen, mixingGen: length of three stages
        growth, mumSubPop, migrModel: migration related parameters
        migrRate, muaRate, recRate:     rates
        sampleSize, numSample, sampleName: sample related
        outputDir:    output result to this directory. (*.bin, *.seq, *.loc)
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
    # simulator
    simu = simulator( pop, 
        randomMating(subPopSize=popSizeFunc ),
        rep = 1)
    # evolve! If --dryrun is set, only show info
    simu.evolve( 
        preOps = [
            # initialize all loci with two haplotypes (111,222)
            initByValue(value=[[x]*numLoci for x in range(1,3)],
                proportions=[.5]*2)
            ],
        ops = [
            # k-allele model for mutation of SNP
            kamMutator(rate=mutaRate, maxAllele=2),
            # recombination rate
            recombinator(rate=recRate),
            # split population after burnin, to each sized subpopulations
            splitSubPops(0, proportions=[1./numSubPop]*numSubPop, at=[split]),
            # migration
            migrOp,
            # report statistics
            stat(popSize=True, alleleFreq=[0,1], LD=[0,1], step=10),
            # report progress
            pyEval(r'"Generation %d, population size %d, allelefre=%.3g, %.3g, LD=%.3g\n" % (gen, popSize, alleleFreq[0][1], alleleFreq[1][1], LD_prime[0][1])', step=10),
            # show elapsed time
            ticToc(at=[ split, mixing, endGen])            
            ],
        gen = endGen
    )
    # get the population
    pop = simu.population(0)
    # save population 
    pop.save(os.path.join(outputDir, sampleName + ".pop"))
    # random sample
    sample = RandomSample(pop, size = sampleSize, times = numSample)
    for s in range(len(sample)):
        SaveLDhat(sample[s], os.path.join(outputDir, sampleName) + str(s))


if __name__ == '__main__':
        # get parameters
    par = simuOpt.simuOpt(options,
        '''This program simulates the evolution of a set SNP loci, subject 
         to the impact of mutation, migration, recombination and population size change. 
         Click 'help' for more information about the evolutionary scenario.''',
        __doc__)
   
    if not par.getParam():
        sys.exit(1)
    #
    par.saveConfig(par.configFile)
    simuHotSpot(*par.asList()[:-1])
    
