#!/usr/bin/env python
#
# Purpose:
#     This python module provides function simuGWAS to evolve a population
#     forward in time, subject to rapid population expansion, mutation,
#     recombination and natural selection.
#
# License:
#     This program is freely available in the simuPOP's online cookbook
#     (http://simupop.sourceforge.net/cookbook). You can redistribute it and/or
#     modify it freely, with or without this license notice. However, this
#     license notice should be present in the online version of this file.
#
#     This program is NOT part of simuPOP and is NOT protected by simuPOP's GPL
#     license. It is provided in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY. If you notice any bug or have some new ideas, you can
#     modify this file and, as a courtesy to other simuPOP users, incoporate
#     your changes to the online version of this file. If you are uncertain
#     about your changes, please feel free to discuss your changes in the
#     simuPOP mailinglist (simupop-list@lists.sourceforge.net, subscription
#     required).
#
# Change Log:
#     2010-03-01 Bo Peng <bpeng@mdanderson.org>
#         Add script.

'''
This script evolves a population forward in time, subject to rapid population
expansion, mutation, recombination and natural selection. A trajectory simulation
method is used to control the allele frequency of optional disease predisposing
loci. A scaling approach can be used to improve efficiency when weak, additive
genetic factors are used.
'''

import simuOpt
simuOpt.setOptions(alleleType='binary', version='1.0.1', optimized=False, quiet=True)
from simuPOP import *

from simuPOP.utils import simulateForwardTrajectory, simulateBackwardTrajectory, migrSteppingStoneRates

import os, sys, math, exceptions, types, time

options = [
    {'longarg': 'initPop=',
     'default': 'init.pop',
     'label': 'Initia population',
     'allowedTypes': types.StringType,
     'validate': simuOpt.valueValidFile(),
     'description': 'An initial population created using script selectedMarkers.py'
    },
    {'longarg': 'dumpRec=',
     'default': '',
     'allowedTypes': types.StringType,
     'description': '''This is a hidden option, if set, all recombination events will
            be dumped to the specified file. Please refer to the simuPOP user's guide
            for explanation of this file.'''
    },
    {'longarg': 'trajFile=',
     'default': '',
     'allowedTypes': types.StringType,
     'description': '''This is a hidden option which, if set from command line, will
            output the simulated trajectory to the specified file. The format
            is generation freq.'''
    },
    {'longarg': 'haploCount=',
     'default': [0,0],
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''This is an hidden option, if set, the number of unique haplotypes
            between specified loci will be outputted.'''
    },
    {'longarg': 'filename=',
     'label': 'File to save expanded population',
     'default': 'expanded.pop',
     'useDefault': True,
     'description': '''File to save the expanded population.''',
     'allowedTypes': [types.StringType],
    },
    {'separator': 'Mutation, recombination, etc'},
    {
     'longarg': 'mutaRate=',
     'default': 1e-8,
     'useDefault': True,
     'label': 'Mutation rate',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': '''Mutation rate at all markers. This value will be scaled
            by the scaling parameter.''',
     'validate': simuOpt.valueBetween(0,1),
    },
    {'longarg': 'recIntensity=',
     'default': 1e-8,
     'useDefault': True,
     'allowedTypes': [types.FloatType, types.IntType, types.LongType],
     'label': 'Recombination intensity',
     'description': '''Recombination intensity. The actual recombination rate
            between two adjacent loci will be this intensity times loci distance
            in basepair, and times the scaling parameter. For example, if two
            loci are 1M bp apart, the recombination rate will be 1e6 x 1e-8 (default
            value), namely 1e-2 times scaling parameter. (This is roughly 0.01 
            per cM). If a variable 'geneticMap' that stores the genetic map of each
            locus exists in this population, the genetic map will be used.''', 
     'validate': simuOpt.valueGE(0),
    },
    {'longarg': 'scale=',
     'default': 1,
     'useDefault': True,
     'allowedTypes': [types.IntType, types.LongType, types.FloatType],
     'label': 'Acceleration scale',
     'description': '''This parameter is used to speed up recombination, mutation
                and selection. Briefly speaking, certain parts of the evolutionary
                process is accelerated but random genetic drift is kept. Please
                refer to Peng 2008 for more details''',
     'validate': simuOpt.valueGT(0),
    },
    {'separator': 'Demographic model'},
    {'longarg': 'expandGen=',
     'default': 500,
     'useDefault': True,
     'label': 'Generations to expand',
     'description': '''Number of generations to evolve during the population
                expansion stage. The actual evolved population is scaled down by
                parameter --scale. (If scale==10, expandGen=1000, the actually
                evolved generation is 100).''',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'initSize=',
     'default': 0,
     'useDefault': False,
     'description': '''Size of first generation of the evolution. If this is zero, the
                size of the initial population will be used. Otherwise, the first
                generation will have specified number of individuals. This is to avoid
                losing genetic information if the initial population is very small.''',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': simuOpt.valueGE(100)
    },
    {'longarg': 'expandSize=',
     'default': 50000,
     'useDefault': True,
     'label': 'Expanded population size',
     'description': '''Size of the expanded population. The default value is the recommended
                value when all hapmap populations are used (60+60+90)*200 and the simulation
                runs with a scaling factor of 1 (no scaling). When scaling is used, the
                expanded population size should also be reduced (divided by the scaling factor)
                to simulate a population with comparable properties of the unscaled one.
                Otherwise, the simulation will appear to go through faster population expansion
                which leads to reduced genetic drift.''',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': simuOpt.valueGE(100)
    },
    {'longarg': 'migrRate=',
     'default': 0.001,
     'useDefault': True,
     'label': 'Migration rate ',
     'description': '''If there are more than one subpopulations, this parameter
                controlls the rate at which migration is allowed.
                A circular stepping stone migration model will be used. ''',
     'allowedTypes': [types.IntType, types.FloatType],
     'validate':    simuOpt.valueBetween(0,1)
    },
    {'separator': 'Disease information'},
    {'longarg': 'DPL=',
     'default': [],
     'useDefault': True,
     'label': 'Names of disease predisposing loci',
     'description': '''A list of names of disease predisposing loci. These loci
                must exist in the population.''',
     'allowedTypes': [types.TupleType, types.ListType],
     'validate': simuOpt.valueListOf(types.StringType)
    },
    {'longarg': 'curAlleleFreq=',
     'default': [0.05],
     'useDefault': True,
     'label': 'Final allele frequencies',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Current allele frequencies for each DPL.
                If a number is given, it is assumed to be the frequency
                for all DPL.''',
     'validate': simuOpt.valueListOf( simuOpt.valueBetween(0,1))
    },
    {'longarg': 'trajectory=',
     'default': 'forward',
     'useDefault': True,
     'label': 'Trajectory simulation method',
     'allowedTypes': types.StringType,
     'chooseOneOf': ['forward', 'backward'],
     'description': '''Trajectory simulation type. A forward method assumes that the mutants
                are old so the trajectory will start from the existing allele frequency.
                A backward method assumes that the mutants are relatively new (less than
                expandGen) so the exiting alleles are cleared and will be introduced later.''',
    },
    {'longarg': 'fitness=',
     'default': [1, 1, 1],
     'useDefault': True,
     'label': 'Fitness of genotype AA,Aa,aa',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Fitness of genotype, can be:
                f1, f2, f3: if one DPL, the fitness for genotype AA, Aa and aa
                f1, f2, f3: if multiple DPL, the same fitness for each locus
                [a1, a2, a3, b1, b2, b3, ...] if mlSelModel = 'additive' 
                    or multiplicative, fitness at each locus. The overall fitness
                    is determined by mlSelModel
                [a1, a2, a3, b1, b2, b3, c1, c2, c3, ...] and mlSelModel = interaction.
                    For example, in the 2-DPL case, the numbers are (by row)
                        BB Bb bb
                    AA  a1 a2 a3
                    Aa  b1 b2 b3
                    aa  c1 c2 c3
                3^n numbers are needed for n DPL.
        ''',
     'validate':    simuOpt.valueListOf(simuOpt.valueGE(0.)),
    },
    {'longarg': 'mlSelModel=',
     'default': 'none',
     'label': 'Multi-locus selection model',
     'useDefault': True,
     'description': '''Model of overall fitness value given fitness values for each DPL.
                multiplicative: f =  Prod(f_i)
                additive: f = 1-Sum(1-f_i)
                interaction: the intepretation of fitness parameter is different.
                    see fitness.
                ''',
     'allowedTypes': [types.StringType],
     'chooseOneOf': [ 'additive', 'multiplicative', 'interaction', 'none']
    },
]


def expDemoFunc(N0, N1, expandGen):
    '''
    Return an exponential population expansion demographic function.

    N0: a list of initial subpopulation sizes.
    N1: final population size.
    
    gen: generations to evolve.
    '''
    if type(N1) in [types.IntType, types.LongType]:
        N1 = [int(N1*1.0*x/sum(N0)) for x in N0]
    elif len(N1) != len(N0):
        raise exceptions.ValueError("Number of subpopulations should be the same")
    #
    step = [float(x-y) / expandGen for x,y in zip(N1, N0)]
    def func(gen):
        # try to make the exact population size
        if gen == expandGen - 1:
            return N1
        else:
            return [int(x + (gen + 1) * y) for x, y in zip(N0, step)]
    return func

def Ne(demoFunc, pars):
    ne = 0
    for i in range(pars.expandGen):
        ne += 1./sum(demoFunc(i))
    return int(pars.expandGen / ne)

def simuGWAS(pars, pop, logger=None):
    '''The main program'''
    # record the starting time:if expand("%") == ""|browse confirm w|else|confirm w|endif
    #
    startTime = time.time()
    #
    pars.mutaRate *= pars.scale
    pars.recIntensity *= pars.scale
    pars.migrRate *= pars.scale
    pars.expandGen = int(pars.expandGen / pars.scale)
    pars.fitness = [1 + (x-1) * pars.scale for x in pars.fitness]
    pop.dvars().scale = pars.scale
    #
    if pars.initSize == 0:
        demoFunc = expDemoFunc(pop.subPopSizes(), pars.expandSize, pars.expandGen)
    else:
        if type(pars.initSize) in [type(()), type([])] and len(pars.initSize) != pop.numSubPop():
            raise ValueError('If initial population size is specified, it should be specified for all subpopulations.')
        if pop.numSubPop() == 1:
            if type(pars.initSize) not in  [type(()), type([])]:
                pars.initSize = [pars.initSize]
        else:
            if type(pars.initSize) in [type(0), type(0L)]:
                pars.initSize = [pars.initSize*1.0*x/pop.popSize() for x in pop.subPopSizes()]
        demoFunc = expDemoFunc(pars.initSize, pars.expandSize, pars.expandGen)
    if logger:
        logger.info('Expected effective population size at generation %d is %d' %
            (pars.expandGen, Ne(demoFunc, pars)))

    #
    # define a trajectory function
    pars.DPL_idx = pop.lociByNames(pars.DPL)
    if len(pars.DPL) > 0:
        numDPL = len(pars.DPL)
        # curAlleleFreq expand 0.5 -> [0.5,0.5,...]
        if len(pars.curAlleleFreq) == 1:
            pars.curAlleleFreq = pars.curAlleleFreq * numDPL
        # fitness
        if pars.mlSelModel == 'none':
            if numDPL > 1:
                raise exceptions.ValueError('A multi-locus selection model is needed if multiple disease predisposing loci are specified.')
            if len(pars.fitness) != 3:
                raise exceptions.ValueError('Fitness value should have three elemenst. %s given.' % (pars.fitness,))
            fitness = pars.fitness
        elif pars.mlSelModel == 'interaction':
            if numDPL == 1:
                raise exceptions.ValueError("Interaction model can only be used with more than one DPL");
            if len(pars.fitness) != 3**numDPL:
                raise exceptions.ValueError("Please specify 3^n fitness values for n DPL");
            fitness = pars.fitness
        else:
            if pars.fitness == []:    # neutral process
                fitness = [1,1,1]*numDPL
            else:
                # for a single DPL
                if len(pars.fitness) == 3:
                    fitness = pars.fitness*numDPL
                elif len(pars.fitness) == numDPL*3:
                    fitness = pars.fitness
                else:
                    raise exceptions.ValueError("Please specify fitness for each DPL")
        #
        if pars.trajectory == 'forward':
            # has to be all backward or all forward
            stat(pop, alleleFreq=pars.DPL_idx, vars='alleleFreq_sp')
            currentFreq = []
            # in the order: LOC0: sp0, sp1, sp2, LOC1: sp0, sp1, sp2, ...
            for sp in range(pop.numSubPop()):
                for idx,loc in enumerate(pars.DPL_idx):
                    currentFreq.append(pop.dvars(sp).alleleFreq[loc][1])
            #
            endFreq=[(x - min(0.01, x/5.), x + min(0.01, x/5., (1-x)/5.)) for x in pars.curAlleleFreq]
            if logger:
                logger.info('Simulating trajectory forward in time, from allele frequency %s to frequency (range) %s' % \
                    (currentFreq, endFreq))
            traj = simulateForwardTrajectory(
                N=demoFunc,
                beginGen=0,
                endGen=pars.expandGen,
                beginFreq=currentFreq,
                endFreq=endFreq,
                nLoci=len(pars.DPL),
                fitness=fitness,
                maxAttempts=1000,
                logger=logger
            )
            introOps = [] 
        else:
            # clear existing mutants at these loci
            pop.recodeAlleles(alleles=[0,0], loci=pars.DPL_idx)
            if logger:
                logger.info('Simulating trajectory backward in time, with ending allele frequency %s' % pars.curAlleleFreq)
            traj = simulateBackwardTrajectory(
                N=demoFunc,
                endGen=pars.expandGen,
                endFreq=pars.curAlleleFreq,
                nLoci=len(pars.DPL),
                fitness=fitness,
                minMutAge=1,
                maxMutAge=pars.expandGen,
                logger=logger
            )
            introOps = traj.mutators(loci=pars.DPL_idx)
        if traj is None:
            raise SystemError('Failed to generated trajectory after 10000 attempts. '
                'This usually means that the demographic and genetic settings are '
                'very extreme which makes it very likely for an allele to reach designed'
                'allele frequency. Please adjust your parameters and try again.')
        trajFunc = traj.func()
        if pars.trajFile:
            out = open(pars.trajFile, 'w')
            for gen in range(pars.expandGen):
                print >> out, gen, ' '.join(['%.4f' % x for x in trajFunc(gen)])
            out.close()
    else:
        trajFunc = None
        introOps = []
    # recombination
    try:
        # try to use a genetic map
        pos = [pop.dvars().geneticMap[pop.locusName(x)] for x in range(pop.totNumLoci())]
        loc = []
        rate = []
        for ch in range(pop.numChrom()):
            beg = pop.chromBegin(ch)
            end = pop.chromEnd(ch)
            loc.extend(range(beg, end - 1))
            # genetic map use cM/Mb as combined rate so we need to 
	    # 1e-6 to reset recIntensity.
            rate.extend([(pos[x+1] - pos[x])*1e6*pars.recIntensity for x in range(beg, end - 1)])
        if pars.dumpRec:
            recOp = Recombinator(rates=rate, loci=loc, infoFields='ind_id', output='>>%s' % pars.dumpRec)
        else:
            recOp = Recombinator(rates=rate, loci=loc)
        print 'Scaled recombination at %.3f cM/Mb over %.2f centiMorgan genetic (%.0f bp physical) distance (first chromosome)' % \
            (pars.recIntensity*1e6, (pop.dvars().geneticMap[pop.locusName(pop.numLoci(0)-1)] - \
                pop.dvars().geneticMap[pop.locusName(0)]),
                int(pop.locusPos(pop.numLoci(0)-1) - pop.locusPos(0)))
    except:
        # if not possible, use a physical map
        if pars.dumpRec:
            recOp = Recombinator(intensity=pars.recIntensity, infoFields='ind_id', output='>>%s' % pars.dumpRec)
        else:
            recOp = Recombinator(intensity=pars.recIntensity)
    #
    if pars.dumpRec != '':
        pop.addInfoFields('ind_id')
    if pop.numSubPop() > 1:
        pop.addInfoFields('migrate_to')
        migrOp = Migrator(rate=migrSteppingStoneRates(pars.migrRate, pop.numSubPop()))
    else:
        migrOp = NoneOp()
    # recombination rate at the end of each chromosome will be invalid
    # but this does not matter
    #
    pop.dvars().DPL = pars.DPL_idx
    pop.dvars().scale = pars.scale
    pop.evolve(
        initOps = [
            InitSex(),
            IfElse(pars.dumpRec != '', IdTagger()),
        ],
        preOps = [
            SNPMutator(u=pars.mutaRate, v=pars.mutaRate, loci=range(pop.totNumLoci())),
            migrOp,
        ] + introOps,
        matingScheme = ControlledRandomMating(
            loci = pars.DPL_idx,
            alleles = [1]*len(pars.DPL_idx),
            freqFunc = trajFunc,
            ops = [
                IfElse(pars.dumpRec != '', IdTagger()),
                recOp,
            ],
            subPopSize = demoFunc),
        postOps = [
            Stat(popSize = True, alleleFreq=pars.DPL_idx, structure=range(pop.totNumLoci())),
            IfElse(len(pars.DPL) != 0,
                PyEval(r'"After %3d generations, size=%s, freq=%s\n" % ((gen + 1 ) * scale, subPopSize, ", ".join(["%.4f" % alleleFreq[x][1] for x in DPL]))'),
                PyEval(r'"After %3d generations, size=%s\n" % ((gen + 1 )* scale, subPopSize)')
            ),
            IfElse(pop.numSubPop() > 1,
                PyEval(r"'F_st = %.3f\n' % F_st", step=10)),
            IfElse(pars.haploCount[1] > pars.haploCount[0],
                ifOps=[
                    Stat(haploFreq=range(pars.haploCount[0], pars.haploCount[1])),
                    PyEval(r'"Number of haplotypes between loci %d and %d is %%d\n" %% len(haploNum[%s])' % \
                        (pars.haploCount[0], pars.haploCount[1], tuple(range(pars.haploCount[0], pars.haploCount[1]))))
                ],
                step=10,
            ),
        ],
        gen = pars.expandGen
    )
    if logger:
        logger.info('Simulation finishes in %d seconds.' % (time.time() - startTime))
    return pop


short_desc = '''This program evolves a subset of the HapMap dataset
forward in time in order to simulate case control samples with realistic
allele frequency and linkage disequilibrium structure.'''

if __name__ == '__main__':
    # get all parameters
    pars = simuOpt.Params(options, short_desc, __doc__)
    # when user click cancel ...
    if not pars.getParam():
       sys.exit(1)
    # create a logger
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()
    logger.info('Loading initial population from file %s. ' % pars.initPop)
    init = loadPopulation(pars.initPop)
    pop = simuGWAS(pars, init, logger)
    if pars.filename:
        logger.info('Saving expanded population to %s' % pars.filename)
        pop.vars().clear()
        pop.save(pars.filename)
