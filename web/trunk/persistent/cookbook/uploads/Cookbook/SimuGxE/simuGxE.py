#!/usr/bin/env python
'''
This script simulates a population based on the HapMap dataset.
'''

from simuOpt import *
setOptions(alleleType='binary')
from simuPOP import *

import os, sys, math
from types import *
from exceptions import *
from simuPOP.utils import simuProgress, SaveCSV

# This does not yet exist
from simuCommon import *

HapMap_pops = ['CEU', 'YRI', 'JPT+CHB']

options = [
    param_name(default='Forward'),
    {'separator': 'Progress record and report'},
    {'longarg': 'step=',
     'default': 100,
     'label': 'Progress report interval',
     'useDefault': True,
     'allowedTypes': [IntType, LongType],
     'description': '''Gap between generations at which population statistics are
                calculated and reported. (This parameter is affected by --scale)'''
    },
    {'longarg': 'saveStep=',
     'default': 0,
     'label': 'Save population interval',
     'useDefault': True,
     'allowedTypes': [IntType, LongType],
     'description': '''Initial, expanded and admixed populations will be saved. This option
                allows you to save populations every --saveStep generations, starting
                from population expansion. If saveStep = 0 (default), no intermediate population
                will be saved. Otherwise, populations at generation 0, saveStep, 2*saveStep, ...
                will be saved as expand_xx.pop where xx is generation number.''',
    },
    {'longarg': 'saveName=',
     'default': 'expand_',
     'useDefault': True,
     'description': '''Prefix of the intermediately saved populations, relative to simulation path''',
     'allowedTypes': [StringType],
    },
    {'separator': 'Populations and markers to use'},
    param_HapMap_dir(),
    param_HapMap_pop(),
    param_markerList(),
    param_markerListCols(),
    param_chroms(default=[2]),
    param_numMarkers(default=[1000]),
    param_startPos(default=[30]),
    param_endPos(),
    param_minAF(),
    param_initName(default='initPop.pop'),
    #
    {'separator': 'Mutation, recombination, etc'},
    param_mutaRate(default=4e-6),
    param_recMap(),
    param_recIntensity(1e-9),
    param_fitness(),
    param_mlSelModel(),
    param_scale(default=10),
    {'separator': 'Population expansion'},
    {'longarg': 'initCopy=',
     'default': 20,
     'useDefault': True,
     'label': 'Initial propagation',
     'description': '''How to expand the initial small HapMap sample to
                 avoid quick loss of heterogenity. By default, each individual
                 is copied 10 times.''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'initGen=',
     'default': 20,
     'useDefault': True,
     'label': 'Generations to evolve',
     'description': '''Number of generations to evolve to get the seed
                population. The actual evolved population is scaled down by
                parameter --scale. (If scale==10, initGen=1000, the actually
                evolved generation is 100).''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'initSize=',
     'default': 5000,
     'useDefault': True,
     'label': 'Size of the seed population',
     'description': '''Size of the seed population. The default value is the recommended
                value when all hapmap populations are used (60+60+90)*20. You may want
                to reduce it according to the populations used.''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGE(100)
    },
    {'longarg': 'expandGen=',
     'default': 500,
     'useDefault': True,
     'label': 'Generations to expand',
     'description': '''Number of generations to evolve during the population
                expansion stage. The actual evolved population is scaled down by
                parameter --scale. (If scale==10, expandGen=1000, the actually
                evolved generation is 100).''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'expandSize=',
     'default': 50000,
     'useDefault': True,
     'label': 'Expanded population size',
     'description': '''Size of the expanded population. The default value if the recommended
                value when all hapmap populations are used (60+60+90)*200. You may want to
                reduce it according to the population used, or increase it if disease
                prevalence if low and insufficient cases are generated.''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGE(100)
    },
    # draw sample using a penetrance model
    param_model(''),
    param_cases(0),
    param_controls(1000),
    param_caseProp([]),
    param_controlProp([]),
    param_sampleName('sample'),
    param_saveFormats(['simuPOP'], chooseFrom=['simuPOP', 'csv']),
]



def adjustParam(par):
    #
    par.initFile = paramSetFile(par.name, par.initName)
    par.sampleFile = paramSetFile(par.name, par.sampleName + '.csv')
    par.samplePop = paramSetFile(par.name, par.sampleName + '.pop')
    par.chosenMarkerFile = paramSetFile(par.name, 'markers.lst')
    #
    if len(par.chroms) == 0:
        raise ValueError('Please specify one or more chromosomes')
    # in case that chrom is a tuple
    par.chroms = list(par.chroms)
    numChrom = len(par.chroms)
    par.numMarkers = paramExpandList(par.numMarkers, numChrom,
        'Please specify number of marker for each chromosome')
    par.startPos = paramExpandList(par.startPos, numChrom,
        'Wrong starting positions')
    par.endPos = paramExpandList(par.endPos, numChrom,
        'Wrong endinging positions')
    # this parameter does not need to be configurable.
    par.ldSampleSize = 200
    par.markerListFile = os.path.join(par.name, 'markers.lst')
    #
    # parameters for fitness...
    par.mlSelModel = {
        'additive': Additive,
        'multiplicative': Multiplicative,
        'interaction': 'interaction',
        'none': None
        }[par.mlSelModel]
    #
    par.mutaRate *= par.scale
    par.recIntensity *= par.scale
    par.step = int(par.step / par.scale)
    par.saveStep = int(par.saveStep / par.scale)
    par.initGen = int(par.initGen / par.scale)
    par.expandGen = int(par.expandGen / par.scale)
    if par.initGen == 0:
        raise ValueError('No initial stage. It is possible that your scale parameter' + \
            ' (%s) is larger than initGen (%s).' % (par.scale, par.initGen))
    if par.expandGen == 0:
        raise ValueError('No expansion stage. It is possible that your scale parameter ' + \
            ' (%s) is larger than expandGen (%s) ' % (par.scale, par.expandGen))


######################################################################

def expDemoFunc(N0, N1, N2, gen1, gen2):
    '''
    Return an exponential population expansion demographic function that has
    two stages of expansion.

    N0: a list of initial subpopulation sizes.
    N1: middle population size.
    N2: ending population size.
    
    gen1: generations to evolve in the slow-expansion stage.
    gen2: generations to evolve in the fast-expansion stage.
    '''
    if type(N1) in [IntType, LongType]:
        midSize = [int(N1*1.0*x/sum(N0)) for x in N0]
    elif len(N1) != len(N0):
        raise exceptions.ValueError("Number of subpopulations should be the same")
    else:
        midSize = N1
    if type(N2) in [IntType, LongType]:
        endSize = [int(N2*1.0*x/sum(midSize)) for x in midSize]
    elif len(N2) != len(N0):
        raise exceptions.ValueError("Number of subpopulations should be the same")
    else:
        endSize = N2
    #
    rate1 = [(math.log(midSize[x]) - math.log(N0[x]))/gen1 for x in range(len(N0))]
    rate2 = [(math.log(endSize[x]) - math.log(midSize[x]))/gen2 for x in range(len(N0))]
    def func(gen, oldSize=[]):
        if gen < gen1:
            return [int(N0[x]*math.exp(gen*rate1[x])) for x in range(len(N0))]
        else:
            return [int(midSize[x]*math.exp((gen-gen1)*rate2[x])) for x in range(len(N0))]
    return func


def effPopSize(func, gen):
    '''Estimate effective population size for a given demographic function'''
    nSP = len(func(0))
    Ne = [0]*nSP
    for i in range(gen):
        for j in range(nSP):
            Ne[j] += 1./func(i)[j]
    return [gen/x for x in Ne]


def writeTrajectory(popFunc, trajFunc, gen, file):
    'Save trajectory to a file'
    sz = popFunc(0)
    t = trajFunc(0)
    numSP = len(sz)
    numLoci = len(t) / len(sz)
    file = open(file, 'w')
    print >> file, 'gen, %s, %s' % (', '.join(['sp_%d' % x for x in range(numSP)]),
        ' ,'.join(['traj_loc%d_sp%d' % (x, y) for x in range(numSP) for y in range(numLoci)]))
    for g in range(gen):
        print >> file, '%d, %s, %s' % (g, ', '.join([str(x) for x in popFunc(g)]),
            ', '.join(['%.4f' % x for x in trajFunc(g)]))
    file.close()


def getOperators(pop, par, progress=False, savePop=False, vsp=False, mutation=False,
        recombination=False, selection=False):
    '''Return mutation and recombination operators'''
    preOps = []
    postOps = []
    duringOps = []
    if progress:
        # statistics calculation and display
        exp = ['gen %3d', 'size=%s']
        preGen = 'gen*scale'
        postGen = '(gen+1)*scale-1'
        var = ['%s', 'subPopSize']
        #if len(par.DSL) > 0:
        #    exp.append('alleleFreq=%s')
        #    var.append('", ".join(["%%.3f" %% alleleFreq[x][1] for x in DSL])')
        if vsp:
            exp.append('VSP=%s')
            var.append('virtualPopSize')
        keyGens = [par.initGen - 1, -1]
        postOps.extend([
            stat(popSize = True, Fst = range(pop.totNumLoci()),
                step = par.step, stage=PreMating),
            pyEval(r'"At the beginning of %s\n" %% (%s)' % (', '.join(exp), ', '.join(var) % preGen),
                step=par.step, stage=PreMating),
            stat(popSize = True, Fst = range(pop.totNumLoci()),
                at = keyGens),
            pyEval(r'"At the end of %s\n" %% (%s)' % (', '.join(exp), ', '.join(var) % postGen),
                at = keyGens)
        ])
    if savePop and par.saveStep > 0:
        postOps.extend([
            pyEval(r"'Saving current generation to %s%%d.pop\n' %% (gen*scale)" % par.saveName,
                step=par.saveStep, stage=PreMating),
            savePopulation(outputExpr="'%s/%s%%d.pop' %% (gen*scale)" % (par.name, par.saveName),
                step=par.saveStep, stage=PreMating),
        ])
    if mutation:
        preOps.append(kamMutator(rate=par.mutaRate, loci=range(pop.totNumLoci())))
    if recombination:
        if par.recMap == 'physical':
            print 'Scaled recombination at %.3f cM/Mb over %.2f Mb physical distance (first chromosome)' % \
                (par.recIntensity * 1e6, pop.lociDist(0, pop.numLoci(0)-1))
            duringOps.append(recombinator(intensity=par.recIntensity))
        else: # use map distance
            try:
                pos = [pop.dvars().genDist[pop.locusName(x)] for x in range(pop.totNumLoci())]
            except Exception,e:
                print e
                print 'Invalid or incomplete population variable genDist'
                print 'Please run loadHapMap again to set up genetic distance'
            loc = []
            rate = []
            for ch in range(pop.numChrom()):
                beg = pop.chromBegin(ch)
                end = pop.chromEnd(ch)
                loc.extend(range(beg, end - 1))
                rate.extend([(pos[x+1] - pos[x])*par.recIntensity for x in range(beg, end - 1)])
            print 'Scaled recombination at %.3f cM/Mb over %.2f centiMorgan genetic (%.2f Mb physical) distance (first chromosome)' % \
                (par.recIntensity*1e6, (pop.dvars().genDist[pop.locusName(pop.numLoci(0)-1)] - \
                    pop.dvars().genDist[pop.locusName(0)]), pop.lociDist(0, pop.numLoci(0)-1))
            # recombination rate at the end of each chromosome will be invalid
            # but this does not matter
            duringOps.append(recombinator(rate=rate, loci = loc))
    if selection:
        if par.mlSelModel in [Additive, Multiplicative]:
            preOps.append(mlSelector(
                # with five multiple-allele selector as parameter
                [ maSelector(locus=par.DSL_Idx[x], wildtype=[0],
                    fitness=[par.fitness[3*x], par.fitness[3*x+1], par.fitness[3*x+2]]) \
                        for x in range(len(par.DSL_Idx)) ],
                mode=par.mlSelModel))
        elif par.mlSelModel == 'interaction':
            # multi-allele selector can handle multiple DSL case
            preOps.append(maSelector(loci=par.DSL_Idx, fitness=par.fitness, wildtype=[0]))
    return {'preOps': preOps, 'duringOps': duringOps, 'postOps': postOps}


def createInitialPopulation(par, logger=None):
    '''Create an initial population, with parameters (from the par structure)
    HapMap_dir:     directory that stores hapmap data.
    chrom:          chromosomes to use
    markerList:     list of markers to use. It can be a map file of illumina or
        affymetrix chips.
    numMarkers:     number of markers per chromosome
    startPos:       starting position on each chromosome
    endPos:         ending position on each chromosome
    minAF:          minimal allele frequency
    minDist:        minimal distance between adjecent markers
    pops:           hapmap populations to use
    '''
    # load markers!
    for ch in par.chroms:
        if not os.path.isfile(os.path.join(par.HapMap_dir, 'HapMap_%s_chr%d.pop' % (par.HapMap_pop, ch))):
            raise ValueError('''Failed to load or download hapmap data for chromosome %d
                    Please copy script loadHapMap.py to the current directory, or add
                    path to this script to environmental variable$PYTHONPATH,
                    or run this script manually to download, import, and save HapMap
                    data in simuPOP format''' % ch)
    if par.markerList != '':
        names = getMarkersFromAnnotation(markerList=par.markerList,
            markerListCols=par.markerListCols, chroms=par.chroms,
            startPos=par.startPos, endPos=par.endPos,
            logger=logger)
    else:
        names = []
    #
    pop = getHapMapMarkers(HapMap_dir=par.HapMap_dir, names=names,
        chroms=par.chroms, HapMap_pops=[par.HapMap_pop], startPos=par.startPos,
        endPos=par.endPos, numMarkers=par.numMarkers, minAF=par.minAF,
        logger=logger)
    # if this population fine?
    if pop.numChrom() != len(par.chroms):
        raise ValueError('Something wrong. The population does not have enough chromosomes')
    # write marker and map file
    writeMarkerInfo(par.chosenMarkerFile, pop, logger=logger)
    #
    if par.initFile is not None and not os.path.isfile(par.initFile):
        if logger is not None:
            logger.info('Saving initial population to %s ' % par.initFile)
        pop.save(par.initFile)
    #
    return pop

caseCnt = 0
controlCnt = 0
totalCnt = 0

def filterInd(ind, par):
    '''
    Select individuals according to their affection status
    '''
    global caseCnt, controlCnt, totalCnt
    totalCnt += 1
    if par.model != '' and not par.affected(ind, par):
        return False
    if ind.affected() and caseCnt < par.cases:
        caseCnt += 1
        par.progress.update(caseCnt + controlCnt + 1)
        return True
    elif (not ind.affected()) and controlCnt < par.controls:
        controlCnt += 1
        par.progress.update(caseCnt + controlCnt + 1)
        return True
    return False
    

def simuForward(par, pop, logger=None):
    '''The main program'''
    #
    if par.model == '':
        if par.cases != 0:
            if logger is not None:
                logger.error('Cannot simulate cases without a disease model')
            sys.exit(1)
        if logger is not None:
            logger.info('No model is given so only controls will be simulated')
        par.DSL = []
        par.DSL_idx = []
        par.infoFields = []
    else:
        execfile(par.model + '.py', globals(), globals())
        par.DSL = getDSL(pop)
        if len(par.DSL) == 0:
            if logger is not None:
                logger.error('Could not locate any DSL')
            sys.exit(1)
        par.DSL_idx = pop.lociByNames(par.DSL)
        par.infoFields = getInfoFields()
        par.affected = affected
        for loc in par.DSL_idx:
            if logger is not None:
                logger.info('DSL: %s, chrom %s, index %d, pos: %df, Frequency of allele 1 = %.4f%%' % \
                    (pop.locusName(loc), pop.chromLocusPair(loc)[0], loc,
                pop.locusPos(loc), 100 * pop.dvars().alleleFreq[loc][1]))
    #
    pop.dvars().scale = par.scale
    newSize = [x*par.initCopy for x in pop.subPopSizes()]
    if logger is not None:
        logger.info('Propagating population to size %s' % newSize)
    pop.resize(newSize, propagate=True)
    #
    popSizeFunc = expDemoFunc(pop.subPopSizes(), par.initSize, par.expandSize, 
        par.initGen, par.expandGen)
    if logger is not None:
        logger.info('Estimated effective population size is %s' % \
            effPopSize(popSizeFunc, par.initGen + par.expandGen))
    #
    if par.model == '':
        if logger is not None:
            logger.info("Evolving the population freely...")
        simu = simulator(pop, randomMating(subPopSize = popSizeFunc))
        simu.evolve(
            initOps = initSex(),
            gen = par.initGen + par.expandGen,
            **getOperators(pop, par,
                progress=True,
                savePop=False,
                selection=False,
                mutation=True,
                recombination=True)
            )
        expandedPop = simu.extract(0)
    else:
        if logger is not None:
            logger.info("Using controlled random mating...")
        # define a trajectory function
        Stat(pop, alleleFreq=par.DSL_idx)
        currentFreq = []
        # in the order: LOC0: sp0, sp1, sp2, LOC1: sp0, sp1, sp2, ...
        for idx,loc in enumerate(par.DSL_idx):
            for sp in range(pop.numSubPop()):
                currentFreq.append(pop.dvars(sp).alleleFreq[loc][1])
        #
        traj = ForwardFreqTrajectory(
            curGen = 0,
            endGen = par.initGen + par.expandGen,
            curFreq = currentFreq,
            freq = [(max(0, x - 0.01), x + 0.01) for x in currentFreq],
            #fitness = par.fitness,
            #migrRate = par.backMigrRate,
            NtFunc = popSizeFunc,
            maxAttempts = 10000
        )
        if len(traj) == 0:
            raise SystemError('Failed to generated trajectory after 10000 attempts. '
                'This usually means that the demographic and genetic settings are '
                'very extreme which makes it very likely for an allele to reach designed'
                'allele frequency. Please adjust your parameters and try again.')
        # define a trajectory function
        def trajFunc(gen):
            return [x[gen] for x in traj]
        #for i in range(par.initGen + par.expandGen + 1):
        #    print trajFunc(i)
        # record trajectory
        if logger is not None:
            logger.info('Using controlled random mating with trajectory: %s' % traj)
        #
        simu = simulator(pop,
            controlledRandomMating(
                loci = par.DSL_idx,
                alleles = [1]*len(par.DSL_idx),
                freqFunc = trajFunc,
                subPopSize = popSizeFunc)
        )
        simu.evolve(
            initOps = initSex(),
            gen = par.initGen + par.expandGen,
            ** getOperators(pop, par,
                progress=True,
                savePop=False,
                selection=False,
                mutation=True,
                recombination=True
            )
        )
        expandedPop = simu.extract(0)
    if logger is not None:
        logger.info('Start drawing samples')
    # evolve one last generation to get samples
    expandedPop.mergeSubPops(range(expandedPop.numSubPop()))
    simu = simulator(expandedPop, randomMating(subPopSize=[par.cases + par.controls]))
    par.progress = simuProgress('Generating case (%d) and control (%d) sample' % \
        (par.cases, par.controls), par.cases + par.controls)
    simu.evolve(
        duringOps = pyOperator(func=filterInd, param=par, offspringOnly=True),
        gen = 1
    )
    par.progress.done()
    return simu.extract(0)


short_desc = '''This program simulates an admixed population based on
two or more HapMap populations. Please follow the intructions
of the help message to prepare HapMap population.'''

# determine which script to run.
if __name__ == '__main__':
    #
    # get all parameters
    par = simuOpt(options, short_desc, __doc__)
    # when user click cancel ...
    if not par.getParam():
       sys.exit(1)
    # simulation directory
    if not os.path.isdir(par.name):
        print 'Creating directory', par.name
        os.makedirs(par.name)
    if not os.path.isdir(par.name):
        raise SystemError('Can not create directory %s, exiting' % par.name)
    #
    logFile = os.path.join(par.name, 'simulation.log')
    logger = getLogger(logFile)
    #
    cfgFile = 'simulation.cfg'
    logger.info('Save configuration to %s' % cfgFile)
    # save current configuration
    par.saveConfig(os.path.join(par.name, cfgFile))
    adjustParam(par)
    #
    if os.path.isfile(par.initFile):
        logger.info('Loading existing initial population from file %s. '
            'Please remove this file if parameters are changed.' % par.initFile)
        init = LoadPopulation(par.initFile)
    else:
        logger.info('Creating initial population')
        init = createInitialPopulation(par, logger)
    sample = simuForward(par, init.clone(), logger)
    if par.sampleName not in [None, ''] and 'simuPOP' in par.saveFormats:
        logger.info('Saving sample to %s in simuPOP format' % par.samplePop)
        sample.vars().clear()
        sample.save(par.samplePop)
    # save in csv format
    def addMarker(markers):
        return markers[0] + markers[1] + 1
    if par.sampleName not in [None, ''] and 'csv' in par.saveFormats:
        logger.info('Saving sample to %s in csv format' % par.sampleFile)
        SaveCSV(sample, output=par.sampleFile, fields=['sex', 'affection'] + par.infoFields,
            combine=addMarker)
