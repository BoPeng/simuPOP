#!/usr/bin/env python
#
# simulation for Reich's paper:
#     On the allelic spectrum of human disease
# plus migration and complexity of the common disease
#
# Bo Peng (bpeng@mdanderson.org)
#
# $LastChangedDate: 2009-03-11 11:55:59 -0500 (Wed, 11 Mar 2009) $
# $Rev: 2530 $
# 
#

"""

Introduction
=============


Reich and Lander's 2001 paper "On the allelic spectrum of human disease"
proposed a population genetics framework to model the evolution of 
allelic spectra. The model is based on the fact that human Population
grew quickly from around 10,000 to 6 billion in 18,000 -150,000 years.
Their analysis showed that at the founder Population, both common and
rare diseases have simple spectra. After the sudden expansion of 
population size, the allelic spectra of simple diseases become complex;
while those of complex diseases remained simple.

This script simulates the evolution of a disease using a more complicated
version of this evolutionary process. The allelic spectrum is observed, 
recorded and plotted (if R/Rpy is available).

"""
import simuOpt
simuOpt.setOptions(alleleType='long')
from simuPOP import *

import os, sys, types, exceptions, math

# declare all options, getParam will use these info to get parameters
# from a tk/wx-based dialog, command line, options etc.
options = [
    {'name': 'numDSL',
     'default': 1,
     'label': 'Number of DSL',
     'description': '''Number of disease susecptibility loci for the disease.''',
     'type': [int],
     'validator':    simuOpt.valueGT(0)
    },
    {'name': 'initSpec',
     'default': [0.9]+[0.02]*5,
     'label': 'Initial allelic spectrum',
     'description': '''Initial allelic spectrum for the disease. It will 
                be the same for all disease susceptibility loci of the disease.
                The first element should be the allele frequency of the wild type (allele 1)
                and the rest are frquencies of disease alleles. These frequencies should
                add up to one. You can also specify spectra for all DSL using a 2d list''',
     'type': [types.ListType, types.TupleType],
    },
    {'name': 'selModel',
     'default': ['recessive'],
     'label': 'Selection model(s)',
     'description': '''Given selection coef, for recessive model,
                the fitness values are [1, 1, 1-s] for genotype AA,Aa,aa (A
                is the wild type). For additive model, the fitness values are
                [1,1-s/2,1-s]. You can also specify a list of models for each DSL.''' ,
     'type': [types.ListType, types.TupleType],
     'validator': simuOpt.valueListOf( simuOpt.valueOneOf(['recessive','additive'])),
    },
    {'name': 'selModelAllDSL',
     'default': 'additive',
     'label': 'Multi-locus selection model',
     'description': '''Overall fitness values given fitness values for each DSL,
                fitness values are Prod(f_i) for multiplicative model and
                1-Sum(1-f_i) for additive model. You can also choose customized. In this
                case, you will have to provide your own fitness table in selCoef.
                For example, for a 2-locus model, you need to provide 
                    [w11,w12,w13,w21,w22,w23,w31,w32,w33] 
                where, for example, w21 is the fitness value for genotype AaBB.''',
     'type': [str],
     'chooseOneOf': [ 'additive', 'multiplicative', 'customized']
    }, 
    {'name': 'selCoef',
     'default': [0.02],
     'label': 'Selection coefficient(s)',
     'description': '''Selection coef (s) for the disease.
                Fitness values will be [1,1,1-s] for the receissive model and
                [1,1-s/2,1-s] for the additive model. You can specify
                different values for each DSL. In the case of customized selection
                model, this coefficient is a fitness table, NOT selection coefficient(s).''',
     'type': [types.ListType, types.TupleType],
     'validator':    simuOpt.valueListOf( simuOpt.valueGT(0.))
    },
    {'name': 'mutaModel',
     'default': 'k-allele',
     'label': 'Mutation model',
     'type': [str],
     'description': '''Microsatellite markers are mutated using    
                symmetric stepwise mutation wile SNP markers are mutaed
                using a 2-allele model (kam) DSL are not mutated unless in disease
                introduction stage.''',
     'validator':    simuOpt.valueOneOf(['k-allele', 'stepwise']),
     'chooseOneOf': ['k-allele', 'stepwise']
    },
    {'name': 'maxAllele',
     'default': 255,
     'label': 'Max allele at each DSL',
     'type': [int],
     'description': '''Maximum allele state at each DSL. For stepwise
                model, 255 is more than enough. For k-allele model, you may
                want to have more than 255 alleles especially when expected 
                effective allele numbers are big. Note that using large maxAllele
                to mimic infinite allele model will slow down the program 
                significantly.
                NOTE: very large maxAllele will significantly slow down the simulation
                due to the creation of large array as population variable.''',
     'validator':    simuOpt.valueGT(1),
    }, 
    {'name': 'mutaRate',
     'default': [0.0001],
     'label': 'Mutation rate(s)',
     'type': [types.ListType, types.TupleType],
     'description': '''Mutation rate for all DSL. Can be different for each DSL.''',
     'validator':    simuOpt.valueListOf( simuOpt.valueBetween(0,1))
    },
    {'name': 'initSize',
     'default': 10000,
     'label': 'Initial population size',
     'type': [int, long],
     'description': '''Initial population size. This size will be maintained
                till the end of burnin stage''',
     'validator':    simuOpt.valueGT(0)
    },
    {'name': 'finalSize',
     'default': 1000000,
     'label': 'Final population size',
     'type': [int, long],
     'description': 'Ending population size (after expansion.',
     'validator':    simuOpt.valueGT(0)
    }, 
    {'name': 'burnin',
     'default': 2000,
     'label': 'Length of burn-in stage',
     'type': [int],
     'description': 'Number of generations of the burn in stage.',
     'validator':    simuOpt.valueGT(0)
    },
    {'name': 'noMigrGen',
     'default': 400,
     'label': 'Length of no-migration stage',
     'type': [int, long],
     'description': '''Number of generations when migration is zero. This stage
                is used to build up population structure.''',
     'validator':    simuOpt.valueGT(0)
    },
    {'name': 'mixingGen',
     'default': 100,
     'label': 'Length of mixing stage',
     'type': [int, long],
     'description': '''Number of generations when migration is present. This stage
                will mix individuals from subpopulations using an circular stepping stone
                migration model.''',
     'validator':    simuOpt.valueGT(0)
    },    
    {'name': 'growth',
     'default': 'exponential',
     'label': 'Population growth model',
     'description': '''How Population is grown from initSize to finalSize.
                Choose between instant, linear and exponential''',
     'chooseOneOf': ['exponential', 'linear', 'instant'],
    },
    {'name': 'numSubPop',
     'default': 1,
     'label': 'Number of split subpops',
     'type': [int],
     'description': 'Number of subpopulations to be split into after burnin stage.',
     'validator':    simuOpt.valueGT(0)
    },
    {'name': 'migrModel',
     'default': 'none',
     'label': 'Migration model',
     'type': [str],
     'description': '''Migration model. Choose between stepping stone and island.    
                A stepping stone model will be a circular model.''',
     'validator':    simuOpt.valueOneOf(['island', 'stepping stone', 'none']),
     'chooseOneOf': ['stepping stone', 'island', 'none']
    }, 
    {'name': 'migrRate',
     'default': 0.05,
     'label': 'Migration rate',
     'description': '''Migration rate during mixing stage. 
                Island or circular stepping stone migration model can be used. ''',
     'type': [float, int],
     'validator':    simuOpt.valueBetween(0,1)
    },
     {'name': 'update',
     'default': 100,
     'label': 'Update figure every # gen',
     'type': [int],
     'description': '''Update figure every some generation.''',
     'validator':    simuOpt.valueGE(1)
    },    
    {'name': 'dispPlot',
     'default': True,
     'label': 'Display plot?',
     'type': [types.BooleanType],
     'description': 'If false, do not disply figure.',
    },
    {'name': 'saveAt',
     'default': [x*100 for x in range(1,11)],
     'label': 'Save figure at generation',
     'type': [types.ListType, types.TupleType],
     'description': '''At these generations, figures will be saved, instead of displayed.
                Note that:
                    1): generations at which figure is not updated will not be saved. (update parameter)
                    2); figures will not be displayed at these generations
                    3): file name will be cdcvXXX.eps
                ''',
     'validator':    simuOpt.valueListOf(simuOpt.valueGE(0))
    },
    {'name': 'savePop',
     'default': '',
     'type': [str], 
     'description': 'save final population in a file. This feature is only available from command line.',
    },
    {'name': 'resume',
     'default': '',
     'type': [str], 
     'description': 'resume from a saved population. This option is only available from command line.'
    },
    {'name': 'resumeAtGen',
     'default': 0,
     'type': [int, long],
     'description': 'when resume is not empty. Start from this generation. This option is only available from command line.'
    },
    {'name': 'name',
     'default': 'cdcv',
     'type': [str],
     'label': 'Name of the simulation',
     'description': 'Base name for configuration (.cfg) log file (.log) and figures (.eps)'
    },
]

def LinearExpansion(initSize, endSize, end, burnin=0, split=0, numSubPop=1, bottleneckGen=-1, bottleneckSize=0):
    '''
    Linearly expand population size from intiSize to endSize
    after burnin, split the population at generation split.
    '''
    inc = (endSize-initSize)/float(end-burnin)
    def func(gen):
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
    def func(gen):
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
    def func(gen):
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

# 
# calculate statistics from a populations
# 0. percentage
# 1. number of alleles
# 2. effective number of alleles
# 3. overall disease allele frequency
# 4. percentage of most frequent alleles among disease alleles
# 5. percentage of five most frequent alleles among disease alleles
# 6. percentage of alleles from ancestral population
# 
def getStats(v, highest, numDSL, allelesBeforeExpansion):
    # the following are statistics for each DSL
    perc = []                                 # percentage
    numAllele = []                        # number of alleles
    effNumAllele = []                 # effective number of alleles
    overallFreq = []                    # size of disease
    percMostCommon = []             # percentage of most common
    perc5MostCommon = []            # percentage of five most common
    percAncestralAllele = []    # percentage of alleles from before expansion
    for d in range(numDSL):
        # overall DSA frequency
        overallFreq.append( 1 - v.alleleFreq[d][0])
        if overallFreq[-1] == 0. : # no disease allele
            numAllele.append(0)
            effNumAllele.append(0)
            perc.append([])
            percMostCommon.append(0)
            perc5MostCommon.append(0)
            percAncestralAllele.append(0)
            continue
        # number of alleles
        numAllele.append(len(v.alleleFreq[d]))
        # effective number of alleles
        daf = [v.alleleFreq[d][x] for x in v.alleleFreq[d].keys() if x != 0]
        effNumAllele.append(overallFreq[-1]*overallFreq[-1]/sum( [x*x for x in daf]))
        #effNumAllele.append(1./sum( [x*x for x in pop.dvars().alleleFreq[d][1:] ])-1)
        # count percentage of alleles derived from alleles before expansion
        num = [v.alleleNum[d][x] for x in v.alleleNum[d].keys() if x != 0]
        allNum = sum(num)
        ancNum = 0
        if len(allelesBeforeExpansion) > 0:
            for al in allelesBeforeExpansion[d]:
                ancNum += v.alleleNum[d][al]
        percAncestralAllele.append(ancNum*1./allNum)
        # percentage of each one
        sumAllele = sum(num)*1.
        percNonZero = [0]*min(highest, len(num))
        # python 2.4 can do sort(reverse=True)
        num.sort()
        num.reverse()
        for n in range(min(len(num), highest)):
            percNonZero[n] = num[n]/sumAllele
        perc.append(percNonZero)
        percMostCommon.append(perc[-1][0]*100)
        perc5MostCommon.append(sum(perc[-1][0:5])*100)
    return( [perc, numAllele, effNumAllele, overallFreq, percMostCommon, perc5MostCommon, percAncestralAllele] )
 
     
# This function plot the allelic spectrum
# using given allele frequencies
#
def PlotSpectra(pop, param):
    " swtich from old-style python plotter to new style, use PyOperator"
    (numDSL, saveAt, highest, plot, plotLabel, name) = param
    # use global logOutput handle
    # this is less efficient but make sure we can see partial results
    logOutput = open(os.path.join(name, name+'.log'), "a")
    logOutput.write("%d\t" % pop.dvars().gen)
    for sp in range(pop.numSubPop()):
        # unpack result
        [perc, numAllele, effNumAllele, overallFreq, percMostCommon, perc5MostCommon, percAncestralAllele] \
            = getStats(pop.dvars(sp), highest, numDSL, pop.dvars().allelesBeforeExpansion)
        for d in range(numDSL):
            logOutput.write( '%.5f\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t' % (effNumAllele[d], \
                numAllele[d], overallFreq[d], percMostCommon[d], perc5MostCommon[d], percAncestralAllele[d]))
    # if there are subpop
    if pop.numSubPop() > 1: # then write overall state
        [perc, numAllele, effNumAllele, overallFreq, percMostCommon, perc5MostCommon, percAncestralAllele] \
            = getStats(pop.dvars(), highest, numDSL, pop.dvars().allelesBeforeExpansion)
        for d in range(numDSL):
            logOutput.write( '%.5f\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t' % (effNumAllele[d], \
                numAllele[d], overallFreq[d], percMostCommon[d], perc5MostCommon[d], percAncestralAllele[d]))
    logOutput.write("\n")
    logOutput.close()
    # 
    # if do not plot, return
    if not plot:
        return True
    # record global ne history
    global NeHist, NeMax, FHist, FMax
    NeHist[0].append(pop.dvars().gen)
    FHist[0].append(pop.dvars().gen)
    for d in range(numDSL):
        NeHist[d+1].append( effNumAllele[d] )
        FHist[d+1].append( overallFreq[d] )
    FMax = max( FMax, max( overallFreq) )
    NeMax = max( NeMax, max( effNumAllele))
    # now, variable like perc holds results for the whole populations.
    # plot something:
    if pop.dvars().gen in saveAt:
        print "Saving figure in %s%d.eps (instead of displaying it)" % (name,pop.dvars().gen)
        r.postscript(file=os.path.join(name, '%s%d.eps' % (name, pop.dvars().gen)), width=6, height=8 )
    #
    # set no conversion mode to save execution time
    set_default_mode(NO_CONVERSION)
    # a two column layout
    r.layout( r.matrix(
            [1, 1,        # head
             2, 3,        # label: 
             4, 5,        # history of ne
             6, 7,        # history of f
             8, 9] +    # label:
             range(10,10+numDSL*2)                 # spectra at each DSL
             + [10+numDSL*2, 10+numDSL*2], # final legend
        byrow=True,
        ncol=2),
        width=[1/6.,5/6.],
        height=['1 cm', '1 cm', str(1./(numDSL+2)), str(1./(numDSL+2)), '1 cm'] + [str(1./(numDSL+2))]*numDSL+['2 cm'])
    # par...
    r.par(mar=[0.1]*4)
    r.par(oma=range(1,5))
    # windows 1
    r.ltitle('Generation %d, PopSize %d' % (pop.dvars().gen, pop.popSize()),
        backcolor='white', forecolor='darkblue', cex=2.5)
    # windows 2
    r.ltitle(''),
    # windows 3
    r.ltitle('History of effective number of alleles'),
    # adjust font size for statistics
    if numDSL < 3:
        cx = 2
    elif numDSL <= 10:
        cx = 1.3
    else:
        cx = 1
    # windows 4
    r.ltitle('Effective\nnumber\n of \nalleles', backcolor='white', forecolor='blue', cex=cx)
    # windows 5: history
    r.plot(NeHist[0], NeHist[1], ylim=[0, NeMax*5./4.], col=1,
        axes = False, type='n', xlab='', ylab='')
    for d in range(0, numDSL):
        r.lines(NeHist[0], NeHist[d+1], type='l', col=d+1)
    r.axis(1)
    r.axis(2, r.pretty([0, NeMax], n=2))
    r.box()
    # windows 6
    r.ltitle('Total\nDisease\nAllele\nFreq.', backcolor='white', forecolor='blue', cex=cx)
    # windows 7
    r.plot(FHist[0], FHist[1], ylim=[0, FMax*5./4.], col=1,
        axes = False, type='n', xlab='', ylab='')
    for d in range(0, numDSL):
        r.lines(FHist[0], FHist[d+1], type='l', col=d+1)
    r.axis(2, r.pretty([0, FMax], n=2))
    r.box()
    # window 8
    r.ltitle('DSL')
    # windows 9
    maxPerc = 0
    for d in range(1, numDSL):
        if len(perc[d]) > 0 and perc[d][0] > maxPerc:
            maxPerc = perc[d][0]
    r.ltitle('Allelic Spectrum (Max Perc=%.1f%%)' % (maxPerc*100))
    # windows 10 ~ 10 + 2*numDSL 
    for d in range(0, numDSL):
        # left labels
        r.ltitle('%.2f, %d\n%.3f%%\n%.1f%%(1)\n%.1f%%(5)\n%.1f%%' \
        % (effNumAllele[d], numAllele[d], overallFreq[d]*100,    \
             percMostCommon[d], perc5MostCommon[d], percAncestralAllele[d]*100),
        backcolor='white', forecolor='blue', cex=cx)    
        # 4 - all DSL of disease
        if len(perc[d]) == 0:
            r.plot(0, type='n', axes=False)
        else:
            r.plot(perc[d], type='n', axes=False, ylim=[0,maxPerc], xlab='', ylab='',
                xlim=[0,10], cex=5)
            l = len(perc[d])
            r.rect(range(l), [0]*l, range(1,l+1), perc[d], col=d+1, border='white')
            r.text(4, maxPerc*2./3., plotLabel[d], cex=1.5, adj=0)
            r.axis(1)
            r.axis(2, r.pretty([0, maxPerc], n=2), 
                r.paste(r.pretty([0, 100*maxPerc], n=2)))
            r.box()
    # the last legend
    r.ltitle('''Effective number of disease alleles, number of alleles
total disease allele frequency
Percent of 1 (5) most allele among all disease alleles''',
        cex=1.5, forecolor='blue')
    if pop.dvars().gen in saveAt:
        r.dev_off()
    return True


def simuCDCV(numDSL, initSpec, selModel,
        selModelAllDSL, selCoef, mutaModel, maxAllele, mutaRate, 
        initSize, finalSize, burnin, noMigrGen,
        mixingGen, growth, numSubPop, migrModel, migrRate,
        update, dispPlot, saveAt, savePop, resume,
        resumeAtGen, name):
    ''' parameters are self-expanary. See help info for
        detailed simulation scheme. '''
    # generations
    mixing = burnin + noMigrGen
    end = mixing + mixingGen
    # pop size
    if growth == 'linear':
        incFunc = LinearExpansion(initSize, finalSize, end,
            burnin, burnin, numSubPop)
    elif growth == 'exponential':
        incFunc = ExponentialExpansion(initSize, finalSize, end,
            burnin, burnin, numSubPop)
    elif growth == 'instant':
        incFunc = InstantExpansion(initSize, finalSize, end,
            burnin, burnin, numSubPop)
    #
    # create a simulator, if not in resume mode
    if resume == '':
        pop = Population(size=incFunc(0), loci=[1]*(numDSL), infoFields=['fitness'])
    else:
        try:
            print "Resuming simulation from file ", resume, " at generation ", resumeAtGen
            pop = LoadPopulation(resume)
        except exceptions.Exception, e:
            print "Can not resume from population "+ resume + ". Aborting."
            raise e

    # determine mutation etc
    if mutaModel == 'k-allele':
        mutation = KAlleleMutator(rates=mutaRate, loci=range(numDSL), k=maxAllele+1)
    else:
        mutation = StepwiseMutator(rates=mutaRate, loci=range(numDSL), maxAllele=maxAllele)
    # determine selection
    #
    if selModelAllDSL == 'customized':
        selection = MaSelector(loci=range(numDSL), fitness=selCoef, wildtype=[0] )
    else:
        sel = []
        for d in range(numDSL):
            if selModel[d] == 'recessive':
                sel.append(MaSelector(loci=d, fitness=[1,1,1-selCoef[d]], wildtype=[0]))
            else: 
                sel.append(MaSelector(loci=d, fitness=[1,1-selCoef[d]/2.,1-selCoef[d]], wildtype=[0]))
        # now, the whole selector
        if selModelAllDSL == 'additive':
            selection = MlSelector(sel, mode=ADDITIVE)
        elif selModelAllDSL == 'multiplicative':
            selection = MlSelector(sel, mode=MULTIPLICATIVE)
    # migration
    if numSubPop == 1 or migrModel == 'none':
        # no migration
        migration = NoneOp()
    else:
        if migrModel == 'island':
            migration = Migrator(MigrIslandRates(migrRate, numSubPop), begin=mixing)
        else:
            migration = Migrator(MigrSteppingStoneRates(migrRate, numSubPop, circular=True), begin=mixing)
    #            
    # prepare log file, if not in resume mode
    if resume == '':    # not resume
        logOutput = open(os.path.join(name, name+'.log'), 'w')
        logOutput.write("gen\t")
        for d in range(numDSL):
            logOutput.write( 'ne\tn\tf0\tp1\tp5\tanc\t')
        logOutput.write("\n")
        logOutput.close()
    # use global
    global NeHist, NeMax, FHist, FMax
    # determine plot label
    plotLabel = []
    for i in range(numDSL):
        if selModelAllDSL == 'customized':
            plotLabel.append('mu=%g, s=customized' % mutaRate[i])
        else:
            plotLabel.append('mu=%g, s=%g' % (mutaRate[i], selCoef[i]))
    # history of Ne, the first one is gen
    NeHist = []
    FHist = []
    for i in range(numDSL+1):
        NeHist.append( [] )
        FHist.append( [] )
    NeMax = 0
    FMax = 0
    # start evolution
    pop.evolve(                            # start evolution
        initOps=
            [InitSex()] + 
            # initialize DSL 
            [InitGenotype(loci=[x], freq=initSpec[x]) for x in range(numDSL)] +
            [PyExec('allelesBeforeExpansion=[]')],
        preOps = [
            #
            SplitSubPops(0, proportions=[1./numSubPop]*numSubPop, at=[burnin]),
            # mutate
            mutation,
            # selection
            selection,
            # migration
            migration, 
        ],
        matingScheme = RandomMating(subPopSize=incFunc),
        postOps = [         
            # report population size, for monitoring purpose only
            # count allele frequencies at both loci
            Stat(popSize=True, alleleFreq=range(numDSL), vars=['alleleFreq', 'alleleFreq_sp', 'alleleNum', 'alleleNum_sp'], step=update),
            # report generation and popsize
            PyEval(r"'%d\t%d\t%f\n' % (gen, popSize, alleleFreq[0][0])", step=50),
            #
            # record alleles before expansion, used to count percentage of alleles derived
            # from before expansion.
            PyExec( '''for i in range(%d):
                allelesBeforeExpansion.extend([alleleFreq[i].keys()])
                print "Ancestral alleles before expansion: ", allelesBeforeExpansion[i]''' % \
                numDSL, at=[burnin]),
            # visualizer
            PyOperator(func=PlotSpectra, param=(numDSL, saveAt, 50, dispPlot, plotLabel, name), step=update ),
            # monitor execution time
            TicToc(step=100),
            ## Pause at any user key input (for presentation purpose)
            ## Pause(stopOnKeyStroke=1)
            SavePopulation(output='!"%s/%s-%%d.txt"%% gen' % (name, name), step=3000),
        ],
        gen=end
    )
    #
    if savePop != '':
        pop.SavePopulation(savePop)


if __name__ == '__main__':
    # get parameters
    par = simuOpt.Params(options, 
        'This program simulates the evolution of a mono or polygenic disease using a three-stage\n' +
        'evolutionary scenario. The allelic spectrum at all disease susceptibility loci are recorded\n' +
        'and plotted (if R/RPy is available).', __doc__)
    if not par.getParam():
        sys.exit(1)
    #
    if not os.path.isdir(par.name):
        os.makedirs(par.name)
        par.saveConfig(os.path.join(par.name, par.name + '.cfg'))
    #
    # unpack parameters
    (numDSL, initSpecTmp,
        selModelTmp, selModelAllDSL, selCoefTmp, mutaModel, maxAllele, mutaRateTmp, 
        initSize, finalSize, burnin, noMigrGen, mixingGen, growth, numSubPop, 
        migrModel, migrRate, update, dispPlot, saveAt, savePop, resume,
        resumeAtGen, name) = par.asList()

    if dispPlot:
        try:
            from simuPOP.plotter import *
            from rpy import *
            r.library('lattice')
        except:
            print "RPy module not found. Can not view spectrum history."
            hasRPy = False
        else:
            hasRPy = True
    else:
        hasRPy = False
        
    # verify parameters, basically, all parameters will be a vector
    # for each DSL. Single values will be copied.
    initSpec = []
    if type(initSpecTmp) not in [types.ListType, types.TupleType]:
        raise exceptions.ValueError("Expecting a spectrum or a list of spectra")
    if type(initSpecTmp[0]) == float:
        initSpec = [initSpecTmp]*numDSL
    elif type(initSpecTmp[0]) not in [types.TupleType, types.ListType]:
        raise exceptions.ValueError("Expecting a spectrum or a list of spectra")
    elif len( initSpecTmp ) not in [1, numDSL]:
        raise exceptions.ValueError("Length of initSpec should equal to numDSL")
    elif len(initSpecTmp) == 1:
        initSpec = initSpecTmp*numDSL
    else:
        initSpec = initSpecTmp
    for spec in initSpec:
        if abs(sum(spec) - 1.0) > 1.e-7 :
            raise exceptions.ValueError("Allele frequencies should add up to one.")
    #
    if len(selModelTmp) == 1:
        selModel = selModelTmp*numDSL
    else:
        if len(selModelTmp) != numDSL:
            raise exceptions.ValueError("Number of models should match number of DSL")
        selModel = selModelTmp
    #
    if len(selCoefTmp) == 1:
        selCoef = selCoefTmp*numDSL
    else:
        if selModelAllDSL != 'customized' and len(selCoefTmp) != numDSL:
            raise exceptions.ValueError("Number of selection coef should match number of DSL")
        selCoef = selCoefTmp
    #
    if len(mutaRateTmp) == 1:
        mutaRate = mutaRateTmp*(numDSL)
    elif len(mutaRateTmp) == numDSL:
        mutaRate = mutaRateTmp
    else:
        raise exceptions.ValueError("Mutation rate, if more than one are given, should have length numDSL")
    #
    #define a function, copied from R homepage figure script
    # this will be used to draw labels of the R-based figures
    if hasRPy:
        r('''ltitle= function(x,backcolor="#e8c9c1",forecolor="darkred",cex=2,ypos=0.4){
            plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
            polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col=backcolor,border=NA)
            text(x=0,y=ypos,pos=4,cex=cex,labels=x,col=forecolor)}''')
    # run simulation
    simuCDCV(numDSL, initSpec,
        selModel, selModelAllDSL, selCoef, 
        mutaModel, maxAllele, mutaRate, 
        initSize, finalSize, burnin, noMigrGen,
        mixingGen, growth, numSubPop, migrModel, migrRate,
        update, dispPlot, saveAt, savePop,
        resume, resumeAtGen, name)
    #
    # write report? What kind?
    # raw_input("Press any key to close this program and the R plot")        
        


