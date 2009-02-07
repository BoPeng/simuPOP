#!/usr/bin/env python
#
# simulation for Reich's paper:
#     On the allelic spectrum of human disease
# plus migration and complexity of the common disease
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedDate$
# $Rev$
# 
# Usage:
#     simuReich2001.py 
#
#

"""

Introduction
=============


Reich and Lander's 2001 paper "On the allelic spectrum of human disease"
proposed a population genetics framework to model the evolution of 
allelic spectra. The model is based on the fact that human population
grew quickly from around 10,000 to 6 billion in 18,000 -150,000 years.
Their analysis showed that at the founder population, both common and
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
from simuUtil import *

import os, sys, types, exceptions, math

# declare all options, getParam will use these info to get parameters
# from a tk/wx-based dialog, command line, options etc.
options = [
    {'arg': 'h',
     'longarg': 'help',
     'default': False, 
     'description': 'Print this usage message.',
     'allowedTypes': [types.NoneType, type(True)],
     'jump': -1                    # if -h is specified, ignore any other parameters.
    },
    {'longarg': 'numDSL=',
     'default': 1,
     'label': 'Number of DSL',
     'description': '''Number of disease susecptibility loci for the disease.''',
     'allowedTypes': [types.IntType],
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'initSpec=',
     'default': [0.9]+[0.02]*5,
     'label': 'Initial allelic spectrum',
     'description': '''Initial allelic spectrum for the disease. It will 
                be the same for all disease susceptibility loci of the disease.
                The first element should be the allele frequency of the wild type (allele 1)
                and the rest are frquencies of disease alleles. These frequencies should
                add up to one. You can also specify spectra for all DSL using a 2d list''',
     'allowedTypes': [types.ListType, types.TupleType],
    },
    {'longarg': 'selModel=',
     'default': ['recessive'],
     'label': 'Selection model(s)',
     'description': '''Given selection coef, for recessive model,
                the fitness values are [1, 1, 1-s] for genotype AA,Aa,aa (A
                is the wild type). For additive model, the fitness values are
                [1,1-s/2,1-s]. You can also specify a list of models for each DSL.''' ,
     'allowedTypes': [types.ListType, types.TupleType],
     'validate': simuOpt.valueListOf( simuOpt.valueOneOf(['recessive','additive'])),
    },
    {'longarg': 'selModelAllDSL=',
     'default': 'additive',
     'label': 'Multi-locus selection model',
     'description': '''Overall fitness values given fitness values for each DSL,
                fitness values are Prod(f_i) for multiplicative model and
                1-Sum(1-f_i) for additive model. You can also choose customized. In this
                case, you will have to provide your own fitness table in selCoef.
                For example, for a 2-locus model, you need to provide 
                    [w11,w12,w13,w21,w22,w23,w31,w32,w33] 
                where, for example, w21 is the fitness value for genotype AaBB.''',
     'allowedTypes': [types.StringType],
     'chooseOneOf': [ 'additive', 'multiplicative', 'customized']
    }, 
    {'longarg': 'selCoef=',
     'default': [0.02],
     'label': 'Selection coefficient(s)',
     'description': '''Selection coef (s) for the disease.
                Fitness values will be [1,1,1-s] for the receissive model and
                [1,1-s/2,1-s] for the additive model. You can specify
                different values for each DSL. In the case of customized selection
                model, this coefficient is a fitness table, NOT selection coefficient(s).''',
     'allowedTypes': [types.ListType, types.TupleType],
     'validate':    simuOpt.valueListOf( simuOpt.valueGT(0.))
    },
    {'longarg': 'mutaModel=',
     'default': 'k-allele',
     'label': 'Mutation model',
     'allowedTypes': [types.StringType],
     'description': '''Microsatellite markers are mutated using    
                symmetric stepwise mutation wile SNP markers are mutaed
                using a 2-allele model (kam) DSL are not mutated unless in disease
                introduction stage.''',
     'validate':    simuOpt.valueOneOf(['k-allele', 'stepwise']),
     'chooseOneOf': ['k-allele', 'stepwise']
    },
    {'longarg': 'maxAllele=',
     'default': 255,
     'label': 'Max allele at each DSL',
     'allowedTypes': [types.IntType],
     'description': '''Maximum allele state at each DSL. For stepwise
                model, 255 is more than enough. For k-allele model, you may
                want to have more than 255 alleles especially when expected 
                effective allele numbers are big. Note that using large maxAllele
                to mimic infinite allele model will slow down the program 
                significantly.
                NOTE: very large maxAllele will significantly slow down the simulation
                due to the creation of large array as population variable.''',
     'validate':    simuOpt.valueGT(1),
    }, 
    {'longarg': 'mutaRate=',
     'default': [0.0001],
     'label': 'Mutation rate(s)',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Mutation rate for all DSL. Can be different for each DSL.''',
     'validate':    simuOpt.valueListOf( simuOpt.valueBetween(0,1))
    },
    {'longarg': 'initSize=',
     'default': 10000,
     'label': 'Initial population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial population size. This size will be maintained
                till the end of burnin stage''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'finalSize=',
     'default': 1000000,
     'label': 'Final population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Ending population size (after expansion.',
     'validate':    simuOpt.valueGT(0)
    }, 
    {'longarg': 'burnin=',
     'default': 2000,
     'label': 'Length of burn-in stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of generations of the burn in stage.',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'noMigrGen=',
     'default': 400,
     'label': 'Length of no-migration stage',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Number of generations when migration is zero. This stage
                is used to build up population structure.''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'mixingGen=',
     'default': 100,
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
                Choose between instant, linear and exponential''',
     'chooseOneOf': ['exponential', 'linear', 'instant'],
    },
    {'longarg': 'numSubPop=',
     'default': 1,
     'label': 'Number of split subpops',
     'allowedTypes': [types.IntType],
     'description': 'Number of subpopulations to be split into after burnin stage.',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'migrModel=',
     'default': 'none',
     'label': 'Migration model',
     'allowedTypes': [types.StringType],
     'description': '''Migration model. Choose between stepping stone and island.    
                A stepping stone model will be a circular model.''',
     'validate':    simuOpt.valueOneOf(['island', 'stepping stone', 'none']),
     'chooseOneOf': ['stepping stone', 'island', 'none']
    }, 
    {'longarg': 'migrRate=',
     'default': 0.05,
     'label': 'Migration rate',
     'description': '''Migration rate during mixing stage. 
                Island or circular stepping stone migration model can be used. ''',
     'allowedTypes': [types.FloatType, types.IntType],
     'validate':    simuOpt.valueBetween(0,1)
    },
     {'longarg': 'update=',
     'default': 10,
     'label': 'Update figure every # gen',
     'allowedTypes': [types.IntType],
     'description': '''Update figure every some generation.''',
     'validate':    simuOpt.valueGE(1)
    },    
    {'longarg': 'dispPlot=',
     'default': True,
     'label': 'Display plot?',
     'allowedTypes': [types.BooleanType],
     'description': 'If false, do not disply figure.',
    },
    {'longarg': 'saveAt=',
     'default': [x*100 for x in range(1,11)],
     'label': 'Save figure at generation',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''At these generations, figures will be saved, instead of displayed.
                Note that:
                    1): generations at which figure is not updated will not be saved. (update parameter)
                    2); figures will not be displayed at these generations
                    3): file name will be cdcvXXX.eps
                ''',
     'validate':    simuOpt.valueListOf(simuOpt.valueGE(0))
    },
    {'longarg': 'savePop=',
     'default': '',
     'allowedTypes': [types.StringType], 
     'description': 'save final population in a file. This feature is only available from command line.',
    },
    {'longarg': 'resume=',
     'default': '',
     'allowedTypes': [types.StringType], 
     'description': 'resume from a saved population. This option is only available from command line.'
    },
    {'longarg': 'resumeAtGen=',
     'default': 0,
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'when resume is not empty. Start from this generation. This option is only available from command line.'
    },
    {'longarg': 'dryrun',
     'default': False,
     'allowedTypes': [types.IntType],
     'validate':    simuOpt.valueOneOf([True, False]),
     'description':    'Only display how simulation will perform.'
     # do not save to config, do not prompt, so this appeared to be an undocumented option.
    },
    {'longarg': 'name=',
     'default': 'cdcv',
     'allowedTypes': [types.StringType],
     'label': 'Name of the simulation',
     'description': 'Base name for configuration (.cfg) log file (.log) and figures (.eps)'
    },
    {'arg': 'v',
     'longarg': 'verbose',
     'default': False,
     'allowedTypes': [types.NoneType, types.IntType],
     'description': 'Verbose mode.'
    }
]

def getOptions(details=__doc__):
    ''' get options from options structure,
        if this module is imported, instead of run directly,
        user can specify parameter in some other way.
    '''
    #
    # get all parameters, __doc__ is used for help info
    opt = simuOpt.simuOpt(options, 
        'This program simulates the evolution of a mono or polygenic disease using a three-stage\n' +
        'evolutionary scenario. The allelic spectrum at all disease susceptibility loci are recorded\n' +
        'and plotted (if R/RPy is available).\n', details)
    if not opt.getParam(nCol=2):
        sys.exit(1)
    # -h or --help
    if opt.help:    
        print opt.usage()
        sys.exit(0)
    #
    # --name
    if opt.name is None: 
        print "Name of the simulation is invalid"
        sys.exit(1)
    else:
        if not os.path.isdir(opt.name):
            os.makedirs(opt.name)
        opt.save(os.path.join(opt.name, opt.name + '.cfg'))
    #
    # --verbose or -v 
    if opt.verbose:                 # verbose
        print opt.listOptions
    # return the rest of the parameters
    return opt.asList()[1:-1]

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
        num = list(v.alleleNum[d][1:])
        # number of alleles
        numAllele.append(len(num) - num.count(0) )
        # effective number of alleles
        effNumAllele.append(overallFreq[-1]*overallFreq[-1]/sum( [x*x for x in v.alleleFreq[d][1:] ]))
        #effNumAllele.append(1./sum( [x*x for x in pop.dvars().alleleFreq[d][1:] ])-1)
        # count percentage of alleles derived from alleles before expansion
        allNum = sum(num)
        ancNum = 0
        if len(allelesBeforeExpansion) > 0:
            for al in allelesBeforeExpansion[d]:
                try:
                    ancNum += num[al-1]
                except:
                    pass
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
    " swtich from old-style python plotter to new style, use pyOperator"
    (numDSL, saveAt, highest, plot, plotLabel, name) = param
    # use global logOutput handle
    # this is less efficient but make sure we can see partial results
    logOutput = open(os.path.join(name, name+'.log'), "a")
    logOutput.write("%d\t" % pop.gen())
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
    NeHist[0].append(pop.gen())
    FHist[0].append(pop.gen())
    for d in range(numDSL):
        NeHist[d+1].append( effNumAllele[d] )
        FHist[d+1].append( overallFreq[d] )
    FMax = max( FMax, max( overallFreq) )
    NeMax = max( NeMax, max( effNumAllele))
    # now, variable like perc holds results for the whole populations.
    # plot something:
    if pop.gen() in saveAt:
        print "Saving figure in %s%d.eps (instead of displaying it)" % (name,pop.gen())
        r.postscript(file=os.path.join(name, '%s%d.eps' % (name, pop.gen())), width=6, height=8 )
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
    r.ltitle('Generation %d, PopSize %d' % (pop.gen(), pop.popSize()),
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
    if pop.gen() in saveAt:
        r.dev_off()
    return True


def simuCDCV(numDSL, initSpec, selModel,
        selModelAllDSL, selCoef, mutaModel, maxAllele, mutaRate, 
        initSize, finalSize, burnin, noMigrGen,
        mixingGen, growth, numSubPop, migrModel, migrRate,
        update, dispPlot, saveAt, savePop, resume,
        resumeAtGen, name, dryrun):
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
        simu = simulator(     
            population(size=incFunc(0), loci=[1]*(numDSL),
                maxAllele = maxAllele, infoFields=['fitness']),    
            randomMating(newSubPopSizeFunc=incFunc)
        )
    else:
        try:
            print "Resuming simulation from file ", resume, " at generation ", resumeAtGen
            pop = LoadPopulation(resume)
            simu = simulator(pop, randomMating(newSubPopSizeFunc=incFunc))
            simu.setGen(resumeAtGen)
        except exceptions.Exception, e:
            print "Can not resume from population "+ resume + ". Aborting."
            raise e

    # determine mutation etc
    if mutaModel == 'k-allele':
        mutation = kamMutator(rate=mutaRate, loci=range(numDSL), maxAllele=maxAllele)
    else:
        mutation = smmMutator(rate=mutaRate, loci=range(numDSL), maxAllele=maxAllele)
    # determine selection
    #
    if selModelAllDSL == 'customized':
        selection = maSelector( loci=range(numDSL), fitness=selCoef, wildtype=[0] )
    else:
        sel = []
        for d in range(numDSL):
            if selModel[d] == 'recessive':
                sel.append(maSelector(locus=d, fitness=[1,1,1-selCoef[d]], wildtype=[0]))
            else: 
                sel.append(maSelector(locus=d, fitness=[1,1-selCoef[d]/2.,1-selCoef[d]], wildtype=[0]))
        # now, the whole selector
        if selModelAllDSL == 'additive':
            selection = mlSelector(sel, mode=SEL_Additive)
        elif selModelAllDSL == 'multiplicative':
            selection = mlSelector(sel, mode=SEL_Multiplicative)
    # migration
    if numSubPop == 1 or migrModel == 'none':
        # no migration
        migration = noneOp()
    else:
        if migrModel == 'island':
            migration = migrator(MigrIslandRates(migrRate, numSubPop), begin=mixing)
        else:
            migration = migrator(MigrSteppingStoneRates(migrRate, numSubPop, circular=True), begin=mixing)
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
    simu.evolve(                            # start evolution
        preOps=
            # initialize DSL 
            [initByFreq(atLoci=[x], alleleFreq=initSpec[x]) for x in range(numDSL)] +
            [pyExec('allelesBeforeExpansion=[]')],
        ops=[         
            # report population size, for monitoring purpose only
            # count allele frequencies at both loci
            stat(popSize=True, alleleFreq=range(numDSL), step=update),
            # report generation and popsize
            pyEval(r"'%d\t%d\t%f\n' % (gen, popSize, alleleFreq[0][0])", step=50),
            #
            # record alleles before expansion, used to count percentage of alleles derived
            # from before expansion.
            pyExec( '''for i in range(%d):
                allelesBeforeExpansion.append([])
                for a in range(1,len(alleleNum[i])):
                    if alleleNum[i][a] != 0:
                        allelesBeforeExpansion[i].append(a)
                print "Ancestral alleles before expansion: ", allelesBeforeExpansion[i]''' % \
                numDSL, at=[burnin]),
            #
            splitSubPop(0, proportions=[1./numSubPop]*numSubPop, at=[burnin]),
            # mutate
            mutation,
            # selection
            selection,
            # migration
            migration, 
            # visualizer
            pyOperator(func=PlotSpectra, param=(numDSL, saveAt, 50, dispPlot, plotLabel, name), step=update ),
            # monitor execution time
            ticToc(step=100),
            ## pause at any user key input (for presentation purpose)
            ## pause(stopOnKeyStroke=1)
            savePopulation(outputExpr='"%s/%s-%%d.txt"%% gen' % (name, name), step=3000),
        ],
        end=end,
        dryrun = dryrun
    )
    #
    if savePop != '':
        simu.population(0).savePopulation(savePop)


if __name__ == '__main__':
    # get parameters
    allParam = getOptions()
    #
    # unpack parameters
    (numDSL, initSpecTmp,
        selModelTmp, selModelAllDSL, selCoefTmp, mutaModel, maxAllele, mutaRateTmp, 
        initSize, finalSize, burnin, noMigrGen, mixingGen, growth, numSubPop, 
        migrModel, migrRate, update, dispPlot, saveAt, savePop, resume,
        resumeAtGen, dryrun, name) = allParam

    if dispPlot:
        try:
            from simuRPy import *
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
    if type(initSpecTmp[0]) == types.FloatType:
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
        resume, resumeAtGen, name, dryrun)
    #
    # write report? What kind?
    # raw_input("Press any key to close this program and the R plot")        
        


