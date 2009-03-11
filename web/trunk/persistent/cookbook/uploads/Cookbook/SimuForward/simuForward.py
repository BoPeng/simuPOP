#!/usr/bin/env python
#
# $LastChangedDate: 2006-11-08 15:52:33 -0600 (Wed, 08 Nov 2006) $
# $Rev: 529 $
#

"""
Evolve a population forward in time, introduce a disease mutant
at given generation.
"""

import simuOpt
# load simuPOP libraries
from simuPOP import *
from simuUtil import MigrSteppingStoneRates, MigrIslandRates
import os, sys, exceptions, types, math

options = [
    {'separator': 'Genotype structure:'},    
    {'longarg': 'numChrom=',
     'default': 10,
     'label': 'Number of chromosomes',
     'description': 'Number of chromosomes.',
     'allowedTypes': [types.IntType],
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'numLoci=',
     'default': 20,
     'label': 'Number of loci on each chrom',
     'description': '''Number of loci on each chromosome, current there 
             only equal number of markers on each chromosome is supported.''',
     'allowedTypes': [types.IntType],
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'markerType=',
     'default': 'microsatellite',
     'label': 'Marker type',
     'description': '''Type of markers. Can be microsatellite or SNP.
                Microsatellite markers will be mutated using a symmetric
                stepwise model while SNP markers will be mutated using a
                2-state model. Mutation rate should be much smaller for
                SNP markers than that of microsatellite markers.''',
     'allowedTypes': [types.StringType],
     'validate':    simuOpt.valueOneOf(['microsatellite', 'SNP']),
     'chooseOneOf': ['microsatellite', 'SNP']
    },
    {'longarg': 'DSL=',
     'default': [45, 85, 125],
     'label': 'DSL after marker (0-indexed)',
     'description': '''A list of loci *after* a marker. For example, 
                35 means a disease locus after the 16th marker on chromosome 2,
                (if numChrom=numLoci=20). The number of DSL is important since
                it determines the complexity of the disease.
                A single number is allowed and implies a simple disease with
                one disease locus.''',
     'allowedTypes': [types.TupleType, types.ListType],
     'validate':    simuOpt.valueListOf( simuOpt.valueGT(0) )
    },
    {'longarg': 'DSLLoc=',
     'default': [.5, .5, .5],
     'label': 'DSL location between markers',
     'description': '''A list of loci location between two markers.
                Since all disease loci will be *between* equal spaced markers,
                the location should be between 0 and 1. A single value is acceptable
                as the location of all DSL.''',
     'allowedTypes': [types.TupleType, types.ListType],
     'validate':    simuOpt.valueListOf( simuOpt.valueBetween(0,1))
    },
    #
    #
    {'separator': 'Demographic model:'},
    {'longarg': 'initSize=',
     'default': 10000,
     'label': 'Initial population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial population size. This size will be maintained
                till the end of burnin stage''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'endingSize=',
     'default': 200000,
     'label': 'Final population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Final population size after population expansion.',
     'validate':    simuOpt.valueGT(0)
    }, 
    {'longarg': 'growthModel=',
     'default': 'exponential',
     'label': 'Population growth model',
     'description': '''How population is grown from initSize to endingSize.
                Choose between linear and exponential. Population expansion
                starts from after burninGen''',
     'chooseOneOf': ['exponential', 'linear'],
    },
    {'longarg': 'burninGen=',
     'default': 3000,
     'label': 'Length of burnin stage',
     'allowedTypes': [types.IntType],
     'description': '''Length of the burnin stage, also the generation to introduce mutants,
                population starts to expand after burnin stage.''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'introLen=',
     'default': 50,
     'label': 'Length of mutant introduction stage',
     'allowedTypes': [types.IntType],
     'description': '''Disease mutants are introduced after burnin, for introLen
                generations. This should be done *before* splitGen, since otherwise
                all mutants will be introduced to subpopulation 0.''',
     'validate':    simuOpt.valueGT(0),
    },
    {'longarg': 'splitGen=',
     'default': 5000,
     'label': 'When to split population',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''At which generation to split the population. 
                The population will start to grow after disease introduction stage but will
                not split till this generation. Note that if the disease is 
                introduced after this stage, it will be in one of subpopulations.
     ''',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'mixingGen=',
     'default': 8000,
     'label': 'When to mix population',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''At which generation to start mixing (allow migration.
                This number should be greater than or equal to split gen.''',
     'validate':    simuOpt.valueGE(0)
    },
    {'longarg': 'endingGen=',
     'default': 10000,
     'label': 'Ending generation number',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''At which generation to stop the simulation.
                This is the total generation number.''',
     'validate':    simuOpt.valueGE(0)
    },
    #
    #
    {'separator': 'Diease introduction:'},
    {'longarg': 'introSel=',
     'default': 0,
     'label': 'Selection coef of mutants',
     'allowedTypes': [types.FloatType, types.IntType],
     'description': '''Selection coef of mutants during disease-introduction stage.
                Individuals with these mutants will have fitness values 1,1-s/2,1-s for genotype
                AA,Aa,aa (a is mutatnt). The deault value is -0.5 meaning strong positive 
                selection on these mutants. You can set this value to zero to disable such
                artificial boost of allele frequency.''',
    },
    {'longarg': 'minAlleleFreq=',
     'default': 0.05,
     'label': 'minimal Allele Frequency',
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''Mininal allele frequencies required for each DSL,
                The simulation will restart if allele frequency is lower than
                this number. Can be given as a number (for all DSL), or a list.''',
     'validate':    simuOpt.valueListOf(simuOpt.valueBetween(0,1))
    },
    {'longarg': 'maxAlleleFreq=',
     'default': 0.20,
     'label': 'maximum Allele Frequency',
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''Maximum allele frequencies required for each DSL,
                The simulation will restart if allele frequency is greater than
                this number.''',
     'validate':    simuOpt.valueListOf(simuOpt.valueBetween(0,1))
    },
    #
    #
    {'separator': 'Migration parameters:'},    
    {'longarg': 'numSubPop=',
     'default': 1,
     'label': 'Number of subpopulations to split',
     'allowedTypes': [types.IntType],
     'description': 'Number of subpopulations to be split into after burnin stage.',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'migrModel=',
     'default': 'stepping stone',
     'useDefault': True,
     'label': 'Migration model',
     'allowedTypes': [types.StringType],
     'description': '''Migration model. Choose between stepping stone (circular),
                island and none. ''',
     'validate':    simuOpt.valueOneOf(['island', 'stepping stone', 'none']),
     'chooseOneOf': ['stepping stone', 'island', 'none']
    }, 
    {'longarg': 'migrRate=',
     'default': 0.0001,
     'label': 'Migration rate',
     'description': '''Migration rate during mixing stage. 
                A circular stepping stone migration model will be used. ''',
     'allowedTypes': [types.IntType, types.FloatType],
     'validate':    simuOpt.valueBetween(0,1)
    },
    #
    #
    {'separator': 'Disease model:'},
    {'longarg': 'fitness=',
     'default': [1, 1.0001, 1.0002],
     'label': 'Fitness of genotype AA,Aa,aa',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Fitness of genotype, can be:
                f1, f2, f3: if one DSL, the fitness for genotype AA, Aa and aa
                f1, f2, f3: if multiple DSL, the same fitness for each locus
                [a1, a2, a3, b1, b2, b3, ...] if selMultiLocusModel = 'additive' 
                    or multiplicative, fitness at each locus. The overall fitness
                    is determined by selMultiLocusModel
                [a1, a2, a3, b1, b2, b3, c1, c2, c3, ...] and selMultiLocusModel = interaction.
                    For example, in the 2-DSL case, the numbers are (by row)
                        BB Bb bb
                    AA  a1 a2 a3
                    Aa  b1 b2 b3
                    aa  c1 c2 c3
                3^n numbers are needed for n DSL.
        ''',
     'validate':    simuOpt.valueListOf(simuOpt.valueGE(0.)),
    },
    {'longarg': 'selMultiLocusModel=',
    'default': 'none',
    'label': 'Multi-locus selection model',
    'description': '''Model of overall fitness value given fitness values for each DSL.
                multiplicative: f =  Prod(f_i)
                additive: f = 1-Sum(1-f_i)
                interaction: the intepretation of fitness parameter is different.
                    see fitness.
                ''',
    'allowedTypes': [types.StringType],
    'chooseOneOf': [ 'additive', 'multiplicative', 'interaction', 'none']
    },
    #
    #
    {'separator': 'Evolutionary parameters:'},
    {'longarg': 'mutaRate=',
     'default': 0.0001,
     'label': 'Mutation rate',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': '''Microsatellite markers are mutated using    
                symmetric stepwise mutation wile SNP markers are mutaed
                using a 2-allele model (kam) DSL are not mutated unless in disease
                introduction stage.''',
     'validate':    simuOpt.valueBetween(0,1)
    },
    {'longarg': 'recRate=',
     'default': [0.0005],
     'label': 'Recombination rate',
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''If a number is given, it is the recombination rate between
                 adjacent markers. If a list is given, it should be the recombination
                 rate between all markers and DSL. For example, if you have two chromosome
                 with 5 markers, and with DSL after 3. The marker/DSL layout is
                     0 1 2 3 x 5        6 7 8 9 10    
                 The specified recombination rate should be those between 
                     0-1, 1-2, 2-3, 3-x, x-5, 6-7, 7-8, 8-9, 9-10.
                 Note that the distance between 3-x, x-5 is smaller than distance
                 between markers. If the list is too long, remember, that simuPOP dislog 
                 accept any python expression like [0.0005]*25 + [0.0009]*30.
     ''',
     'validate': simuOpt.valueListOf(simuOpt.valueBetween(0,1))
    },
    #
    {'separator': 'Final population preparation:'},
    {'longarg': 'savedGen=',
     'default': 2,
     'useDefault': True,
     'label': 'Generations to save',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': simuOpt.valueBetween(1, 3),
     'description': '''How many generations to save in the final population. 1 means no parental 
                generations, and 3 means grandparent, parent and current generations. Default
                is two, which is good for sampling schemes like affected sibpair sampling.''',
    },
    {'longarg': 'numOffspring=',
     'default': 2,
     'useDefault': True,
     'allowedTypes': [types.FloatType, types.IntType, types.LongType],
     'label': 'Number of offspring per mating', 
     'description': '''Number of offspring in these last generations. 2 is good for affected sibpair. More
                is needed for the sampling of large pedigrees. The value can be the number of offspring
                per mating when numOffMode=constant, but can also be parameters for a distribution. Details
                see description of numOffMode.''',
     'validate': simuOpt.valueGT(0),
    },
    {'longarg': 'numOffMode=',
     'default': 'constant',
     'useDefault': True,
     'label': 'Mode to determine number of offspring', 
     'chooseOneOf': [ 'constant', 'geometric'],
     'allowedTypes': [types.StringType],
     'description': '''Two ways to determine number of offspring. If constant, numOffspring will be used
                to generate the same number of offspring per mating event (in the final generations). If
                'geometric', P(k) = p*(1-p)^(k-1) where p=numOffspring. 0<p<1. The mean number of offspring
                is 1/p. This mode can be used when you need some large pedigrees.''',
    },
    #
    # 
    {'separator': 'Miscellaneous:'},
    {'longarg': 'savePop=',
     'default': [],
     'useDefault': True,
     'label': 'Save population at generations',
     'allowedTypes': [types.ListType, types.TupleType],
     'validate': simuOpt.valueListOf(simuOpt.valueGE(0)),
     'description': '''A list of generations at which populations are saved, 
                default to []. It can be used for intermediate analysis, but simulations can
                not be resumed from this population due to the re-generation of allele 
                frequency trajectories.''',
    },
    {'longarg': 'name=',
     'default': 'simu',
     'allowedTypes': [types.StringType],
     'label': 'Simulation name',
     'description': '''Name of simulation, files saved will be 
                    name + '.log': statistics output
                    name + '.cfg': configuration
                    name + .bin/txt/xml: saved popuation''',
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

    
def outputStatistics(pop, args): 
    ''' This function will be working with a pyOperator
        to output statistics. Many parameters will be passed,
        packed as a tuple.
    
     We need to output
     1. Fst before and after migration
     2. allele frequency at DSL 
     3. Mean observed heterogeneity at the marker loci
    '''
    # unwrap parameter
    (introLen, splitGen, mixingGen, endGen, outfile) = args
    # 
    gen = pop.gen()
    # see how long the simulation has been running
    if gen == introLen:
        print "Start population growth\t\t"
    elif gen == splitGen:
        print "Start no-migration stage\t\t"
    elif gen == mixingGen:
        print "Start mixing\t\t"
    elif gen == endGen:
        print "End of simulation. \n"
    #
    #### preparations 
    # number of loci when not counting DSL
    numLoci = pop.dvars().numLoci
    # non DSL loci, for convenience
    nonDSL = range( pop.totNumLoci() )
    DSL = pop.dvars().DSL
    for loc in DSL:
        nonDSL.remove(loc)
    #
    #@@@ NOW, output statistics
    # append to output file.
    output = open(outfile, 'a')
    # first, calculate LD and other statistics
    Stat(pop, alleleFreq=DSL, Fst = nonDSL, heteroFreq = range(pop.totNumLoci()))
    # output D', allele frequency at split, mixing and endGen
    print >> output, "Average Fst estimated from non-DSL at gen %d: %.4f \n" % (gen, pop.dvars().AvgFst)
    print >> output, "\n\nAllele frequencies\nall\t",
    for d in DSL:
        print >> output, '%.4f ' % (1. - pop.dvars().alleleFreq[d][0]),
    for sp in range(pop.numSubPop()):
        print >> output, "\n%d\t" % sp,
        for d in DSL:
            print >> output, '%.4f ' % (1. - pop.dvars(sp).alleleFreq[d][0]),
    print >> output, "\n",
    # hetero frequency
    AvgHetero = 0
    for d in range(pop.totNumLoci()):
        AvgHetero += pop.dvars().heteroFreq[d][0]
    AvgHetero /= pop.totNumLoci()
    # save to pop
    pop.dvars().AvgHetero = AvgHetero
    # output it
    print >> output, '\nAverage counted heterozygosity is %.4f.\n' % AvgHetero
    output.close()
    return True


def dynaAdvSelector(pop, param):
    ''' This selector will apply advantage/purifying selection to DSl according
        to allele frequency at each DSL. minAlleleFreq and maxAlleleFreq
        are stored in pop.dvars().
    '''
    (minAlleleFreq, maxAlleleFreq, DSL, introSel) = param;
    # get allele frequencies
    Stat(pop, alleleFreq=DSL)
    # gives 1,1.25,1.5 to promote allele if freq < lower cound
    # gives 1,0.9,0.8 to select agsinst DSL with freq > upper
    sel = []
    freq = pop.dvars().alleleFreq
    # print +- symbol for each DSL to visualize how frequencies are manipulated
    for i in range(len(DSL)):
        # positive selection (promote allele)
        if 1-freq[DSL[i]][0] < minAlleleFreq[i]:
            print '+',
            sel.append( maSelector(locus=DSL[i], wildtype=[0], fitness=[1,1-introSel/2., 1-introSel]) )
        # negative selection (select against allele)
        elif 1-freq[DSL[i]][0] > maxAlleleFreq[i]:
            print '-',
            sel.append( maSelector(locus=DSL[i], wildtype=[0], fitness=[1,0.9,0.8]) )
        # encourage slightly towards upper bound
        else:
            print ' ',
            sel.append( maSelector(locus=DSL[i], wildtype=[0], fitness=[1,1.002,1.004]) )
        # apply multi-locus selector, note that this operator will only
        # set a variable fitness in pop, actual selection happens during mating.
        if len(sel ) > 0:    # need adjustment (needed if 'else' part is empty)
            pop.turnOffSelection()
            MlSelect(pop, sel, mode=SEL_Multiplicative)
    print ' '.join(['%.4f' % (1-freq[x][0]) for x in DSL])
    return True


def simuForward(numChrom, numLoci, markerType, DSLafter, DSLdist, 
        initSize, endingSize, growthModel, 
        burninGen, introLen, splitGen, mixingGen, endingGen, 
        introSel, minAlleleFreq, maxAlleleFreq, 
        numSubPop, migrModel, migrRate,
        fitness, mlSelModel, 
        mutaRate, recRate, savedGen, numOffspring, numOffMode,
        savePop, filename):
    ''' run a simulation of complex disease with given parameters. 
    '''    
    ###
    ### CHECK PARAMETERS
    ###
    if initSize < len(DSLafter):
        raise exceptions.ValueError("Initial population size is too small. (Less than number of DSL)")
    # expand .5 -> [0.5, 0.5...]
    if len(DSLdist) == 1:
        DSLdist = DSLdist * len(DSLafter)
    if len( DSLafter ) != len(DSLdist):
        raise exceptions.ValueError("Please specify DSL distance for each DSL.")
    numDSL = len(DSLafter)
    if burninGen  + introLen > splitGen or splitGen > mixingGen or splitGen > endingGen:
        raise exceptions.ValueError("Generations should in the order of burninGen, introLen, splitGen, mixingGen and ending")
    # fitness
    if mlSelModel == 'none':
        fitness = []
    elif mlSelModel == 'interaction':
        if numDSL == 1:
            raise exceptions.ValueError("Interaction model can only be used with more than one DSL");
        if len(fitness) != 3**numDSL:
            raise exceptions.ValueError("Please specify 3^n fitness values for n DSL");
    else:
        if fitness == []:    # neutral process
            fitness = [1,1,1]*numDSL
        else:
            # for a single DSL
            if len(fitness) == 3:
                fitness = fitness*numDSL
            elif len(fitness) != numDSL*3:
                raise exceptions.ValueError("Please specify fitness for each DSL")
    # range
    if len(minAlleleFreq) == 1:
        minAlleleFreq = minAlleleFreq * len(DSLafter)
    elif len(minAlleleFreq) != len(DSLafter):
        raise exceptions.ValueError("min allele frequency should be specified for each DSL")
    if len(maxAlleleFreq) == 1:
        maxAlleleFreq = maxAlleleFreq * len(DSLafter)
    elif len(maxAlleleFreq) != len(DSLafter):
        raise exceptions.ValueError("max allele frequency should be specified for each DSL")

    # number of offspring
    if numOffMode == 'geometric' and numOffspring > 1:
        raise exceptions.ValueError("numOffspring is p for a geometric distribution when numOffMode='geometric'. It should be elss than 1.") 
    ###
    ### translate numChrom, numLoci, DSLafterLoci to
    ### loci, lociPos, DSL, nonDSL in the usual index
    ### (count DSL along with markers)
    ###
    if markerType == 'microsatellite':
        maxAle = 99                                     # max allele
    else:
        maxAle = 1                                        # SNP (0 and 1)
    numDSL = len(DSLafter)
    if len(DSLdist) != numDSL:
        raise exceptions.ValueError("--DSL and --DSLloc has different length.")
    # real indices of DSL
    DSL = list(DSLafter)
    for i in range(0,numDSL):
        DSL[i] += i+1
    # adjust loci and lociPos
    loci = [0]*numChrom
    lociPos = []
    i = 0    # current absolute locus 
    j = 0    # current DSL
    for ch in range(0, numChrom):
        lociPos.append([])
        for loc in range(0,numLoci):
            # DSL after original loci indices
            if j < numDSL and i == DSLafter[j]:
                loci[ch] += 2
                lociPos[ch].append(loc+1)
                lociPos[ch].append(loc+1 + DSLdist[j])
                j += 1
            else:
                loci[ch] += 1
                lociPos[ch].append(loc+1)
            i += 1
    # non-DSL loci, for convenience
    nonDSL = range(0, sum(loci))
    for loc in DSL:
        nonDSL.remove(loc)
    ###
    ### initialization 
    ###
    if maxAle > 1:    # Not SNP
        preOperators = [
            # initialize all loci with 10 haplotypes
            initByValue(value=[[x]*sum(loci) for x in range(50, 60)],
                proportions=[.1]*10), 
            # and then init DSL with all wild type alleles
            initByValue([0]*len(DSL), loci=DSL)
        ]
    else: # SNP
        preOperators = [
            # initialize all loci with two haplotypes (000, 111)
            initByValue(value=[[x]*sum(loci) for x in [0,1] ],
                proportions=[.5]*2),
            # and then init DSL with all wild type alleles
            initByValue([0]*len(DSL), loci=DSL)
        ]
    ###
    ### mutation, start from gen 0,
    ###
    if maxAle > 1:    # Not SNP
        # symmetric mutation model for microsatellite
        mutator = smmMutator(rate=mutaRate, maxAllele=maxAle, loci=nonDSL)
    else:
        # k-allele model for mutation of SNP
        mutator = kamMutator(rate=mutaRate, maxAllele=1, loci=nonDSL)
    ###
    ### Recombination
    ###
    if len(recRate) == 1: 
        rec = recombinator(intensity=recRate[0])
    else:
        # rec after ...
        if len(recRate) != sum(loci) - numChrom:
            raise exceptions.ValueError("Recombination rate specification is wrong. Please see help file")
        al = []
        start = 0
        for ch in range(numChrom):
            al.extend([x+start for x in range(loci(ch))] )
            start += loci[ch]
        rec = recombinator(rate=recRate, afterLoci=al)
    ###
    ### skip burnin stage when the simulation is restarted.
    ###
    burnin_pop = 'burnin.bin'
    # clear burnin pop from last time.
    if os.path.isfile(burnin_pop):
        os.remove(burnin_pop)
    ###
    ### output progress
    ###
    operators = [
        mutator, 
        rec, 
        savePopulation(burnin_pop, at=[burninGen-1]),
        stat(alleleFreq = DSL, popSize =True, end = burninGen, step = 100),
        stat(alleleFreq = DSL, popSize =True, begin = burninGen, end = burninGen + introLen),
        stat(alleleFreq = DSL, popSize = True, begin = burninGen + introLen, step = 100),
        # output to screen
        pyEval( expr=r'"%d(%d): "%(gen, popSize) + " ".join(["%.5f"%(1-alleleFreq[x][0]) for x in DSL])+"\n"',
            step=100, end = burninGen),
        pyEval( expr=r'"%d(%d): "%(gen, popSize) + " ".join(["%.5f"%(1-alleleFreq[x][0]) for x in DSL])+"\n"',
            step=100, begin = burninGen + introLen),
        # output to file (append)
        #pyEval( expr=r'"%d %d " %(gen, popSize) + " ".join(["%.5f"%(1-alleleFreq[x][0]) for x in DSL])+"\n"',
        #    output='>>>'+filename+'.traj')
    ]
    ###
    ### introduction of disease mutants
    ###
    ###                                 endingGen
    # need several conditional point mutator to
    # introduce disease. It does not matter that we introduce
    # mutant specifically to individuel i since individuals
    # are unordered
    #
    for i in range(numDSL):
        # this operator literally says: if there is no DS allele,
        # introduce one. Note that operator stat has to be called
        # before this one.
        operators.append( 
            ifElse("alleleFreq[%d][0]==1." % DSL[i],
                pointMutator(loci=[DSL[i]], toAllele=1, inds=[i]),
            begin=burninGen, end = burninGen + introLen) 
        )
    # optionally, the mutants will be given some selective pressure
    # if introSel < 0, mutants will have selective advantage and will
    # reach high allele frequency quickly.
    if introSel != 0:
        operators.append(
            pyOperator(func=dynaAdvSelector, param = (minAlleleFreq, maxAlleleFreq, DSL, introSel), 
                begin=burninGen, end = burninGen + introLen) )
    #
    operators.extend([
        # the simulation will stop if the disease allele frequencies
        # are not within range at the end of this stage
        terminateIf("True in [(1.-alleleFreq[DSL[i]][0] < minAlleleFreq[i] or 1.-alleleFreq[DSL[i]][0] > maxAlleleFreq[i]) for i in range(len(DSL))]",
            at = [burninGen + introLen]), 
        # terminate if any disease allele get lost
        terminateIf("True in [alleleFreq[x][0] == 1. for x in DSL]",
            begin = burninGen + introLen),
    ])  
    ### 
    ### split to subpopulations
    ### 
    if numSubPop > 1:
        operators.append( 
            splitSubPop(0, proportions=[1./numSubPop]*numSubPop, at=[splitGen]),
        )
    ###
    ### selection 
    ###
    if mlSelModel in ['additive', 'multiplicative']:
        mlSelModel = {'additive':SEL_Additive, 
            'multiplicative':SEL_Multiplicative}[mlSelModel]
        operators.append( mlSelector(
            # with five multiple-allele selector as parameter
            [ maSelector(locus=DSL[x], wildtype=[0], 
                fitness=[fitness[3*x],fitness[3*x+1],fitness[3*x+2]]) for x in range(len(DSL)) ],
            mode=mlSelModel, begin= burninGen + introLen),
        )
    elif mlSelModel == 'interaction':
        # multi-allele selector can handle multiple DSL case
        operators.append( maSelector(loci=DSL, fitness=fitness, wildtype=[0], begin = burninGen + introLen) )
    ###
    ### migration
    ###
    if numSubPop > 1 and migrModel == 'island' and migrRate > 0:
        operators.append( migrator(MigrIslandRates(migrRate, numSubPop),
            mode=MigrByProbability, begin=mixingGen) )
    elif numSubPop > 1 and migrModel == 'stepping stone' and migrRate > 0:
        operators.append( migrator(MigrSteppingStoneRates(migrRate, numSubPop, 
            circular=True),    mode=MigrByProbability, begin=mixingGen) )
    ###
    ### output statistics, track performance
    ###
    operators.extend([
        pyOperator(func=outputStatistics, 
            param = (burninGen, splitGen, mixingGen, endingGen, filename+'.log'),
            at = [burninGen, splitGen, mixingGen, endingGen ] ), 
        # show elapsed time
        ticToc(at=[burninGen, burninGen + introLen, splitGen, mixingGen, endingGen]),
        ]
    )
    if len(savePop) > 0:
        # save population at given generations
        operators.append(savePopulation(outputExpr='os.path.join(name, "%s_%d.pop" % (name, gen))', 
            at=savePop))
    ###
    ###  demographic model
    ###
    if growthModel == 'linear':
        popSizeFunc = LinearExpansion(initSize, endingSize, endingGen,
            burninGen, splitGen, numSubPop)
    elif growthModel == 'exponential':
        popSizeFunc = ExponentialExpansion(initSize, endingSize, endingGen,
            burninGen, splitGen, numSubPop)
    else:
        raise exceptions.ValueError("Growth model can be one of linear and exponential. Given " + growth)   
    ### 
    ### prepare mating scheme
    ###
    # clear log file if it exists
    try:
        os.remove(filename+'.log')
    except:
        pass
    #
    # may need to run several times to
    # get a valid population (requirements on allele frequency
    fixedCount = 1
    while(True):
        # create a simulator
        if os.path.isfile(burnin_pop):
            pop = LoadPopulation(burnin_pop)
            simu = simulator(pop, randomMating(), rep=1)
            # skip the burnin stage.
            print "Loading %s and skip burnin stage (gen=%d)" % (burnin_pop, burninGen)
            simu.setGen(burninGen)
            simu.evolve(ops = operators, gen = endingGen)
        else:
            pop = population(size=popSizeFunc(0), ploidy=2,
                    loci = loci,
                    lociPos = lociPos,
    				infoFields = ['fitness', 'father_idx', 'mother_idx'])
            # save DSL info, some operators will use it.
            pop.dvars().DSL = DSL
            pop.dvars().minAlleleFreq = minAlleleFreq
            pop.dvars().maxAlleleFreq = maxAlleleFreq
            pop.dvars().numLoci = numLoci
            simu = simulator(pop, randomMating(), rep=1)
            #
            simu.evolve( preOps = preOperators, ops = operators, gen=endingGen - savedGen)
        if simu.gen() != endingGen - savedGen:
            print "Population restarted at gen ", simu.gen(), endingGen-savedGen
            print "Overall fixed population ", fixedCount
            print "Allelefreq ", ( "%.3f " * numDSL + "\n\n") % \
                tuple([1-simu.dvars(0).alleleFreq[x][0] for x in DSL]) 
            fixedCount += 1
        else:
            # prepare for the last several generations
            # change mating scheme to random mating.
            #
            # save saveGen-1 ancestral generations
            simu.population(0).setAncestralDepth(savedGen-1)
            simu.population(0).addInfoFields(['father_idx', 'mother_idx'])
            operators.append(parentsTagger())
            simu.setMatingScheme(
                randomMating(
                    newSubPopSizeFunc=popSizeFunc,            # demographic model
                    numOffspring = numOffspring,
                    mode = {'constant': MATE_NumOffspring, 
                            'geometric': MATE_GeometricDistribution
                           }[numOffMode]
                )
            )
            # different evolution method for the last few generations
            simu.evolve(ops=operators, gen=endingGen)
            # succeed save information
            pop = simu.population(0)
            # we want to save info on how this population is generated.
            # This is not required but is a good practise
            pop.dvars().DSLafter = DSLafter
            pop.dvars().DSLdist = DSLdist
            pop.dvars().initSize = initSize
            pop.dvars().endingSize = endingSize
            pop.dvars().burninGen = burninGen
            pop.dvars().burninGen = burninGen
            pop.dvars().splitGen = splitGen
            pop.dvars().mixingGen = mixingGen
            pop.dvars().endingGen = endingGen
            pop.dvars().numSubPop = numSubPop
            pop.dvars().mutaRate = mutaRate
            pop.dvars().mutaModel = "symmetric stepwise"
            pop.dvars().migrRate = migrRate
            pop.dvars().migrModel = "circular stepping stone"
            pop.dvars().recRate = recRate
            pop.dvars().numOffspring = numOffspring
            pop.dvars().numOffMode = numOffMode
            print "Saving population to " + filename + '\n'
            #TurnOnDebug(DBG_UTILITY)
            simu.population(0).savePopulation(filename)
            break;
    return True

if __name__ == '__main__':
    par = simuOpt.simuOpt(options, 
        '''This program simulates the evolution of a complex common disease, subject 
     to the impact of mutation, migration, recombination and population size change. 
     Click 'help' for more information about the evolutionary scenario.''', __doc__)
    if not par.getParam():
        sys.exit(1)
    #
    if not os.path.isdir(par.name):
        os.mkdir(par.name)
    par.saveConfig(os.path.join(par.name, par.name + '.cfg'))
    # unpack options
    (numChrom, numLoci, markerType, DSLafter, DSLdist, 
        initSize, endingSize, growthModel, 
        burninGen, introLen, splitGen, mixingGen, endingGen, 
        introSel, minAlleleFreq, maxAlleleFreq, 
        numSubPop, migrModel, migrRate,
        fitness, selMultiLocusModel,
        mutaRate, recRate, savedGen, numOffspring, numOffMode,
        savePop, name) = par.asList()
    #
    ################## RUN THE SIMULATION ###############
    simuForward(numChrom, numLoci, markerType, DSLafter, DSLdist, 
        initSize, endingSize, growthModel, 
        burninGen, introLen, splitGen, mixingGen, endingGen, 
        introSel, minAlleleFreq, maxAlleleFreq, 
        numSubPop, migrModel, migrRate,  
        fitness, selMultiLocusModel, 
        mutaRate, recRate, savedGen, numOffspring, numOffMode,
        savePop, os.path.join(name, name))
    print "Done!"
