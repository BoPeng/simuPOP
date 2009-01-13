#!/usr/bin/env python
#
# Purpose:    generate dataset for common complex disease 
#                     with certain number of disease susceptibility
#                     loci. 
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedDate$
# $Rev$
#
# Known limitations/bugs:
# 
#    * Currently, this script only handles constant selection pressure
#        and independent multi-locus selection model. 
#        

"""

Introduction
=============

This script uses a significantly different mechanism to control the allele 
frequency of disease susceptibility loci than simuForward.py. I will describe 
the method briefly here. More details please see 

    Peng B, Amos CI, Kimmel M (2007) Forward-Time Simulations of Human 
    Populations with Complex Diseases. PLoS Genet 3(3): e47

This program simulates the evolution of a complex common disease under the 
influence of mutation, migration, recombination and population size change. 
Starting from a small founder population, each simulation will go through
the following steps:

    1. Simulate the trajectory of allele frequency using specified disease model.
    2. Burn-in the population with mutation and recombination
    3. Introduce disease alleles and evolve the population with pre-simulated 
       allele frequency.
    4. Population structure and migration are specified along with demographic
       models. The population can be split into several equally-sized subpopulations
       and then evolve independently, or with migration. 

The result of the simulation is a large multi-generation population. To analyze 
the population, you will typically need to 

    1. Apply a penetrance function to the population and determine the affectedness
       for each individual

    2. Draw Population and/or pedigree based samples and save in popular 
       formats so that the samples can be analyzed by other software like
       genehunter.

Another script analComplexDisease gives excellent examples on how to perform these tasks.
Please follow that script to perform your own analyses.

The program is written in Python using the simuPOP modules. For more information,
please visit simuPOP website http://simupop.sourceforge.net .


Genotype structure and disease
===============================

With typical settings, each individual will have 10 chromosomes, each having
20 equal spaced microsatellite or SNP markers. A disease will be caused by 
several disease susceptibility loci (DSL) between markers. For example, a
DSL may be .3 unit to the right of marker 25 (the six marker on the second
chromosome). Since we assume that fitness is only determined by genotype, 
not affectedness status or trait value, we do not assign individual 
affectedness till the sampling stage.


Evolutionary Scenario
=====================

The evolutionary process can be divided into three stages:

Burn-in stage
------------

A founder population will be initialized with a small number of haplotypes.
All DSL will be initialized with wild type alleles ( no disease). This leads
to complete linkage disequilibrium between markers and no LD with DSL.
The population will evolve for a number of generations, subject to mutation 
and recombination. This will break down LD between markers and let (hopefully)
each marker reach a mutation-drift equilibrium.

The mutation happens only at non-DSL markers, following a symmetric
stepwise mutation model for microsattelite and a Juke-Cantor model for
SNP markers. Recombination is uniform across all markers.


No-migration stage
-------------------

The population will be split into 10 subpopulations and starts to grow,
aiming at 100,000 or more individuals at the final generation. No migration 
is allowed between subpopulations so population heterogeneity will build up.

Mixing stage
-------------

Individuals from different subpopulations will be able to migrate following
a circular step-stone model. Population heterogeneity will be reduced to
a level depending on migration rate and length of this stage.

Introduction of disease
-----------------------

The disease allele frequency is simulated before the simulation is performed.
A single disease mutant is introduce to each DSL at simulated mutant-introduction
generation. The allele frequency then evolve according to the simulated frequency
trajectory.

However, you can also specify some free DSL, who will evolve in a different manner.
Namely, they will be brought to high allele frequency in your specified time frame,
and then evolve freely. This method may fail due to extinction of disease alleles,
but it has the advantage of being able to simulate linked DSL.


Statistics Monitored
====================

A number of statistics will be measured and saved. They are:

    1. Fst before and after mixing
    2. Observed heterogeneity before and after mixing
    3. Disease allele frequency trajectory.

You can add more stat() operator if you need to preserve more information.

"""

import simuOpt
import os, sys, exceptions, types, math

#
# declare all options. getParam will use these information to get parameters
# from a tk/wxPython-based dialog, command line, config file or user input
#
options = [
    {'arg': 'h',
     'longarg': 'help',
     'useDefault': True,
     'default': False, 
     'description': 'Print this usage message.',
     'allowedTypes': [types.NoneType, type(True)],
     'jump': -1                    # if -h is specified, ignore any other parameters.
    },
    #
    # 
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
                Choose between linear and exponential''',
     'chooseOneOf': ['exponential', 'linear'],
    },
    {'longarg': 'burninGen=',
     'default': 3000,
     'label': 'Length of burn-in stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of generations of the burn in stage.',
     'validate':    simuOpt.valueGT(0)
    },
    {'longarg': 'splitGen=',
     'default': 5000,
     'label': 'When to split population',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''At which generation to split the population. 
                The population will start to grow after burnin stage but will
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
    {'longarg': 'alleleDistInSubPop=',
     'useDefault': True,
     'default': 'uneven',
     'label': 'Allele distribution in subpopulations',
     'description': '''If disease allele frequencies in each subpopulation are
                not specified, how to distribute disease alleles in subpopulations.
                If 'even' is chose, number of disease alleles will be distributed
                according to a multinomial distribution so is largely even among 
                subpopulations. If 'exponential' is chosen, the number of disease alleles
                will be proportional to the interval lengths of 0 x x x 1 while x are uniform
                [0,1]. The distribution of interval lengths, are roughly exponential 
                (conditional on overall length 1). ''',
     'allowedTypes': [types.StringType],
     'validate':    simuOpt.valueOneOf(['even', 'uneven']),
     'chooseOneOf': ['even', 'uneven']
    },
    #
    #
    {'separator': 'Disease model:'},
    {'longarg': 'curAlleleFreq=',
     'default': [0.05, 0.05, 0.05],
     'label': 'Final allele frequencies',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Current allele frequencies for each DSL.
                If a number is given, it is assumed to be the frequency
                for all DSL.''',
     'validate': simuOpt.valueListOf( simuOpt.valueBetween(0,1))
    },
    {'longarg': 'minMutAge=',
     'default': 0,
     'label': 'Minimum mutant age',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Minimum mutation age. Default to 0. Because the population
                may be split into subpopulations at splitGen, if the mutation age is too
                short, it may fall in one of the subpopulations. To void this, you can 
                set minMutAge to be endingGen - splitGen. Note that all mutants have the same
                minimum mutant age regardless final allele frequency.'''
    },
    {'longarg': 'maxMutAge=',
     'default': 0,
     'label': 'Maximum mutant age',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Maximum mutant age. Default to 0, means no max and the real
                maximum mutant age will be endingGen. However, if you do not want mutant to
                be generated in the burnin stage. You can set maxMutAge to endingGen-burnin.
                Note that all mutants have the same maximum mutant age regardless final 
                allele frequency.'''
    },
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
                 between markers.
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
    {'longarg': 'simuAge=',
     'default': 0,
     'useDefault': True,
     'allowedTypes': [types.IntType],
     'validate': simuOpt.valueGE(0),
     'description': '''Simulate the age of disease mutant given times, and print out the ages,
                The simulation is not run if a positive number is given. This option is useful
                when one want to determine minMutAge and maxMutAge for a specific demographic 
                and fitness model.'''
    },
    {'longarg': 'dryrun',
     'useDefault': True,
     'default': False,
     'description':    'Only display how simulation will perform.'
     # do not save to config, do not prompt, so this appeared to be an undocumented option.
    },
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
    {'longarg': 'simuName=',
     'default': 'simu',
     'allowedTypes': [types.StringType],
     'label': 'Simulation name',
     'description': '''Name of simulation, files saved will be 
                    name + '.log': statistics output
                    name + '.cfg': configuration
                    name + .bin/txt/xml: saved popuation''',
    },
    {'longarg': 'saveFormat=',
     'default': 'txt',
     'useDefault': True,
     'allowedTypes': [types.StringType],
     'label': 'Format to save population',
     'description': '''Format to save population, can be 
                 text, bin or xml. Note that the binary format, although
                 smallest, may not be portable between different machines.''',    
     'chooseOneOf': [ 'txt', 'bin', 'xml']
    },
    {'arg': 'v',
     'longarg': 'verbose',
     'useDefault': True,
     'default': False,
     'description': 'Verbose mode.'
    },
]

# __doc__ + parameter definitions use >500 lines, 
# more than one third of the total length of the script.

def getOptions(details=__doc__):
    ''' get options from options structure,
        if this module is imported, instead of ran directly,
        user can specify parameter in some other way.
    '''
    # get all parameters, __doc__ is used for help info
    allParam = simuOpt.getParam(options, 
        '''    This program simulates the evolution of a complex common disease, subject 
     to the impact of mutation, migration, recombination and population size change. 
     Click 'help' for more information about the evolutionary scenario.''', details, nCol=2)
    # when user click cancel ...
    if len(allParam) == 0:
        sys.exit(1)
    # -h or --help
    if allParam[0]:    
        print simuOpt.usage(options, __doc__)
        sys.exit(0)
    # change current directory, all files wil be saved here.
    try:
        os.mkdir(allParam[-3])
    except:
        pass
    # --saveConfig
    if allParam[-3] != None: # saveConfig
        simuOpt.saveConfig(options, os.path.join(allParam[-3], allParam[-3]+'.cfg'), allParam)
    # --verbose or -v (these is no beautifying of [floats]
    if allParam[-1]:                 # verbose
        simuOpt.printConfig(options, allParam)
    # return the rest of the parameters
    return allParam[1:-1]


def trajFunc(endingGen, traj):
    '''HIDDEN
    return freq at each generation from a simulated trajctories. '''
    def func(gen):
        freq = []
        for tr in traj:
            if gen < endingGen - len(tr) + 1:
                freq.append( 0 )
            else:
                freq.append( tr[ gen - (endingGen - len(tr) + 1) ] )
        return freq
    return func


def FreqTrajectoryMultiStochWithSubPop(
        curGen,
        numLoci,
        freq,
        NtFunc,
        minMutAge,
        maxMutAge,
        fitness=[],
        mode = 'uneven',
        ploidy=2,
        restartIfFail=True,
        fitnessFunc=None):
    ''' HIDDEN
    Simulate frequency trajectory with subpopulation structure,
    migration is currently ignored. The essential part of this
    script is to simulate the trajectory of each subpopulation
    independently by calling FreqTrajectoryMultiStoch with properly
    wrapped NtFunc function.

    If mode = 'even' (default)
        When freq is the same length
        of the number of loci. The allele frequency at the last
        generation will be multi-nomially distributed. If freq
        for each subpop is specified in the order of loc1-sp1, loc1-sp2, ..
        loc2-sp1, .... This freq will be used directly.

    If mode = 'uneven'.
        The number of disease alleles
        will be proportional to the interval lengths of 0 x x x 1 while x are
        uniform [0,1]. The distribution of interval lengths, are roughly
        exponential (conditional on overall length 1). '

    If mode = 'none'
        subpop will be ignored.

    This script assume a single-split model of NtFunc
    '''
    numSP = len(NtFunc(curGen))
    if maxMutAge == 0:
        maxMutAge = endGen
    TurnOnDebug(DBG_GENERAL)
    if numSP == 1 or mode == 'none':
        # given in the form of [[.1,.2]]
        if type(freq[0]) in [type((0,)), type([])]:
            freq = freq[0]
        traj = FreqTrajectoryMultiStoch(
                curGen=curGen,
                freq=freq,
                NtFunc=NtFunc,
                fitness=fitness,
                fitnessFunc=fitnessFunc,
                minMutAge=minMutAge,
                maxMutAge=maxMutAge,
                ploidy=ploidy,
                restartIfFail=True)
        #print traj
        if True in [len(t) < max(2, minMutAge) for t in traj]:
            print "Failed to generate trajectory. You may need to set a different set of parameters."
            print "len: ", [len(t) for t in traj]
            print 'curGen: ', curGen
            print 'freq: ', freq
            print 'begin size: ', NtFunc(0)
            print 'ending size: ', NtFunc(curGen)
            print 'fitness: ', fitness
            print 'minMutAge: ', minMutAge
            print 'maxMutAge: ', maxMutAge
            print 'ploidy: ', ploidy
            sys.exit(1)
        return (traj, [curGen-len(x)+1 for x in traj], trajFunc(curGen, traj))
    # other wise, do it in two stages
    # get the split generation.
    split = curGen;
    while(True):
        if len(NtFunc(split)) == 1:
            break
        split -= 1
    split += 1
    # set default for min/max mutage
    if minMutAge < curGen - split:
        minMutAge = split
    if minMutAge > maxMutAge:
        print "Minimal mutant age %d is larger then maximum age %d" % (minMutAge, maxMutAge)
        sys.exit(1)
    # now, NtFunc(split) has subpopulations
    #
    # for each subpopulation
    # layout for freqALL
    #     [  loc 0   ][loc 1 ][loc 2  ]
    # for each locus
    #     [ loc  0 ] = [at subpop 0, subpop 1, subpop 2...]
    #
    # That is to say, index is for locus i in subpop sp
    #   freqAll[sp+i*numSP],
    if type(freq[0]) in [type((0,)), type([])]:
        # freq is given as [subpop 0][subpop 1] [subpop 2]
        freqAll = [0]*(numLoci*numSP)
        for (sp, f) in enumerate(freq):
            for (loc, s) in enumerate(f):
                freqAll[sp+loc*numSP] = s
    else:
        if len(freq) == numSP*numLoci:
            freqAll = freq
        elif len(freq) == numLoci:
            freqAll = [0]*(numLoci*numSP)
            if mode == 'even':
                for i in range(numLoci):
                    wt = NtFunc(curGen)
                    ps = sum(wt)
                    # total allele number
                    totNum = int(freq[i]*ps)
                    # in subpopulations, according to population size
                    num = rng().randMultinomialVal(totNum, [x/float(ps) for x in wt])
                    for sp in range(numSP):
                        freqAll[sp+i*numSP] = num[sp]/float(wt[sp])
            elif mode == 'uneven':
                for i in range(numLoci):
                    wt = NtFunc(curGen)
                    # total allele number
                    totNum = int(freq[i]*sum(wt))
                    while(True):
                        # get [ x x x x x ] while x is uniform [0,1]
                        num = [0,1]+[rng().randUniform01() for x in range(numSP-1)]
                        num.sort()
                        for sp in range(numSP):
                            freqAll[sp+i*numSP] = (num[sp+1]-num[sp])*totNum/wt[sp]
                        if max(freqAll) < 1:
                            break;
            else:
                print "Wrong mode parameter is used: ", mode
            print "Using ", mode, "distribution of alleles at the last generation"
            print "Frequencies at the last generation: sp0-loc0, loc1, ..., sp1-loc0,..."
            for sp in range(numSP):
                print "SP ", sp, ': ',
                for i in range(numLoci):
                    print "%.3f " % freqAll[sp+i*numSP],
                print
        else:
            raise exceptions.ValueError("Wrong freq length")
    spTraj = [0]*numSP*numLoci
    for sp in range(numSP):
        print "Generting trajectory for subpopulation %d (generation %d - %d), freq=%s" % (sp, split, curGen, [freqAll[sp+x*numSP] for x in range(numLoci)])
        # FreqTraj... will probe Nt for the next geneartion.
        def spPopSize(gen):
            if gen < split:
                return [NtFunc(split-1)[0]]
            else:
                return [NtFunc(gen)[sp]]
        while True:
            t = FreqTrajectoryMultiStoch(
                curGen=curGen,
                freq=[freqAll[sp+x*numSP] for x in range(numLoci)],
                NtFunc=spPopSize,
                fitness=fitness,
                fitnessFunc=fitnessFunc,
                minMutAge=curGen-split,
                maxMutAge=curGen-split,
                ploidy=ploidy,
                restartIfFail=False)
                # failed to generate one of the trajectory
            if 0 in [len(x) for x in t]:
                print "Failed to generate trajectory. You may need to set a different set of parameters."
                sys.exit(1)
            if 0 in [x[0] for x in t]:
                print "Subpop return 0 index. restart "
            else:
                break;
        # now spTraj has LOC0: sp0,1,2, LOC1, sp0,1,2,... LOC2
        for i in range(numLoci):
            spTraj[sp+i*numSP] = t[i]
    # add all trajectories
    traj = []
    for i in range(numLoci):
        traj.append([])
        for g in range(split, curGen+1):
            totAllele = sum( [
                spTraj[sp+i*numSP][g-split] * NtFunc(g)[sp] for sp in range(numSP) ])
            traj[i].append( totAllele / sum(NtFunc(g)) )
    #
    print "Starting allele frequency (at split) ", [traj[i][0] for i in range(numLoci)]
    print "Generating combined trajsctory with range: ", minMutAge, " - ", maxMutAge
    trajBeforeSplit = FreqTrajectoryMultiStoch(
        curGen=split,
        freq=[traj[i][0] for i in range(numLoci)],
        NtFunc=NtFunc,
        fitness=fitness,
        fitnessFunc=fitnessFunc,
        minMutAge=minMutAge-len(traj[0])+1,
        maxMutAge=maxMutAge-len(traj[0])+1,
        ploidy=ploidy,
        restartIfFail=True)
    if 1 in [len(x) for x in trajBeforeSplit]:
        print "Failed to generated trajectory. (Tried more than 1000 times)"
        sys.exit(0)
    def trajFuncWithSubPop(gen):
        if gen >= split:
            return [spTraj[x][gen-split] for x in range(numLoci*numSP)]
        else:
            freq = []
            for tr in trajBeforeSplit:
                if gen < split - len(tr) + 1:
                    freq.append( 0 )
                else:
                    freq.append( tr[ gen - (split - len(tr) + 1) ] )
        return freq
    trajAll = []
    for i in range(numLoci):
        trajAll.append( [] )
        trajAll[i].extend(trajBeforeSplit[i])
        trajAll[i].extend(traj[i][1:])
    # how exactly should I return a trajectory?
    return (trajAll, [curGen-len(x)+1 for x in trajAll ], trajFuncWithSubPop)

    
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
    (burnin, split, mixing, endGen, outfile) = args
    # 
    gen = pop.gen()
    # see how long the simulation has been running
    if gen == burnin:
        print "Start population growth\t\t"
    elif gen == split:
        print "Start no-migration stage\t\t"
    elif gen == mixing:
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
    return True


# simulate function, using a single value of mutation, migration,
# recombination rate
def simuComplexDisease(numChrom, numLoci, markerType, DSLafter, DSLdistTmp, 
        initSize, endingSize, growthModel, 
        burninGen, splitGen, mixingGen, endingGen, 
        numSubPop, migrModel, migrRate, alleleDistInSubPop,
        curAlleleFreqTmp, minMutAge, maxMutAge, fitnessTmp, mlSelModelTmp, 
        mutaRate, recRate, savedGen, numOffspring, numOffMode, simuAge,
        dryrun, savePop, filename, format):
    ''' run a simulation of complex disease with given parameters. 
    '''    
    ###
    ### CHECK PARAMETERS
    ###
    if initSize < len(DSLafter):
        raise exceptions.ValueError("Initial population size is too small. (Less than number of DSL)")
    # expand .5 -> [0.5, 0.5...]
    if len(DSLdistTmp) == 1:
        DSLdist = DSLdistTmp * len(DSLafter)
    else:
        DSLdist = DSLdistTmp
    if len( DSLafter ) != len(DSLdist):
        raise exceptions.ValueError("Please specify DSL distance for each DSL.")
    numDSL = len(DSLafter)
    if burninGen > splitGen or splitGen > mixingGen or splitGen > endingGen:
        raise exceptions.ValueError("Generations should in the order of burnin, split, mixing and ending")
    # curAlleleFreq expand 0.5 -> [0.5,0.5,...]
    if len(curAlleleFreqTmp) == 1:
        curAlleleFreq = curAlleleFreqTmp * len(DSLafter)
    elif len(curAlleleFreqTmp) == len(DSLafter):
        curAlleleFreq = curAlleleFreqTmp
    else:
        raise exceptions.ValueError("min allele frequency should be either a number or a list\n" +
            " of the same length as DSLafter")
    # fitness
    if mlSelModelTmp == 'none':
        fitness = []
    elif mlSelModelTmp == 'interaction':
        if numDSL == 1:
            raise exceptions.ValueError("Interaction model can only be used with more than one DSL");
        if len(fitnessTmp) != 3**numDSL:
            raise exceptions.ValueError("Please specify 3^n fitness values for n DSL");
        fitness = fitnessTmp
    else:
        if fitnessTmp == []:    # neutral process
            fitness = [1,1,1]*numDSL
        else:
            # for a single DSL
            if len(fitnessTmp) == 3:
                fitness = fitnessTmp*numDSL
            elif len(fitnessTmp) != numDSL*3:
                raise exceptions.ValueError("Please specify fitness for each DSL")
            else:
                fitness = fitnessTmp
    # number of offspring
    if numOffMode == 'geometric' and numOffspring > 1:
        raise exceptions.ValueError("numOffspring is p for a geometric distribution when numOffMode='geometric'. It should be elss than 1.") 
    ###
    ### simulating population frequency
    ### 
    print '\n\nSimulating trajectory of allele frequencies'
    # 1. define population size function
    def popSizeFunc(gen, curSize=[]):
        if gen < burninGen:
            return [initSize]
        if growthModel == 'linear':
            inc = (endingSize-initSize)/float(endingGen-burninGen)
            if gen < splitGen:
                return [int(initSize+inc*(gen-burninGen))]
            else:
                return [int(initSize+inc*(gen-burninGen))/numSubPop]*numSubPop
        else:    # exponential
            rate =    (math.log(endingSize)-math.log(initSize))/(endingGen-burninGen)
            if gen < splitGen:
                return [int(initSize*math.exp((gen-burninGen)*rate))]
            else:
                return [int(initSize*math.exp((gen-burninGen)*rate)/numSubPop)]*numSubPop
    # 2. simulating allele frequency trajectory
    if maxMutAge == 0:
        maxMutAge = endingGen
    elif maxMutAge < minMutAge:
        raise exceptions.ValueError('maxMutAge needs to be greater than minMutAge.')
    elif maxMutAge > endingGen:
        print 'maxMutAge should be smaller than endingGen, set it to endingGen.'
        maxMutAge = endingGen
    # this is defined in simuUtil, a wrapper for
    # FreqTrajectoryMultiStoch
    if simuAge > 0:
        mutantAges = []
        for i in range(simuAge):
            (traj, introGens, trajFunc) = FreqTrajectoryMultiStochWithSubPop(
                curGen=endingGen,
                numLoci=len(DSLafter),
                freq=curAlleleFreq, 
                NtFunc=popSizeFunc,
                fitness=fitness, 
                minMutAge=minMutAge, 
                maxMutAge=maxMutAge, 
                mode=alleleDistInSubPop,
                restartIfFail=True)
            mutantAges.append([len(x) for x in traj])
        print "Simulated trajectory length:"
        for i,age in enumerate(mutantAges):
            print 'simu %3d:' % i, ', '.join(['%6d' % x for x in age])
        mean = [sum([x[col] for x in mutantAges])*1./simuAge for col in range(len(DSLafter))]
        print "Mean:    ", ', '.join(['%6.1f' % x for x in mean])
        sys.exit(0)
    (traj, introGens, trajFunc) = FreqTrajectoryMultiStochWithSubPop(
        curGen=endingGen,
        numLoci=len(DSLafter),
        freq=curAlleleFreq, 
        NtFunc=popSizeFunc,
        fitness=fitness, 
        minMutAge=minMutAge, 
        maxMutAge=maxMutAge, 
        mode=alleleDistInSubPop,
        restartIfFail=True)
    #
    # 3. save and plot the simulation scenario
    tfile = open(filename+'.traj', 'w')
    tfile.write('Simulated trajectories:\n')
    for t in traj:
        tfile.write(', '.join([str(x) for x in t ]) + '\n')
    tfile.write('Observed allele frequency (generation, population size, allele frequencies)\n')
    tfile.close()
    print '\n\nTrajectories simulated successfully. Their lengths are %s.\n Check %s.eps for details. \n\n' \
        % ( ', '.join([ str(len(x)) for x in traj]), filename)
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
    #
    # I can calculate the real indices of DSL (after insertion of DSL at DSLafter),
    # but it is easier, and less error prone to use a temporary population to do it.
    #
    # create a population, give each locus a name
    # 0, 1, 2, 3, ...
    tmp = population(1, loci=[numLoci]*numChrom, lociNames=[str(x) for x in range(numLoci*numChrom)])
    # insert DSL, given them names DSL0, DSL1, ...
    tmp.insertAfterLoci(idx=DSLafter, pos=[tmp.locusPos(DSLafter[x]) + DSLdist[x] for x in range(len(DSLafter))], 
        names = ['DSL%d' % x for x in range(len(DSLafter))])
    # now, obtain the indices of these loci using loci names
    # I need these information because the following operators need to know which marker is DSL...
    DSL = tmp.lociByNames(['DSL%d' % x for x in range(len(DSLafter))])
    loci = tmp.numLoci()
    lociPos = tmp.lociPos()
    nonDSL = tmp.lociByNames([str(x) for x in range(numLoci*numChrom)])
    ###
    ### initialization 
    ###
    if maxAle > 1:    # Not SNP
        preOperators = [
            # initialize all loci with 5 haplotypes
            initByValue(value=[[x]*sum(loci) for x in range(48, 53)],
                proportions=[.2]*5), 
            # and then init DSL with all wild type alleles
            initByValue([0]*len(DSL), atLoci=DSL)
        ]
    else: # SNP
        preOperators = [
            # initialize all loci with two haplotypes (0001,111)
            initByValue(value=[[x]*sum(loci) for x in [0,1] ],
                proportions=[.5]*2), 
            # and then init DSL with all wild type alleles
            initByValue([0]*len(DSL), atLoci=DSL)
        ]            
    ###
    ### mutation, start from gen 0,
    ###
    if maxAle > 1:    # Not SNP
        # symmetric mutation model for microsatellite
        mutator = smmMutator(rate=mutaRate, maxAllele=maxAle, atLoci=nonDSL)
    else:
        # k-allele model for mutation of SNP
        mutator = kamMutator(rate=mutaRate, maxAllele=1, atLoci=nonDSL)
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
            al.extend( [x+start for x in range(loci(ch) )] )
            start += loci[ch]
        rec = recombinator(rate=recRate, afterLoci=al)
    ###
    ### output progress
    ###
    operators = [
        mutator, 
        rec, 
        stat(alleleFreq=DSL, popSize=True, LD=[DSL[0], DSL[0]+1], step=1),
        # output to screen
        pyEval( expr=r'"%d(%d): "%(gen, popSize) + " ".join(["%.5f"%(1-alleleFreq[x][0]) for x in DSL])+" %f\n"%LD_prime[DSL[0]][DSL[0]+1]',
            step=100),
        # output to file (append)
        pyEval( expr=r'"%d %d " %(gen, popSize) + " ".join(["%.5f"%(1-alleleFreq[x][0]) for x in DSL])+"\n"',
            output='>>>'+filename+'.traj')
    ]
    ###
    ### introduction of disease mutants
    ###
    ###                                 endingGen
    ###            0 1 ...... i_T
    for i in range(numDSL):
        operators.append( 
            pointMutator(atLoci=[DSL[i]], toAllele=1, inds=[i],
                at = [introGens[i]], stage=PreMating ) ) 
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
    if mlSelModelTmp in ['additive', 'multiplicative']:
        mlSelModel = {'additive':SEL_Additive, 
            'multiplicative':SEL_Multiplicative}[mlSelModelTmp]
        operators.append( mlSelector(
            # with five multiple-allele selector as parameter
            [ maSelector(locus=DSL[x], wildtype=[0], 
                fitness=[fitness[3*x],fitness[3*x+1],fitness[3*x+2]]) for x in range(len(DSL)) ],
            mode=mlSelModel, begin=burninGen),
        )
    elif mlSelModelTmp == 'interaction':
        # multi-allele selector can handle multiple DSL case
        operators.append( maSelector(loci=DSL, fitness=fitness, wildtype=[0],
            begin=burninGen) )
    ###
    ### migration
    ###
    if numSubPop > 1 and migrModel == 'island' and migrRate > 0:
        operators.append( migrator(migrIslandRates(migrRate, numSubPop),
            mode=MigrByProbability, begin=mixingGen) )
    elif numSubPop > 1 and migrModel == 'stepping stone' and migrRate > 0:
        operators.append( migrator(migrSteppingStoneRates(migrRate, numSubPop, 
            circular=True),    mode=MigrByProbability, begin=mixingGen) )
    ###
    ### output statistics, track performance
    ###
    operators.extend([
        pyOperator(func=outputStatistics, 
            param = (burninGen, splitGen, mixingGen, endingGen, filename+'.log'),
            at = [burninGen, splitGen, mixingGen, endingGen ] ), 
        # show elapsed time
        ticToc(at=[burninGen, splitGen, mixingGen, endingGen]),
        ]
    )
    if len(savePop) > 0:
        # save population at given generations
        operators.append(savePopulation(outputExpr='os.path.join(simuName, "%s_%d.%s" % (simuName, gen, format))', 
            at=savePop))
    ### 
    ### prepare mating scheme
    ###
    ### use ancestralDepth, 
    def saveAncestors(gen):
        if gen >= endingGen - savedGen:
            return savedGen
        else:
            return 1
    # create a population, note that I add needed information 
    # fields father_idx, mother_idx later on, with the hope 
    # that simulation can run a bit faster without them. 
    pop = population(subPop=popSizeFunc(0), ploidy=2,
        loci = loci, maxAllele = maxAle, lociPos = lociPos,
        infoFields = ['fitness', 'father_idx', 'mother_idx'])
    # save DSL info, some operators will use it.
    pop.dvars().DSL = DSL
    pop.dvars().numLoci = numLoci
    pop.dvars().curAlleleFreq = curAlleleFreq
    # clear log file if it exists
    try:
        os.remove(filename+'.log')
    except:
        pass
    #
    # start simulation.
    simu = simulator( pop, 
        controlledRandomMating(
            newSubPopSizeFunc=popSizeFunc,            # demographic model
            loci=DSL,                                                     # which loci to control
            alleles=[1]*numDSL,                                 # which allele to control
            freqFunc=trajFunc                                     # frequency control function
        ),
        rep=1)
    # evolve! If --dryrun is set, only show info
    simu.evolve( preOps = preOperators, ops = operators, 
        end=endingGen-savedGen, dryrun=dryrun )
    if dryrun:
        raise exceptions.SystemError("Stop since in dryrun mode.")
    # prepare for the last several generations
    # change mating scheme to random mating.
    #
    # save saveGen-1 ancestral generations
    simu.population(0).setAncestralDepth(savedGen-1)
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
    simu.evolve(ops=operators, end=endingGen)
    # succeed save information
    pop = simu.population(0)
    # we want to save info on how this population is generated.
    # This is not required but is a good practise
    pop.dvars().DSLafter = DSLafter
    pop.dvars().DSLdist = DSLdist
    pop.dvars().initSize = initSize
    pop.dvars().endingSize = endingSize
    pop.dvars().burninGen = burninGen
    pop.dvars().splitGen = splitGen
    pop.dvars().mixingGen = mixingGen
    pop.dvars().endingGen = endingGen
    pop.dvars().numSubPop = numSubPop
    pop.dvars().mutaRate = mutaRate
    pop.dvars().mutaModel = "symmetric stepwise"
    pop.dvars().migrRate = migrRate
    pop.dvars().alleleDistInSubPop = alleleDistInSubPop
    pop.dvars().migrModel = "circular stepping stone"
    pop.dvars().recRate = recRate
    pop.dvars().numOffspring = numOffspring
    pop.dvars().numOffMode = numOffMode
    print "Saving population to " + filename + '.' + format + '\n'
    #TurnOnDebug(DBG_UTILITY)
    simu.population(0).savePopulation(filename+'.'+format)
    return True

if __name__ == '__main__':
    ############## GET OPTIONS ####################################
    allParam = getOptions()
    # unpack options
    (numChrom, numLoci, markerType, DSLafter, DSLdist, 
        initSize, endingSize, growthModel, 
        burninGen, splitGen, mixingGen, endingGen, 
        numSubPop, migrModel, migrRate, alleleDistInSubPop,
        curAlleleFreq, minMutAge, maxMutAge, fitness, selMultiLocusModel,
        mutaRate, recRate, savedGen, numOffspring, numOffMode, simuAge,
        dryrun, savePop, simuName, format) = allParam
    #
    if markerType == 'SNP':
        simuOpt.setOptions(alleleType='binary')
    else:
        simuOpt.setOptions(alleleType='short')
        
    #simuOpt.setOptions(quiet=True)

    # load simuPOP libraries
    from simuPOP import *
    from simuUtil import *
    #
    # check minimal requirement of simuPOP version
    if simuRev() < 383:
        raise exceptions.SystemError('''This scripts requires simuPOP revision >=425. 
            Please upgrade your simuPOP distribution.''' )
    #
    ################## RUN THE SIMULATION ###############
    simuComplexDisease(numChrom, numLoci, markerType, DSLafter, DSLdist, 
        initSize, endingSize, growthModel, 
        burninGen, splitGen, mixingGen, endingGen, 
        numSubPop, migrModel, migrRate, alleleDistInSubPop, 
        curAlleleFreq, minMutAge, maxMutAge, fitness, selMultiLocusModel, 
        mutaRate, recRate, savedGen, numOffspring, numOffMode, simuAge,
        dryrun, savePop, os.path.join(simuName, simuName), format)
    
    print "Done!"
