#!/usr/bin/env python
#
# Purpose:  generate dataset for common complex disease 
#           with certain number of disease susceptibility
#           loci.
#
# Bo Peng (bpeng@rice.edu)
# April, 2005
#
# Last revision date:
#   Sep, 2005
#
# Known bugs:
#   None
# 

"""

Introduction
=============

This program simulates the evolution of a complex common disease under the 
influence of mutation, migration, recombination and population size change. 
Starting from a small founder population, each simulation will go through
the following four stages:

  1. Burn-in the population with mutation and recombination
  2. Introduce disease mutants and make them common with advantagous selection
  3. Split and grow the population without migration, disease alleles
     may be under selection pressure.
  4. Mix subpopulations at given migration level

The result of the simulation is a set of large populations. For each population,

  1. Single or multi-locus penetrance functions are applied to determine 
     affectedness of each individual. Note that the penetrance model would
     better be compatible to the fitness model. You would not want to assign
     affectedness to individuals according to disease susceptibility locus 
     (DSL) one while selection was working on DSL two.
  2. Population and/or pedigree based samples are drawn and saved in popular
     formats.
  3. If geneHunter is available, it will be used to analyze affected sibpair 
     samples using TDT method.

The program is written in Python using simuPOP modules. For more information,
please visit simuPOP website http://simupop.sourceforge.net .


Genotype structure and disease
===============================

With typical settings, each individual will have 20 chromosomes, each having
20 equal spaced microsatellite or SNP markers. A disease will be caused by 
several disease susceptibility loci (DSL) between markers. For example, a
DSL may be .3 unit to the right of marker 25 (the six marker on the second
chromosome using a zero-based index). Since we assume that fitness is 
only determined by genotype, not affectedness status or trait value, we 
do not assign individual affectedness till the sampling stage.


Evolutionary Scenario
=====================

The evolutionary process can be divided into four stages:

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


Disease-introduction stage
--------------------------

During this relatively short stage, mutants will be introduced to 
individuals (one DSL per individual). The mutants will be given strong 
selective advantage so they will reach a high allele frequency quickly. 
A mutant will be re-introduced if it gets lost due to genetic drift.

At the end of this stage, allele frequencies at DSL will be checked.
If they are too low or too high, the simulation will be abandoned and
restart again. Since new mutants get lost quickly, it is unlikely that
allele frequencies at DSL will meet our requirement at the first time.
As a matter of fact, most of the simultion time will be spent on this
trial and restart process.


No-migration stage
-------------------

The population will be split into 10 subpopulations and starts to grow,
aiming at 100,000 or more individuals at the end. No migration is allowed
between subpopulations so population heterogeneity will build up.
Negative selection may decrease the allele frequency of DSL.


Mixing stage
-------------

Individuals from different subpopulations will be able to migrate following
a circular step-stone model. Population heterogeneity will be reduced to
a level depending on migration rate and length of this stage.


Penetrance
==========

Since we assume that fitness only depends on genotype, not affectedness status,
we do not care who are affected during evolution. This has to change at the 
last generation where different sampling schemes will be applied to the 
population. Several customizable penetrance schemes will be used. As a matter
of fact, if there is no selection against any DS allele, we can use any
of the penetrance functions:

  1. recessive single-locus and heterogeneity multi-locus model: 
    At DSL i , penetrance will be computed as
      0, 0, p_i
    for genotype
      aa, aA, AA  (a is wild type)
    where pi are penetrance factors (that can vary between DSL).

    The overall penetrance is 
      1 - Prod(1-P_i)
    where Pi is the penetrance at locus i.
  
  2. additive single-locus and heterogeneity multi-locus model: 
    For each DSL, the penetrance is
      0, p/2, p
    for genotype
      aa, aA, AA
    where the overall penetrance takes the form of
      1 - Prod( 1- P_i)
    This is the heterogeneity model proposed by Neil Risch (1990).
     
  3,4,5,6 are customized penetrances. You will have to modify
    this script to apply them.
  
Statistics Monitored
====================

A number of statistics will be measured and saved. They are:

  1. Fst before and after mixing
  2. Observed heterogeneity before and after mixing
  3. LD (D') between markers and between a DSL and surrounding markers
     at the end of each stage.
  4. Disease allele frequency at the end of each stage.


Samples and Output
==================

Different kinds of samples will be draw from the final large population.

  1. population based case control sample: 
     Regardless of family structure, N cases and N controls will
     be drawn randomly.
  
  2. affected and unaffected sibpairs:
     N/4 affected and N/4 unaffected (sibling) families (two siblings and
     two parents) will be drawn. (Sample size is N cases and N controls
     when counting individuals.)

The datasets will be saved in native simuPOP format and in Linkage format.
DSL markers will be removed so there will be no marker that is directly 
linked to the disease.

All files are put under a specified folder. They are organized by population
and penetrance methods. A html file summary.htm will be automatically 
generated with links to all statistics, datasets etc.


Gene mapping
============

If the location of genehunter is specified. It will be applied to all affected
sibpair samples. This is the raw one-locus TDT method with no correct on things
like multiple testing.

"""

import simuOpt
import os, sys, types, exceptions, os.path, re, math, time

#
# declare all options. getParam will use these information to get parameters
# from a tk/wxPython-based dialog, command line, config file or user input
#
# detailed information about these fields is given in the simuPOP reference
# manual.
options = [
  {'arg': 'h',
   'longarg': 'help',
   'default': False, 
   'description': 'Print this usage message.',
   'allowedTypes': [types.NoneType, type(True)],
   'jump': -1          # if -h is specified, ignore any other parameters.
  },
  {'longarg': 'numChrom=',
   'default': 20,
   'configName': 'Number Of Chromosomes',
   'prompt': 'Number of chromosomes (20):  ',
   'description': 'Number of chromosomes.',
   'allowedTypes': [types.IntType],
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'numLoci=',
   'default': 20,
   'configName': 'Number of loci on each chrom',
   'prompt': 'Number of loci on each chromosome (20):  ',
   'description': 'Number of loci on each chromosome.',
   'allowedTypes': [types.IntType],
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'markerType=',
   'default': 'microsatellite',
   'configName': 'Marker type',
   'prompt': 'Marker type (microsatellite):  ',
   'description': '''Type of markers. Can be microsatellite or SNP.
        Microsatellite markers will be mutated using a symmetric
        stepwise model while SNP markers will be mutated using a
        2-state model. Mutation rate should be much smaller for
        SNP markers than that of microsatellite markers.''',
   'allowedTypes': [types.StringType],
   'validate':  simuOpt.valueOneOf(['microsatellite', 'SNP']),
   'chooseOneOf': ['microsatellite', 'SNP']
  },
  {'longarg': 'DSL=',
   'default': [45, 125, 310],
   'configName': 'DSL After Marker',
   'prompt': 'Enter a list of 0-indexed disease loci ([45,125,310]):  ',
   'description': '''A list of loci *after* a marker. For example, 
        35 means a disease locus after the 16th marker on chromosome 2,
        (if numChrom=numLoci=20). The number of DSL is important since
        it determines the complexity of the disease.
        A single number is allowed and implies a simple disease with
        one disease locus.''',
   'allowedTypes': [types.TupleType, types.ListType],
   'validate':  simuOpt.valueListOf( simuOpt.valueGT(0) )
  },
  {'longarg': 'DSLLoc=',
   'default': [.5, .5, .5],
   'configName': 'DSL Location between markers',
   'prompt': 'Enter the position of each DSL between two markers ([.01,.5,.75]):  ',
   'description': '''A list of loci location between two markers.
        Since all disease loci will be *between* equal spaced markers,
        the location should be between 0 and 1. A single value is acceptable
        as the location of all DSL.''',
   'allowedTypes': [types.TupleType, types.ListType],
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },  
  {'longarg': 'initSize=',
   'default': 1000,
   'configName': 'initial Population Size',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Initial Population size (1000):  ',
   'description': '''Initial population size. This size will be maintained
        till the end of burnin stage''',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'meanInitAllele=',
   'default': 50,
   'configName': 'mean initial alleles for the markers',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Mean initial alleles for the markers, microsatellite only (50):  ',
   'description': '''Initial haplotype for the markers. This option ignored for SNP markers
        since 1111..., 2222.... will be used as initial haplotypes. For 
        microsatellite markers, this number will be the mean of the alleles used.
        In this simulation, markers will be initialized with one of five haplotypes
				given by m-2,m-1,m,m+1,m+2 while m is the mean initial allele. ''',
   'validate':  simuOpt.valueGT(0)
  }, 
  {'longarg': 'burnin=',
   'default': 1000,
   'configName': 'Length of burnin stage',
   'allowedTypes': [types.IntType],
   'prompt': 'Length of burn in stage (1000):  ',
   'description': 'Number of generations of the burn in stage.',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'introGen=',
   'default': 12,
   'configName': 'Length of Disease-Intro stage',
   'allowedTypes': [types.IntType],
   'prompt': 'Length of disease-introduction stage (12):  ',
   'description': '''Number of generations to introduce the disease
        Since the disease will be under positive selection during this
        stage, large introGen will lead to over-common diseases.''',
   'validate':  simuOpt.valueGT(0),
  },
  {'longarg': 'introSel=',
   'default': -1,
   'configName': 'Selection coef of mutants',
   'allowedTypes': [types.FloatType, types.IntType],
   'prompt': 'Selection coef of mutants during disease-introduction stage (-1):  ',
   'description': '''Selection coef of mutants during disease-introduction stage.
        Individuals with these mutants will have fitness values 1,1-s/2,1-s for genotype
        AA,Aa,aa (a is mutatnt). The deault value is -0.5 meaning strong positive 
        selection on these mutants. You can set this value to zero to disable such
        artificial boost of allele frequency.''',
  },
  {'longarg': 'minAlleleFreq=',
   'default': 0.01,
   'configName': 'minimal Allele Frequency',
   'allowedTypes': [types.FloatType],
   'prompt': 'Minimum allele frequency required for each DSL (0.01):  ',
   'description': '''Mininal allele frequencies required for each DSL,
        The simulation will restart if allele frequency is lower than
        this number.''',
   'validate':  simuOpt.valueBetween(0,1)
  },
  {'longarg': 'maxAlleleFreq=',
   'default': 0.2,
   'configName': 'maximum Allele Frequency',
   'allowedTypes': [types.FloatType],
   'prompt': 'Maximum allele frequency required for each DSL (0.2):  ',
   'description': '''Maximum allele frequencies required for each DSL,
        The simulation will restart if allele frequency is greater than
        this number.''',
   'validate':  simuOpt.valueBetween(0,1)
  },
  {'longarg': 'fitness=',
   'default': [1,1,1],
   'configName': 'Fitness of genotype AA,Aa,aa ',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': ':  Fitness of genotype AA,Aa,aa ([1,1,1]): ',
   'description': '''Fitness of genotype AA,Aa,aa for every DSL after disease introduction
        stage. This can be an array (default) of fitness values for genotype AA,Aa,aa ( A
        is wild type) or an array of such array for each DSL.''',
   'validate':  simuOpt.valueGE(0.),
  },
  {'longarg': 'selMultiLocusModel=',
  'default': 'none',
  'configName': 'Multi-locus selection model for all DSL',
  'prompt': 'selection model for the common disease (additive): ',
  'description': '''Model of overall fitness value given fitness values for each DSL.
        fitness values are Prod(f_i) for multiplicative model and
        1-Sum(1-f_i) for additive model. ''',
  'allowedTypes': [types.StringType],
  'chooseOneOf': [ 'additive', 'multiplicative', 'none']
  }, 
  {'longarg': 'numSubPop=',
   'default': 10,
   'configName': 'Number of split subpops',
   'allowedTypes': [types.IntType],
   'prompt': 'Number of subpopulations to split (10):  ',
   'description': 'Number of subpopulations to be split into after burnin stage.',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'finalSize=',
   'default': 200000,
   'configName': 'Final Population Size',
   'prompt': 'Final population size (sum of all subpopulations) (200000):  ',
   'allowedTypes': [types.IntType, types.LongType],
   'description': 'Ending population size (after expansion.',
   'validate':  simuOpt.valueGT(0)
  }, 
  {'longarg': 'noMigrGen=',
   'default': 150,
   'configName': 'Length of no-Migration stage',
   'prompt': 'Length of no-migration stage (build up of population structure (150):  ',
   'allowedTypes': [types.IntType, types.LongType],
   'description': '''Number of generations when migration is zero. This stage
        is used to build up population structure.''',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'mixingGen=',
   'default': 50,
   'configName': 'Length of mixing stage',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Length of mixing stage (population admixing) (50):  ',
   'description': '''Number of generations when migration is present. This stage
        will mix individuals from subpopulations using an circular stepstone
        migration model.''',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'growth=',
   'default': 'exponential',
   'configName': 'Population Growth Model',
   'prompt': 'Population growth style, linear or exponential. (exponential):  ',
   'description': '''How population is grown from initSize to finalSize.
        Choose between linear and exponential''',
   'chooseOneOf': ['exponential', 'linear'],
  },
  {'longarg': 'migrModel=',
   'default': 'stepstone',
   'configName': 'Migration Model',
   'prompt': 'Mutation model. (stepstone):  ',
   'allowedTypes': [types.StringType],
   'description': '''Migration model. Choose between stepstone (circular),
        island and none. ''',
   'validate':  simuOpt.valueOneOf(['island', 'stepstone', 'none']),
   'chooseOneOf': ['stepstone', 'island', 'none']
  }, 
  {'longarg': 'migrRate=',
   'default': [0,0.001,0.01,0.1],
   'configName': 'Migration Rates',
   'prompt': 'Migration rate(s) during mixing stage. A separate dataset\n' +
        'will be genrated for each of the given migration rate. ([0, 0.01, 0.05, 0.1]):  ',
   'description': '''Migration rate during mixing stage. Can be a number or an array.
        A circular stepstone migration model will be used. 
        Separate datasets will be generated for each value of migration rate.''',
   'allowedTypes': [types.ListType, types.TupleType],
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },
  {'longarg': 'mutaRate=',
   'default': [0.001],
   'configName': 'Mutation Rates',
   'prompt': 'Mutation rate at non-DSL markers. (0.001):  ',
   'allowedTypes': [types.ListType, types.TupleType],
   'description': '''Microsatellite markers are mutated using  
        symmetric stepwise mutation wile SNP markers are mutaed
        using a 2-allele model (kam) DSL are not mutated unless in disease
        introduction stage.''',
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },
  {'longarg': 'recRate=',
   'default': [0.1],
   'configName': 'Recombination Rates',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Recombination rate between adjacent markers. Separate datasets\n'+
        'will be generated for each rec values if a list is given. (0.1):  ',
   'description': '''Recombination rate between adjacent markers. if a list is
        given, separate datasets will be given for each value.''',
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },
  {'longarg': 'rep=',
   'default': 2,
   'configName': 'Replicates for each settings',
   'allowedTypes': [types.IntType],
   'prompt': 'Number of populations for each settings (2):  ',
   'description': '''For each migration, mutation .. settings, this is the
        number of replicates.''',
   'validate':  simuOpt.valueGT(0)   
  },
	# temporarily remove randTent format.
  {'longarg': 'saveFormat=',
   'default': ['simuPOP','Linkage'],
   'configName': 'saveInFormat',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': "Save datasets in format (['simuPOP','Linkage','randTent']):  ",
   'description': '''Save generated datasets in specified formats.
        Choosen from simuPOP, Linkage, randTent. ''',
   'validate':  simuOpt.valueListOf( simuOpt.valueOneOf([ 'simuPOP', 'Linkage'])),
   'chooseFrom': [ 'simuPOP', 'Linkage']
  },
  {'longarg': 'peneFunc=',
   'default': ['recessive','additive'],
   'configName': 'penetrance functions',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Penetrance to be used: (recessive, additive):  ',
   'description': '''\
			  Penetrance functions to be applied to the final
        population. Two penetrance fucntions are provided, namely recessive or additive
				single-locus model with heterogeneity multi-locus model. You can define another
				customized penetrance functions by modifying this script. ''',
   'allowedTypes': [types.ListType, types.TupleType],
   'validate':  simuOpt.valueListOf( simuOpt.valueOneOf(['recessive', 'additive', 'custom'])),
   'chooseFrom': [ 'recessive', 'additive', 'custom']
  },
  {'longarg': 'penePara=',
   'default': [0.5],
   'configName': 'penetrance parameter',
   'prompt': 'Penetrance parameter used by penetrance functions. \n' + 
         "Can be an array (for each DSL). (0.5) ",
   'description': '''Penetrance parameter for all DSL. An array of parameter 
        can be given to each DSL. The meaning of this parameter differ by penetrance model.
				For a recessive model, the penetrance is 0,0,p for genotype AA,Aa,aa (a is disease
				allele) respectively. For an additive model, the penetrance is 0,p/2,p respectively.''',
   'allowedTypes': [types.ListType, types.TupleType],
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },
  {'longarg': 'sampleSize=',
   'default': 800,
   'configName':  'Sample Size',
   'allowedTypes':  [types.IntType, types.LongType],
   'prompt':  'Size of the samples (800):  ',
   'description':  '''Size of the samples, that will mean N/4 affected sibpair families (of size 4),
        N/2 cases and controls etc. ''',
   'validate':  simuOpt.valueGT(1)
  },
  {'longarg': 'outputDir=',
   'default': '.',
   'allowedTypes': [types.StringType],
   'configName': 'output Directory',
   'prompt': 'Save datasets into directory (.):  ',
   'description': 'Directory into which the datasets will be saved. ',
   'validate':  simuOpt.valueValidDir()
  },
  {'longarg': 'overwrite=',
   'default': True,
   # for compatibility. python2.2 use IntType, python 2.3 use BooleanType
   'allowedTypes': [type(True)],  
   'description': 'Whether or not overwrite existing data files. ',
   'validate':  simuOpt.valueOneOf([True, False])
  },
  {'longarg': 'geneHunter=',
   'default': '',
   'allowedTypes': [types.StringType],
   'configName': 'Location of gene hunter',
   'prompt': 'Provide location of gene hunter():  ',
   'description': '''Location of gene hunter executable. If provided,
        the TDT method of genehunter will be applied to affected sibpair samples.'''
  },
  {'longarg': 'dryrun',
   'default': False,
   'allowedTypes': [types.IntType],
   'validate':  simuOpt.valueOneOf([True, False]),
   'description':  'Only display how simulation will perform.'
   # do not save to config, do not prompt, so this appeared to be an undocumented option.
  },
  {'longarg': 'saveConfig=',
   'default': sys.argv[0].split('.')[0]+'.cfg',
   'allowedTypes': [types.StringType, types.NoneType],
   'configName': 'saveConfig',
   'prompt': 'Save current configuration to file (' + sys.argv[0].split('.')[0] + '.cfg):  ',
   'description': 'Save current paremeter set to specified file.'
  },
  {'arg': 'v',
   'longarg': 'verbose',
   'default': False,
   'allowedTypes': [types.NoneType, types.IntType],
   'description': 'Verbose mode.'
  },
]


def getOptions(details=__doc__):
  ''' get options from options structure,
    if this module is imported, instead of ran directly,
    user can specify parameter in some other way.
  '''
  #
  # get all parameters, __doc__ is used for help info
  allParam = simuOpt.getParam(options, 
    '''This program simulates the evolution of a complex common disease 
  under the effect of mutation, migration, recombination and population
  size change. Click 'help' for more info.''', details, nCol=2)
  #
  # when user click cancel ...
  if len(allParam) == 0:
    sys.exit(1)
  #  
  # -h or --help
  if allParam[0]:  
    print simuOpt.usage(options, __doc__)
    sys.exit(0)
  #
  # --saveConfig
  if allParam[-2] != None: # saveConfig
    simuOpt.saveConfig(options, allParam[-2], allParam)
  #
  # --verbose or -v (these is no beautifying of [floats]
  if allParam[-1]:         # verbose
    for p in range(0, len(options)):
      if options[p].has_key('configName'):
        if type(allParam[p]) == types.StringType:
          print options[p]['configName'], '\t"'+str(allParam[p])+'"'
        else:
          print options[p]['configName'], '\t', str(allParam[p])
  # return the rest of the parameters
  return allParam[1:-1]

#
# simulate function, using a single value of mutation, migration,
# recombination rate
# 
def simuComplexDisease( numChrom, numLoci, markerType, DSLafter, DSLdist,
    initSize, meanInitAllele, burnin, introGen, introSel, minAlleleFreq,
    maxAlleleFreq, fitness, mlSelModel, numSubPop, finalSize, noMigrGen,
    mixingGen, popSizeFunc, migrModel, mu, mi, rec, dryrun, logFile):
  ''' run a simulation of complex disease with
    given parameters. 
  '''
  #
  if markerType == 'microsatellite':
    maxAle = 99                   # max allele
  else:
    maxAle = 2                    # SNP 
  #
  numDSL = len(DSLafter)
  if len(DSLdist) != numDSL:
    raise exceptions.ValueError("--DSL and --DSLloc has different length.")
  #
  # We have numChrom*numLoci markers + DSL
  # so this loci and lociDist is the *real* loci and lociDist
  loci = [0]*numChrom
  lociDist = []
  i = 0  # current locus
  j = 0  # current DSL
  for ch in range(0, numChrom):
    lociDist.append([])
    for loc in range(0,numLoci):
      if j < numDSL and i == DSLafter[j]:
        loci[ch] += 2
        lociDist[ch].append(loc+1)
        lociDist[ch].append(loc+1 + DSLdist[j])
        j += 1
      else:
        loci[ch] += 1
        lociDist[ch].append(loc+1)
      i += 1
  #
  # adjust DSL indices since DSL[0] will be the one after marker that marker
  # DSL will be used by expressions and has to be global due to python's
  # 2-layer namespace restriction.
  global DSL
  DSL = list(DSLafter)
  for i in range(0,numDSL):
    DSL[i] += i+1
  # 
  # non DSL loci, for convenience
  nonDSL = range(0, sum(loci))
  for loc in DSL:
    nonDSL.remove(loc)
  #
  # We need to have something to hold all genotype info
  # This info can be obtained from any population or simulator
  # but we do not have one yet. 
  #
  # this is a walkaround --- using a dummy population
  gt = population(size=1, ploidy=2, loci = loci, maxAllele=maxAle, lociDist=lociDist)
  # 
  # In this analysis, we also need to know
  # 1. a DSL that is closest to the center of a chromosome
  toCtrDist = [abs(gt.locusDist(x)-numLoci/2) for x in DSL]
  ctrDSL = DSL[ [x==min(toCtrDist) for x in toCtrDist].index(True)]
  ctrChrom = gt.chromLocusPair(ctrDSL)[0]
  # 2. a chromosome without DSL
  noDSLChrom = [x==numLoci for x in loci].index(True)
  #
  # The LD pairs are used elsewhere so I set it here
  i = gt.chromBegin( ctrChrom)
  ctrDSLLD = [ [x, ctrDSL] for x in range(i, ctrDSL)] + \
    [ [ctrDSL, x] for x in range(ctrDSL+1, i+numLoci+1)]
  i = gt.chromBegin( noDSLChrom)
  noDSLLD = [ [i+x, i+numLoci/2] for x in range(numLoci/2)] + \
    [ [i + numLoci/2, i+x] for x in range(numLoci/2+1, numLoci)]
  #
  # generation related
  split  = burnin + introGen
  mixing = split  + noMigrGen
  endGen = mixing + mixingGen  
  #
  # let us go through all stages and meet all the requirements one by one
  #
  # initialization and mutation
  if maxAle > 2:  # Not SNP
    preOperators = [
      # initialize all loci with 5 haplotypes
      initByValue(value=[[x]*gt.totNumLoci() for x in range(meanInitAllele-2,meanInitAllele+2)],
        proportions=[.1]*10), 
      # and then init DSL with all wild type alleles
      initByValue([1]*len(DSL), atLoci=DSL)
    ]
    # symmetric mutation model
    mutator = smmMutator(rate=mu, maxAllele=maxAle,
      atLoci=nonDSL)
  else: # SNP
    preOperators = [
      # initialize all loci with two haplotypes
      initByValue(value=[[x]*gt.totNumLoci() for x in range(1,3)],
        proportions=[.5]*2), 
      # and then init DSL with all wild type alleles
      initByValue([1]*len(DSL), atLoci=DSL)
    ]      
    mutator = kamMutator(rate=mu, maxAllele=gt.maxAllele(),
      atLoci=nonDSL)
  #
  # burn in stage
  #
  # mutation and recombination will always be in effective
  # so we do not have to specify them later
  operators = [
    # mutator, differ by marker type
    mutator, 
    # recombination, use intensity since loci (include DSL)
    # are not equally spaced 
    recombinator(intensity=rec),
  ]
  #
  # disease introduction stage
  operators.append(
    # need to measure allele frequency at DSL
    stat(alleleFreq=DSL, begin = burnin)
    )
  # need five conditional point mutator to
  # introduce disease. It does not matter that we introduce
  # mutant specifically to individuel i since individuals
  # are unordered
  #
  for i in range(numDSL):
    # this operator literally says: if there is no DS allele,
    # introduce one. Note that operator stat has to be called
    # before this one.
    operators.append( 
      ifElse("alleleFreq[%d][1]==1." % DSL[i],
        pointMutator(atLoci=[DSL[i]], toAllele=2, inds=[i]),
      begin=burnin, end=split) 
      )
  #
  # optionally, the mutants will be given some selective pressure
  # if introSel < 0, mutants will have selective advantage and will
  # reach high allele frequency quickly.
  operators.extend([
    mlSelector(
      # with five multiple-allele selector as parameter
      [ maSelector(locus=x, wildtype=[1], fitness=[1, 1-introSel/2, 1-introSel]) for x in DSL ],
      mode=SEL_Multiplicative, begin=burnin, end=split),
    #  
    # and the simulation will stop if the disease allele frequencies
    # are not within range at the end of this stage
    # this process will continue to the end
    terminateIf("True in [(1.-alleleFreq[x][1] < minAlleleFreq or 1.-alleleFreq[x][1] > maxAlleleFreq) for x in DSL]",
      begin=split), 
    ])       
  #
  # no migration stage
  #
  if numSubPop > 1:
    operators.append( 
      # split population after burnin, to each sized subpopulations
      splitSubPop(0, proportions=[1./numSubPop]*numSubPop, at=[split]),
      )
  #
  # start selection
  if mlSelModel != SEL_None:
    operators.append(
      mlSelector(
      # with five multiple-allele selector as parameter
      [ maSelector(locus=DSL[x], wildtype=[1], fitness=fitness[x]) for x in range(len(DSL)) ],
      mode=mlSelModel, begin=split),
     )
  #
  # mixing stage
  #
  # a migrator, stepstone or island
  if numSubPop > 1 and migrModel == 'island':
    operators.append( migrator(migrIslandRates(mi, numSubPop),
      mode=MigrByProbability, begin=mixing) )
  elif numSubPop > 1 and migrModel == 'stepstone':
    operators.append( migrator(migrStepstoneRates(mi, numSubPop, circular=True),
      mode=MigrByProbability, begin=mixing) )
  #
  # prepare for sampling:
  # use last_two numOffspringFunc, simuPOP will produce 2 offspring at 
  # the last two generations, this is when we should store parental 
  # information and trace pedigree info
  #
  # We do this at the end to save some time
  operators.extend([
    # save ancestral populations starting at -2 gen
    setAncestralDepth(2, at=[-2]),
    # track pedigree
    parentsTagger(begin=-2),
  ])
  #
  # Now, most troublesome operators: dealing with output.
  # We need to write expressions that will be evaluated automatically by pyEval
  # operator. 
  # 
  # ** If you are not sure you are getting the right expression (expression error
  #    happend. Use a output() operator to output the expression itself.
  # 
  # We need to track:
  # 1. LD (D') from a central marker to all others on a
  #   non-DSL chromosome, at burnin and before and after migration
  # 2. LD (D') from a central DSL to all others,
  #   at burnin and end
  # 3. Fst before and after migration
  # 4. allele frequency at DSL at the end (always monitored)
  # 5. Mean observed heterogeneity at the marker loci
  #
  #
  DSLLDString =  r"\nD' between DSL %d (chrom %d) and surrounding markers at gen %%d\n" \
    % (ctrDSL, ctrChrom)
  nonDSLLDString = r"\nD' between a center marker %d (chrom %d) and surrounding markers at gen %%d\n" \
    % (gt.chromBegin(noDSLChrom)+numLoci/2, noDSLChrom)
  #
  # each DSL is available as LD_prime[a][b] but we have so many
  # of them, need to write some expressions (that will be evaluated automatically)
  # ''"" etc quickly get complicated
  #
  # What we need is something like this:
  # "%.4f %.4f \n" % (LD_prime[1][2], LD_Prime[3][4])
  # 
  ctrDSLExpr = '"' + ('%.4f ' * len(ctrDSLLD))  + r'\n" % (' 
  for ld in ctrDSLLD:
    ctrDSLExpr += ' LD_prime[%d][%d],' % (ld[0], ld[1])
  ctrDSLExpr += ')'
  #
  nonDSLExpr = '"' + ('%.4f ' * len(noDSLLD))  + r'\n" % ('
  for ld in noDSLLD:
    nonDSLExpr += ' LD_prime[%d][%d],' % (ld[0], ld[1])
  nonDSLExpr += ')'
  #
  # output a table of DSL allele freq like
  # all x x x
  # 1   x x x
  # 2   x x x
  # ...
  alleleFreqExpr = r'"all\t' + (r'%.4f\t' * numDSL) + r'\n'
  for sp in range(numSubPop):
    alleleFreqExpr += str(sp+1) + r'\t' + (r'%.4f\t' * numDSL) + r'\n'
  alleleFreqExpr += '" % ('
  for d in DSL:
    alleleFreqExpr += ' 1.-alleleFreq[%d][1],' % d
  for sp in range(numSubPop):
    for d in DSL:
      alleleFreqExpr += '1.-subPop[%d]["alleleFreq"][%d][1],' % (sp, d)
  alleleFreqExpr += ')'
  # need reduce(....)
  heteroFreqExec = "AvgHetero=("
  for d in range(gt.totNumLoci()):
    heteroFreqExec += 'heteroFreq[%d][0]+' % d
  heteroFreqExec += '0.)/%d' % gt.totNumLoci() 
  heteroFreqExpr = r"'\nAverage counted heterozygosity is %.4f.\n' % AvgHetero"
  #
  # pyEval operators
  operators.extend( [
    # calculate LD
    stat( LD = ctrDSLLD + noDSLLD, at = [burnin, mixing, endGen]),
    # calculate average Fst based on non-DSL markers
    # calculate observed heterogeneity at all markers (including DSL)
    stat( Fst = nonDSL, heteroFreq = range(gt.totNumLoci()), at = [ mixing, endGen]),
    #
    # output statistics
    # LD
    pyEval( '"' + DSLLDString + '" % gen', output = ">>" + logFile,
      at = [split, mixing, endGen]),
    pyEval( ctrDSLExpr,  output = ">>" + logFile, at = [split, mixing, endGen]),
    pyEval( '"' + nonDSLLDString + '" % gen', output = ">>" + logFile,
      at = [split, mixing, endGen]),
    pyEval( nonDSLExpr,  output = ">>" + logFile, at = [split, mixing, endGen]),
    # allele frequency
    pyEval(r'"\nDSL allele frequency at gen %d\n" % gen', output = ">>" + logFile,
      at = [mixing, endGen]), 
    pyEval( alleleFreqExpr, output = ">>" + logFile, at = [mixing, endGen]),
    # heterogeneity
    pyExec( heteroFreqExec, at = [ mixing, endGen]),
    pyEval( heteroFreqExpr, output = ">>" + logFile, at = [ mixing, endGen]),
    # Fst
    pyEval( r'"\nFst at gen %d is %.4f\n" % (gen, AvgFst)', output = ">>" + logFile,
      at = [ mixing, endGen]),
    # monitor progress and 
    # see how long the simulation has been running
    output("Start introducing disease\t\t", at = [burnin] ),
    output("Start no-migration stage\t\t", at = [split ] ),
    output("Start mixing\t\t", at = [mixing]),
    output("End of simulation. \n", at = [endGen]),
    pyEval(r'"Generation %d\n" % gen ', step = 100),
    # Show how long the program has been running.
    ticToc(at=[burnin, split, mixing, endGen]),
    ]
  )
  # with all operators, we can set up a simulator
  # 
  # produce two offsprings only at the last generations.
  def last_two(gen):
    if gen >= endGen -2:
      return 2
    else:
      return 1
  #
  # may need to run several times to
  # get a valid population (requirements on allele frequency
  fixedCount = 1
  while(True):
    # create a simulator
    simu = simulator(
      population(subPop=popSizeFunc(0), ploidy=2,
        loci = [gt.numLoci(x) for x in range(gt.numChrom())],
        maxAllele=gt.maxAllele(),
        lociDist=[gt.locusDist(x) for x in range(gt.totNumLoci())]), 
      randomMating(newSubPopSizeFunc=popSizeFunc,
        numOffspringFunc=last_two),
      rep=1)
    #
    # evolve! If --dryrun is set, only show info
    simu.evolve( preOps = preOperators, ops = operators, end=endGen,
      dryrun=dryrun )
    if dryrun:
      raise exceptions.SystemError("Stop since in dryrun mode.")
    # if succeed
    if simu.gen() == endGen + 1:  # if not fixed
      pop = simu.population(0)
      # we want to save info on how this population is generated.
      # This is not required but is a good practise
      pop.dvars().DSL = DSL
      pop.dvars().DSLAfter = DSLafter
      pop.dvars().DSLdist = DSLdist
      pop.dvars().initSize = initSize
      pop.dvars().meanInitAllele = meanInitAllele
      pop.dvars().finalSize = finalSize
      pop.dvars().burnin = burnin
      pop.dvars().introGen = introGen
      pop.dvars().numSubPop = numSubPop
      pop.dvars().noMigrGen = noMigrGen
      pop.dvars().mixingGen = mixingGen
      pop.dvars().minAlleleFreq = minAlleleFreq
      pop.dvars().maxAlleleFreq = maxAlleleFreq
      pop.dvars().mutationRate = mu
      pop.dvars().mutationModel = "symmetric stepwise"
      pop.dvars().migrationRate = mi
      pop.dvars().migrationModel = "circular step stone"
      pop.dvars().recombinationRate = rec
      pop.dvars().numLoci = numLoci
      # this is for the convenience of plotting LD values
      pop.dvars().ctrDSL = ctrDSL
      pop.dvars().ctrChrom = ctrChrom
      pop.dvars().ctrDSLLD = ctrDSLLD
      pop.dvars().noDSLChrom = noDSLChrom
      pop.dvars().noDSLLD = noDSLLD
      # can not just return the reference...
      # return a copy (and then simulator will be destroyed)
      return simu.getPopulation(0)
    else:
      print "Population restarted at gen ", simu.gen()
      print "Overall fixed population ", fixedCount
      print "Allelefreq ", ( "%.3f " * numDSL + "\n\n") % \
        tuple([1-simu.dvars(0).alleleFreq[x][1] for x in DSL]) 
      fixedCount += 1

#
# plot the LD values for the sample.
def plotLD(pop, epsFile, jpgFile):
  ''' plot LD values in R and convert to jpg if possible '''
  r.postscript(file=epsFile)
  r.par(mfrow=[2,1])
  # return max LD
  res = []
  # dist: distance (location) of marker
  # ldprime: D' value
  dist = []
  ldprime = []
  for ld in pop.dvars().ctrDSLLD:
    if ld[1] == pop.dvars().ctrDSL:
      dist.append(pop.locusDist(ld[0]))
    else:
      dist.append(pop.locusDist(ld[1]))
    ldprime.append(pop.dvars().LD_prime[ld[0]][ld[1]])
  res.append(max(ldprime))
  r.plot( dist, ldprime, main="D' between DSL and other markers on chrom %d" % pop.dvars().ctrChrom,
    xlab="marker location", ylab="D'", type='b')
  r.abline( v = pop.locusDist(pop.dvars().ctrDSL), lty=3 )
  r.axis( 1, [pop.locusDist(pop.dvars().ctrDSL)], ['DSL'])
  dist = []
  ldprime = []
  for ld in pop.dvars().noDSLLD:
    if ld[1] == pop.chromBegin(pop.dvars().noDSLChrom) + numLoci/2:
      dist.append(pop.locusDist(ld[0]))
    else:
      dist.append(pop.locusDist(ld[1]))
    ldprime.append(pop.dvars().LD_prime[ld[0]][ld[1]])    
  res.append(max(ldprime))
  r.plot( dist, ldprime, main="D' between marker %d and other markers on chrom %d" \
    % (numLoci/2, pop.dvars().noDSLChrom),
    xlab="marker location", ylab="D'", type='b')    
  r.abline( v = pop.locusDist(pop.chromBegin(pop.dvars().noDSLChrom)+pop.dvars().numLoci/2), lty=3 )
  r.dev_off()
  # try to get a jpg file
  try:
    if os.system("convert -rotate 90 %s %s " % (epsFile, jpgFile) ) == 0:
      return (2,res)  # success
    else:
      return (1,res) # fail
  except:
    return (1,res)  # fail

# penetrance functions. They will return a penetrance function
# 
def recessive1(pen):
  ''' return a penetrance function (weird? :-) '''
  def func(geno):
    if geno[0] + geno[1] == 4:     # only use the first DSL
      return pen
    else:
      return 0
  return func

def recessive2(pen):
  ''' recessive All '''
  def func(geno):
    val = 1
    numRecessive = 0
    if len(geno) < 4:
      raise exceptions.ValueError("Number of DSL < 2, can not apply this penetrance function.")
    for i in range(2):
      if geno[i*2]+geno[i*2+1] == 4:
        val *= 1 - pen
        numRecessive += 1
    if numRecessive == 0:
      return 0
    else:
      return 1-val
  return func

def recessive3(pen):
  ''' recessive All '''
  def func(geno):
    val = 1
    numRecessive = 0
    if len(geno) < 6:
      raise exceptions.ValueError("Number of DSL < 3, can not apply this penetrance function.")
    for i in range(3):
      if geno[i*2]+geno[i*2+1] == 4:
        val *= 1 - pen
        numRecessive += 1
    if numRecessive == 0:
      return 0
    else:
      return 1-val
  return func
  
def additive1(pen):
  ''' additive2 '''
  def func(geno):
    return (geno[0]+geno[1]-2)*pen/2.
  return func

def additive2(pen):
  ''' additive2 '''
  def func(geno):
    if len(geno) < 4:
      raise exceptions.ValueError("Number of DSL < 2, can not apply this penetrance function.")
    if geno.count(2) < 2:
      return 0.
    # other wise a heterogeneity model
    val = 1
    for i in range(2):
      val *= 1- (geno[i*2]+geno[i*2+1]-2)*pen/2.
    return 1-val
  return func
  
def additive3(pen):
  ''' additive2 '''
  def func(geno):
    if len(geno) < 6:
      raise exceptions.ValueError("Number of DSL < 3, can not apply this penetrance function.")
    if geno.count(2) < 3:
      return 0.
    # other wise a heterogeneity model
    val = 1
    for i in range(3):
      val *= 1- (geno[i*2]+geno[i*2+1]-2)*pen/2.
    return 1-val
  return func


def drawSamples(pop, penFun):
  ''' get samples of different type using a penetrance function '''
  # first, apply peneFunction
  PyPenetrance(pop, loci=DSL, func=penFun)
  #
  # all types of samples
  allSample = []
  #
  # 1. population based case control
  # get number of affected
  Stat(pop, numOfAffected=True)
  nCase = min(pop.dvars().numOfAffected , N/2)
  nControl = min(pop.dvars().numOfUnaffected, N/2)
  try:
    # if N=800, 400 case and 400 controls
    s = CaseControlSample(pop, nCase, nControl)
    # remove DSL
    s[0].removeLoci(remove=DSL)
    allSample.append(s[0])
  except:
    print "Can not draw case control sample. "
    allSample.append(None)
  #
  # 2. affected and unaffected sibpairs
  # this is difficult since simuPOP only has
  # methods to draw one kind of sibpairs and 
  try:
    # get number of affected/unaffected sibpairs
    # There may not be enough to be sampled
    AffectedSibpairSample(pop, countOnly=True)
    nAff = min(pop.dvars().numAffectedSibpairs, N/4)
    AffectedSibpairSample(pop, chooseUnaffected=True, countOnly=True)
    nUnaff = min(pop.dvars().numAffectedSibpairs, N/4)
    #
    affected = AffectedSibpairSample(pop, name='sample1',
       size=nAff)
    # now chose unaffected. These samples will not overlap
    # 
    # NOTE: however, you may not be able to easily merge these two 
    # samples since they may shared parents.
    #
    # Use another name to avoid conflict since these sampled are stored
    # in local namespace
    unaffected = AffectedSibpairSample(pop, chooseUnaffected=True,
      name='sample2', size=nUnaff)
    # remove DSL
    affected[0].removeLoci(remove=DSL)
    unaffected[0].removeLoci(remove=DSL)
    allSample.extend([affected[0], unaffected[0] ])
  except:
    print "Can not draw affected sibpars."
    allSample.extend([None, None])
  return allSample


def TDT(dataDir, data, epsFile, jpgFile):
  ''' use TDT method to analyze the results. Has to have rpy installed '''
  if not hasRPy:
    return (0,[])
  # write a batch file and call gh
  allPvalue = []
  print "Applying TDT method to affected sibpairs "
  for ch in range(numChrom):
    inputfile = dataDir+data+ "_%d" % ch
    if not os.path.isfile(inputfile + ".ped"):
      print "Ped file ", inputfile+".ped does not exist. Can not apply TDT method."
      return (0,[])
    # batch file
    f=open("ghTDT.cmd","w")
    f.write("load markers " + inputfile + ".dat\n")
    f.write("tdt " + inputfile + ".ped\n")
    f.close()
    # run gene hunter
    os.system(geneHunter + " < ghTDT.cmd > res.txt ")
    # get the result
    # ( I was using              
    # | grep '^loc' | tr -d 'a-zA-Z+-' > " + outputfile)
    # but this is not portable
    #
    # use re package
    # get only loc number and p-value
    scan = re.compile('loc(\d+)\s+- Allele \d+\s+\d+\s+\d+\s+[\d.]+\s+([\d.]+)\s*.*')
    minPvalue = [1]*(numLoci-1)
    try:
      res = open("res.txt")  # read result
      for l in res.readlines():
        try:
          # get minimal p-value for all alleles at each locus
          (loc,pvalue) = scan.match(l).groups()
          if minPvalue[int(loc)-1] > float(pvalue):
            minPvalue[int(loc)-1] = float(pvalue)
        except:
          # does not match
          continue
      res.close()
    except:
      print "Can not open result file. TDT failed"
      return (0,[])
    # collect -log10 pvalue
    allPvalue.extend([-math.log10(max(x,1e-6)) for x in minPvalue])
    # There will be numLoci-1 numbers, pad to the correct length
    allPvalue.append(0)
  try:
    os.unlink('res.txt')
    os.unlink('ghTDT.cmd')
  except:
    pass
  # now, we need to see how TDT works with a set of p-values around DSL
  # DSL is global
  res = []
  for d in DSL: 
    res.append( max(allPvalue[(d-2):(d+2)]))
  # use R to draw a picture
  r.postscript(file=epsFile)
  r.plot(allPvalue, main="-log10(P-value) for each marker (TDT)",
    xlab="chromosome", ylab="-log10 p-value", type='l', axes=False)
  r.box()
  r.abline( v = [DSLafter[g]+DSLdist[g] for g in range(len(DSL))], lty=3)
  r.abline( h = -math.log10(0.05))                       
  r.axis(1, [numLoci*x for x in range(numChrom)], [str(x+1) for x in range(numChrom)])
  r.axis(2)
  r.dev_off()
  # try to get a jpg file
  try:
    if os.system("convert -rotate 90 %s %s " % (epsFile, jpgFile) ) == 0:
      return (2,res)  # success
    else:
      return (1,res) # fail
  except:
    return (1,res)  # fail

# create output directory if necessary
# a more friendly version of mkdie
def _mkdir(d):
  try:
    if not os.path.isdir(d):
      os.makedirs(d)
  except:
    print "Can not make output directory ", d
    sys.exit(1)

def processOnePopulation(dataDir, numChrom, numLoci, markerType,
    DSLafter, DSLdist, initSize, meanInitAllele, burnin, introGen, introSel, minAlleleFreq,
    maxAlleleFreq, fitness, mlSelModel, numSubPop, finalSize, noMigrGen,
    mixingGen, popSizeFunc, migrModel, mu, mi, rec, dryrun, popIdx):
  '''
     this function organize all previous functions
     and
     1. generate (or load) a population
     2. apply different kinds of penetrance functions
     3. draw sample
     4. save samples
     5. apply TDT and/or randTent method
     6. return a text summary and a result list
  '''
  # get population
  result = [popIdx, mu, mi, rec]
  logFile = dataDir + "/pop_" + str(popIdx) + ".log"
  popFile = dataDir + "/pop_" + str(popIdx) + ".bin"
  genDataset = True
  if os.path.isfile(popFile) and (not overwrite):
    print "Loading a pre-existing file ", popFile
    pop = LoadPopulation(popFile)
    if abs(pop.dvars().mutationRate - mu) + abs(pop.dvars().migrationRate - mi) \
       + abs( pop.dvars().recombinationRate - rec) < 1e-7:
      genDataset = False
      print "Dataset already exists, load it directly. (use --overwrite True to change this behavior)"
    else:
      print "A dataset is present but was generated with different parameters. Re-generating."
  if genDataset:
    print "Generating dataset ", str(popIdx)
    pop = simuComplexDisease( numChrom, numLoci, markerType, DSLafter, DSLdist,
      initSize, meanInitAllele, burnin, introGen, introSel, minAlleleFreq,
      maxAlleleFreq, fitness, mlSelModel, numSubPop, finalSize, noMigrGen,
      mixingGen, popSizeFunc, migrModel, mu, mi, rec, dryrun, logFile)
    SavePopulation(pop, popFile)
  # have the population, get info
  result.extend( [ pop.dvars().AvgFst, pop.dvars().AvgHetero])
  # log file is ...
  summary = '''<h2><a name="pop_%d">Dataset %d</a></h2>\n''' \
    % (popIdx, popIdx)
  summary += '''<h3>Log file (LD and other statistics): 
    <a href="pop_%d/pop_%d.log">pop_%d.log</a></h3>\n''' \
    % (popIdx, popIdx, popIdx)
  #
  # actually write the log file in the summary page
  try:
    lf = open(logFile)
    summary += "<pre>\n"+lf.read()+"\n</pre>\n"
    lf.close()
  except:
    print "Can not open log file, ignoring. "
  #
  # plot LD, res = 0, fail, 1: eps file, 2: converted to jpg
  epsFile = dataDir + "/LD_" + str(popIdx) + ".eps"
  jpgFile = dataDir + "/LD_" + str(popIdx) + ".jpg"
  #
  if hasRPy:
    (suc,ldres) = plotLD(pop, epsFile, jpgFile)
    if suc > 0 : # eps file successfully generated
      summary += """<p>D' measures on two chromosomes with/without DSL at the last gen: 
      <a href="pop_%d/LD_%d.eps">LD.eps</a></p>\n""" % (popIdx, popIdx)
    if suc > 1 : # jpg file is generated
      summary += '''<img src="pop_%d/LD_%d.jpg"
        width=800, height=600>'''  % (popIdx, popIdx)
    result.extend(ldres)
  result.extend( [1- pop.dvars().alleleFreq[i][1] for i in DSL])
  #
  # apply penetrance and get three samples (different format)
  summary += "<h3>Samples using different penetrance function</h3>\n"
  # if we are going to save in linkage format, and specify
  # allele frequency from the whole population, we need to calculate them
  # now. (Previously, we only have data for DSL
  if 'Linkage' in saveFormat:
    Stat(pop, alleleFreq=range(pop.totNumLoci()))
  for p in range(len(peneFunc)):
    if peneFunc[p] == 'recessive1':
      print "Using recessive1 penetrance function"
      summary += "<h4>Recessive at one locus</h4>\n"
      s = drawSamples(pop, recessive1( penePara[p]))
    elif peneFunc[p] == 'recessive2':
      print "Using recessive2 penetrance function"
      summary += "<h4>Recessive at one locus, 2 DSL</h4>\n"
      s = drawSamples(pop, recessive2(penePara[p]))
    elif peneFunc[p] == 'recessive3':
      print "Using recessive3 penetrance function"
      summary += "<h4>Recessive at one locus, 3 DSL</h4>\n"
      s = drawSamples(pop, recessive3(penePara[p]))
    elif peneFunc[p] == 'additive1':
      print "Using additive1 penetrance function"
      summary += "<h4>Additive at one locus</h4>\n"
      s = drawSamples(pop, additive1(penePara[p]))
    elif peneFunc[p] == 'additive2':
      print "Using recessive3 penetrance function"
      summary += "<h4>Additive at 2 DSL, at least 2 DS alleles</h4>\n"
      s = drawSamples(pop, additive2(penePara[p]))
    else:
      print "Using additive2 penetrance function"
      summary += "<h4>Additive at 3 DSL, at least 3 DS alleles</h4>\n"
      s = drawSamples(pop, additive3(penePara[p]))
    # save these samples
    #
    penDir = dataDir + "/" + peneFunc[p]
    relDir = 'pop_%d/%s/' % (popIdx, peneFunc[p])
    _mkdir(penDir)
    summary += "<p>Case-control, affected and unaffected sibpairs saved in different formats. Sample sizes are "
    if s[0] != None:
      summary += str(s[0].popSize()) + " (case-control), "
    if s[1] != None:
      summary += str(s[1].popSize()) + " (affected sibs), "
    if s[2] != None:
      summary += str(s[2].popSize()) + " (unaffected sibs) "
    summary += '<ul>'
    for format in saveFormat:
      if format == 'simuPOP':
        print "Write to simuPOP format"
        summary += '<li>simuPOP binary format:'
        if s[0] != None: # has case control
          binFile = penDir + "/caseControl.bin"
          s[0].savePopulation(binFile)
          summary += '<a href="%scaseControl.bin"> caseControl.bin</a>, ' % relDir
        if s[1] != None: # has affected sibpairs
          binFile = penDir + "/affectedSibpairs.bin"
          s[1].savePopulation(binFile)
          summary += '<a href="%saffectedSibpairs.bin"> affectedSibpairs.bin</a>, ' % relDir
        if s[2] != None:
          binFile = penDir + "/unaffectedSibpairs.bin"
          s[2].savePopulation(binFile)
          summary += '<a href="%sunaffectedSibpairs.bin"> unaffectedSibpirs.bin</a>, ' % relDir
        summary += '</li>\n'
      if format == 'Linkage':
        summary += '<li>Linkage format by chromosome:'
        linDir = penDir + "/Linkage/"
        _mkdir(linDir)
        if s[1] != None: # has case control
          for ch in range(0,pop.numChrom() ):
            SaveLinkage(pop=s[1], popType='sibpair', output = linDir+"/Aff_%d" % ch,
              chrom=ch, recombination=pop.dvars().recombinationRate,
              alleleFreq=pop.dvars().alleleFreq, daf=0.1)        
          summary +=  '<a href="%sLinkage">affected</a>, ' % relDir
        if s[2] != None:
          for ch in range(0,pop.numChrom() ):
            SaveLinkage(pop=s[2], popType='sibpair', output = linDir+"/Unaff_%d" % ch,
              chrom=ch, recombination=pop.dvars().recombinationRate,                            
              alleleFreq=pop.dvars().alleleFreq,  daf=0.1)        
          summary += '<a href="%sLinkage">unaffected</a>' % relDir
        summary += '</li>\n'
      if format == 'randTent':
        summary += '<li>randTent .csv format:'
        if s[1] != None: # has affected sibpairs
          SaveCSV(s[1], penDir+"/Aff_all.csv")
          summary +=  ( 'affected sibpairs <a href="%sAff_all.csv"> combined</a>, ' ) % relDir
        if s[2] != None:
          SaveCSV(s[2], penDir+"/Unaff_all.csv")
          summary +=  ( ' unaffected sibpairs <a href="%sUnaff_all.csv"> combined</a>' )  % relDir
        summary += '</li>\n'
    summary += '</ul>\n'
    # if there is a valid gene hunter program, run it
    (suc,res) = TDT(penDir, "/Linkage/Aff", penDir + "/TDT.eps", penDir + "/TDT.jpg")
    if suc > 0 : # eps file successfully generated
      summary += """<h4>TDT analysis for affected sibpair data:  <a href="%s/TDT.eps">TDT.eps</a>""" % relDir
    if suc > 1 : # jpg file is generated
      summary += '''<p><img src="%s/TDT.jpg" width=800, height=600></p>''' % relDir
    # keep some numbers depending on the penetrance model
    result.append(res[:int(peneFunc[p][-1])])
  return (summary, result)


def writeReport(content, allParam, results):
  ''' write a HTML file. The parts for each population has
    been written but we need a summary table. '''
  print "Writing a report (saved in summary.htm )"
  try:
    summary = open(outputDir + "/summary.htm", 'w')
  except:
    raise exceptions.IOError("Can not open a summary file : " + outputDir + "/summary.htm to write")
  summary.write('''
  <HTML>
  <HEAD>
  <TITLE>Summary of simulation</TITLE>
  <META NAME="description" CONTENT="summary of simulation">
  <META NAME="keywords" CONTENT="simuPOP">
  <META NAME="resource-type" CONTENT="document">
  <META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=iso-8859-1">
  </HEAD>
  <BODY >
  <h1>Output of %s</h1>
  <p><b>Time of simulation: </b> %s </p>
  <h2>Parameters </h2>
  <ul>''' % (sys.argv[0], time.asctime()) )
  # write out parameters
  # start from options[1]
  for i in range(1, len(options)-1):
    # check type
    if options[i].has_key('configName'):
      if type( allParam[i-1]) in [types.TupleType, types.ListType]:
        summary.write('''<li><b>%s:</b>    %s</li>\n''' % \
          (options[i]['configName'], 
          ', '.join( map(str, allParam[i-1]))))
      else:
        summary.write('''<li><b>%s:</b>    %s</li>\n''' % \
          (options[i]['configName'], 
          str(allParam[i-1])))
  # a table built from res which has
  # idx, muta, migr, rec, af?, Fst, Het, TDT?
  summary.write('''
  </ul>
  <h2>Summary of datasets </h2>
  <p>The following table lists population id, mutation rate, migration
  rate, recombination rate, Fst, average heterozygosity, highest D' between a 
  DSL and all surrounding markers (not necessarily its cloest marker), highest
  D' between a marker with its surrounding markers on a chromsome without
  DSL, allele frequency at DSL, -log10 p-values (TDT method) at
  all relevant DSL. </p>
  <table border="1">
  <tr><th>id </th>
  <th>mu</th>  <th>mi</th>   <th>rec</th>
  <th>Fst</th> <th>Het</th> <th>D'(dsl)</th> <th>D'(non)</th>
  ''' +  ''.join(map(lambda x: '<th>allFrq%d</th>'%(x+1), \
    range(len(allParam[3])) )) )
  #
  # how about TDT?
  if hasRPy:
    summary.write(' '.join(map(lambda x: '<th>%s</th>'%x, allParam[20])))
  summary.write('</tr>')
  # end of headers, now result
  for res in results:
    summary.write('''<tr><td><a href="#pop_%d">%d</a></td> ''' \
      % (res[0], res[0]))    
    for i in range(1,len(res)):
      item = res[i]
      if type(item) == types.FloatType:
        summary.write('<td>%.3g</td>'% item)
      elif type(item) in [types.ListType, types.TupleType]:
        summary.write('<td>' + ', '.join(map(lambda x:'%.2g'%x, item))+'</td>')
      else:
        summary.write('<td>%s</td>' % str(item))
    summary.write('''</tr>''')
  # the middle (big) and last piece) 
  summary.write('''</table>
   %s 
   <h2>Usage of %s</h2>
   <pre>%s </pre>
   </BODY></HTML>
   ''' % ( content, sys.argv[0], simuOpt.usage(options, __doc__)))
  summary.close()


if __name__ == '__main__':
  allParam = getOptions()
  # unpack options
  (numChrom, numLoci, markerType, DSLafter, DSLdist,
    initSize, meanInitAllele, burnin, introGen, introSel,  minAlleleFreq, 
    maxAlleleFreq, fitnessTmp, mlSelModelTmp, 
    numSubPop, finalSize, noMigrGen,
    mixingGen, growth, migrModel, migrRate, mutaRate, recRate,
    numRep, saveFormat, peneFuncTmp, penePara, N, outputDir,
    overwrite, geneHunter, dryrun, saveConfig) = allParam
  #
  # this may not happen at all but we do need to be careful
  if initSize < len(DSLafter):
    raise exceptions.ValueError("Initial population size is too small. (Less than number of DSL)")
  #
  if len( DSLafter ) != len(DSLdist):
    raise exceptions.ValueError("Please specify DSL distance for each DSL.")

  fitness = []
  if type(fitnessTmp[0]) in [types.FloatType, types.IntType]:
    for d in DSLafter:
      # fitness is the same for all DSL
      fitness.append(fitnessTmp)
  else:
    if len(fitnessTmp) != len(DSL):
      raise exceptions.ValueError("Please specify fitness for each DSL")
    fitness = fitnessTmp
  for f in fitness:
    if len(f) != 3:
      raise exceptions.ValueError("Fitness values should be an array of length 3. Given " + str(f) )
     
  # load simuPOP libraries
  from simuPOP import *
  from simuUtil import *

  if mlSelModelTmp == 'additive':
    mlSelModel = SEL_Additive
  elif mlSelModelTmp == 'multiplicative':
    mlSelModel = SEL_Multiplicative
  elif mlSelModelTmp == 'none':
    mlSelModel = SEL_None
  else:
    raise exceptions.ValueError("Wrong multi-locus seleciton model. " + mlSelModelTmp)
    
  #
  # if RPy is available, several figures will be drawn 
  try:
    from simuRPy import *
  except:
    hasRPy = False
  else:
    hasRPy = True
  #
  # determine which function to use
  # define it here to allow more flexiblity of calling
  # simulate function.
  split  = burnin + introGen
  endGen = split + noMigrGen + mixingGen  
  if growth == 'linear':
    popSizeFunc = LinearExpansion(initSize, finalSize, endGen,
      burnin, split, numSubPop)
  elif growth == 'exponential':
    popSizeFunc = ExponentialExpansion(initSize, finalSize, endGen,
      burnin, split, numSubPop)
  else:
    raise exceptions.ValueError("Growth model can be one of linear and exponential. Given " + growth)
  #  
  # Penetrance related !
  peneFunc = []
  for pe in peneFuncTmp:
    if int(pe[-1]) > len(DSLafter):
      print "Number of DSL is ", len(DSLafter), " ignoring penetrance function ", pe
    else:
      peneFunc.append(pe)
      
  # check penePara
  if len(penePara) == 0:
    penePara = [0.25]*len(peneFunc)
  if len(peneFunc) > len(penePara):
    for i in range( len(peneFunc) - len(penePara)):
      penePara.append(penePara[0])
  #
  # everything is ready
  popIdx = 1
  summary = ''
  results = []
  for mu in mutaRate:
    for mi in migrRate:
      for rec in recRate:
        for rep in range(numRep):
          dataDir = outputDir + "/pop_" + str(popIdx)
          # outputDir should already exist
          _mkdir(dataDir)      
          (text, result) =  processOnePopulation(dataDir,
            numChrom, numLoci, markerType, 
            DSLafter, DSLdist, initSize, burnin, introGen, introSel,
            minAlleleFreq, maxAlleleFreq, fitness, mlSelModel, numSubPop, 
            finalSize, noMigrGen, mixingGen, popSizeFunc, migrModel, 
            mu, mi, rec, dryrun, popIdx)
          summary += text
          results.append( result)
          popIdx += 1
  writeReport(summary, allParam, results)
