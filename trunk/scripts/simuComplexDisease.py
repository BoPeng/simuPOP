#!/usr/bin/env python
#
# Purpose:  generate dataset for common complex disease 
#           with certain number of disease susceptibility
#           loci.
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedDate$
# $Rev$
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

Parameters about the lower and upper bound of acceptable allele 
frequency at each DSL can be given. During this stage, the intensity
of selective advantage will vary according to allele frequencies 
at individual DSL. For example, for a DSL that has reached or 
exceeded the upper bound, it will be subject to no or even negative
selection pressure. This mechanism will allow all DSL to reach 
their expected range at the end of this stage.

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
  3. LD (D prime) between markers and between a DSL and surrounding markers
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

import simuOpt, simuUtil
import os, sys, types, exceptions, os.path, re, math, time, copy

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
  # this parameter will not be shown in the dialog, or be prompted
  # It can only be given through comandline parameter.
  {'longarg': 'meanInitAllele=',
   'default': 50,
   'allowedTypes': [types.IntType, types.LongType],
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
        stage, large introGen will lead to over-common diseases. Note 
        that the larger the initial population size is, the longer 
        introGen is required.''',
   'validate':  simuOpt.valueGT(0),
  },
  {'longarg': 'minAlleleFreq=',
   'default': [0.05],
   'configName': 'minimal Allele Frequency',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Minimum allele frequency required for each DSL (0.05):  ',
   'description': '''Mininal allele frequencies required for each DSL.
        If a number is given, it is assumed to be the lower bound for all
        DSL.''',
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },
  {'longarg': 'maxAlleleFreq=',
   'default': [0.10],
   'configName': 'maximum Allele Frequency',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Maximum allele frequency required for each DSL (0.10):  ',
   'description': '''Maximum allele frequencies required for each DSL.
        If a number is given, it is assumed to be the upper bound for all
        DSL.''',
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
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
        will mix individuals from subpopulations using an circular stepping stone
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
   'default': 'stepping stone',
   'configName': 'Migration Model',
   'prompt': 'Mutation model. (stepping stone):  ',
   'allowedTypes': [types.StringType],
   'description': '''Migration model. Choose between stepping stone (circular),
        island and none. ''',
   'validate':  simuOpt.valueOneOf(['island', 'stepping stone', 'none']),
   'chooseOneOf': ['stepping stone', 'island', 'none']
  }, 
  {'longarg': 'migrRate=',
   'default': [0,0.001,0.01,0.1],
   'configName': 'Migration Rates',
   'prompt': 'Migration rate(s) during mixing stage. A separate dataset\n' +
        'will be genrated for each of the given migration rate. ([0, 0.01, 0.05, 0.1]):  ',
   'description': '''Migration rate during mixing stage. Can be a number or an array.
        A circular stepping stone migration model will be used. 
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
   'default': [0.0005],
   'configName': 'Recombination Rates',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Recombination rate between adjacent markers. Separate datasets\n'+
        'will be generated for each rec values if a list is given. (0.0005):  ',
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
   'description': ''' Penetrance functions to be applied to the final
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
  {'longarg': 'numSampleber=',
   'default': 2,
   'configName':  'Sample Number',
   'allowedTypes':  [types.IntType, types.LongType],
   'prompt':  'Number of samples for each penetrance function (2):  ',
   'description':  '''Number of samples to draw for each penetrance function. ''',
   'validate':  simuOpt.valueGT(0)
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
   'prompt': 'Provide location of gene hunter ():  ',
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

# __doc__ + parameter definitions use >500 lines, 
# more than one third of the total length of the script.

def getOptions(details=__doc__):
  ''' get options from options structure,
    if this module is imported, instead of ran directly,
    user can specify parameter in some other way.
  '''
  # get all parameters, __doc__ is used for help info
  allParam = simuOpt.getParam(options, 
    '''This program simulates the evolution of a complex common disease 
  under the effect of mutation, migration, recombination and population
  size change. Click 'help' for more info.''', details, nCol=2)
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
  if allParam[-1]:         # verbose
    for p in range(0, len(options)):
      if options[p].has_key('configName'):
        if type(allParam[p]) == types.StringType:
          print options[p]['configName'], '\t"'+str(allParam[p])+'"'
        else:
          print options[p]['configName'], '\t', str(allParam[p])
  # return the rest of the parameters
  return allParam[1:-1]

# we need to output various statistics, using pyEval operator
# like what I did in the previous versions are too tedious. 
# This version uses a pyOperator to do it all together.
def outputStatistics(pop, args): 
  ''' this function will be working with a pyOperator
    to output statistics. Many parameters will be passed,
    packed as a tuple.
  
   We need to output
   1. LD (D') from a central marker to all others on a
     non-DSL chromosome, at burnin and before and after migration
   2. LD (D') from a central DSL to all others,
     at burnin and end
   3. Fst before and after migration
   4. allele frequency at DSL 
   5. Mean observed heterogeneity at the marker loci
  '''
  # unwrap parameter
  (burnin, split, mixing, endGen, outfile) = args
  # 
  gen = pop.gen()
  # see how long the simulation has been running
  if gen == burnin:
    print "Start introducing disease\t\t"
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
  # In this analysis, we also need to know
  # 1. a DSL that is closest to the center of a chromosome
  toCtrDist = [abs(pop.locusDist(x)-numLoci/2) for x in DSL]
  ctrChromDSL = DSL[ [x==min(toCtrDist) for x in toCtrDist].index(True)]
  ctrChrom = pop.chromLocusPair(ctrChromDSL)[0]
  # 2. first chromosome without DSL
  noDSLChrom = [pop.numLoci(x)==numLoci for x in range(pop.numChrom())].index(True)
  #
  # loci pairs on the selected chromosomes when calculating LD
  i = pop.chromBegin( ctrChrom)
  ctrDSLLD = [ [x, ctrChromDSL] for x in range(i, ctrChromDSL)] + \
    [ [ctrChromDSL, x] for x in range(ctrChromDSL+1, i+numLoci+1)]
  i = pop.chromBegin( noDSLChrom)
  noDSLLD = [ [i+x, i+numLoci/2] for x in range(numLoci/2)] + \
    [ [i + numLoci/2, i+x] for x in range(numLoci/2+1, numLoci)]
  # save these info for later use (a plotting function will
  # use these info.
  pop.dvars().ctrChrom = ctrChrom
  pop.dvars().ctrChromDSL = ctrChromDSL
  pop.dvars().ctrDSLLD = ctrDSLLD
  pop.dvars().noDSLChrom = noDSLChrom
  pop.dvars().noDSLLD = noDSLLD
  #
  #@@@ NOW, output statistics
  # append to output file.
  output = open(outfile, 'a')
  # first, calculate LD and other statistics
  Stat(pop, alleleFreq=DSL, LD = ctrDSLLD + noDSLLD, Fst = nonDSL, 
    heteroFreq = range( pop.totNumLoci() ) )
  # output D', allele frequency at split, mixing and endGen
  if gen in [split, mixing, endGen]:
    print >> output, "Average Fst estimated from non-DSL at gen %d: %.4f \n" % (gen, pop.dvars().AvgFst)
    print >> output, "D between DSL %d (chrom %d) and surrounding markers at gen %d" \
      % (ctrChromDSL, ctrChrom, gen)
    for ld in ctrDSLLD:
      print >> output, '%.4f ' % pop.dvars().LD[ld[0]][ld[1]],
    print >> output, "\n\nD between a center marker %d (chrom %d) and surrounding markers at gen %d" \
      % (pop.chromBegin(noDSLChrom)+numLoci/2, noDSLChrom, gen)
    for ld in noDSLLD:
      print >> output, '%.4f ' % pop.dvars().LD[ld[0]][ld[1]],
    print >> output, "\n\nD' between DSL %d (chrom %d) and surrounding markers at gen %d" \
      % (ctrChromDSL, ctrChrom, gen)
    for ld in ctrDSLLD:
      print >> output, '%.4f ' % pop.dvars().LD_prime[ld[0]][ld[1]],
    print >> output, "\n\nD' between a center marker %d (chrom %d) and surrounding markers at gen %d" \
      % (pop.chromBegin(noDSLChrom)+numLoci/2, noDSLChrom, gen)
    for ld in noDSLLD:
      print >> output, '%.4f ' % pop.dvars().LD_prime[ld[0]][ld[1]],
    print >> output, "\n\nAllele frequencies\nall\t",
    for d in DSL:
      print >> output, '%.4f ' % (1. - pop.dvars().alleleFreq[d][1]),
    for sp in range(numSubPop):
      print >> output, "\n%d\t" % sp,
      for d in DSL:
        print >> output, '%.4f ' % (1. - pop.dvars(sp).alleleFreq[d][1]),
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

# dynamic advantageous selector
def dynaAdvSelector(pop):
  ''' this selector will select give advantage to DSl according to
    allele frequency at each DSL. minAlleleFreq and maxAlleleFreq
    are stored in pop.
  '''
  DSL = pop.dvars().DSL
  # get allele frequencies
  Stat(pop, alleleFreq=DSL)
  # gives 1,1.25,1.5 to promote allele if freq < lower cound
  # gives 1,0.9,0.8 to select agsinst DSL with freq > upper
  sel = []
  freq = pop.dvars().alleleFreq
  # print +- symbol for each DSL to visualize how frequencies are manipulated
  for i in range(len(DSL)):
    # positive selection (promote allele)
    if 1-freq[DSL[i]][1] < pop.dvars().minAlleleFreq[i]:
      print '+',
      sel.append( maSelector(locus=DSL[i], wildtype=[1], fitness=[1,1.5,2]) )
    # negative selection (select against allele)
    elif 1-freq[DSL[i]][1] > pop.dvars().maxAlleleFreq[i]:
      print '-',
      sel.append( maSelector(locus=DSL[i], wildtype=[1], fitness=[1,0.9,0.8]) )
    else:
    # encourage slightly towards upper bound
      sel.append( maSelector(locus=DSL[i], wildtype=[1], fitness=[1,1.02,1.04]) )
      print ' ',
    # apply multi-locus selector, note that this opertor will only
    # set a variable fitness in pop, actual selection happens during mating.
    if len(sel ) > 0:  # need adjustment (needed if 'else' part is empty)
      MlSelect( pop, sel, mode=SEL_Multiplicative)
  return True

# simulate function, using a single value of mutation, migration,
# recombination rate
def simuComplexDisease( numChrom, numLoci, markerType, DSLafter, DSLdist,
    initSize, meanInitAllele, burnin, introGen, minAlleleFreq,
    maxAlleleFreq, fitness, mlSelModel, numSubPop, finalSize, noMigrGen,
    mixingGen, popSizeFunc, migrModel, mu, mi, rec, dryrun, logFile):
  ''' run a simulation of complex disease with given parameters. 
  '''
  if markerType == 'microsatellite':
    maxAle = 99                   # max allele
  else:
    maxAle = 2                    # SNP 
  #### translate numChrom, numLoci, DSLafterLoci to
  #### loci, lociDist, DSL, nonDSL in the usual index
  #### (count markers along with DSL)
  numDSL = len(DSLafter)
  if len(DSLdist) != numDSL:
    raise exceptions.ValueError("--DSL and --DSLloc has different length.")
  # real indices of DSL
  DSL = list(DSLafter)
  for i in range(0,numDSL):
    DSL[i] += i+1
  # adjust loci and lociDist
  loci = [0]*numChrom
  lociDist = []
  i = 0  # current absolute locus 
  j = 0  # current DSL
  for ch in range(0, numChrom):
    lociDist.append([])
    for loc in range(0,numLoci):
      # DSL after original loci indices
      if j < numDSL and i == DSLafter[j]:
        loci[ch] += 2
        lociDist[ch].append(loc+1)
        lociDist[ch].append(loc+1 + DSLdist[j])
        j += 1
      else:
        loci[ch] += 1
        lociDist[ch].append(loc+1)
      i += 1
  # non DSL loci, for convenience
  nonDSL = range(0, sum(loci))
  for loc in DSL:
    nonDSL.remove(loc)
  #### Change generation duration to point  
  split  = burnin + introGen
  mixing = split  + noMigrGen
  endGen = mixing + mixingGen  
  #### initialization and mutation
  if maxAle > 2:  # Not SNP
    preOperators = [
      # initialize all loci with 5 haplotypes
      initByValue(value=[[x]*sum(loci) \
        for x in range(meanInitAllele-2, meanInitAllele+3)],
        proportions=[.2]*5), 
      # and then init DSL with all wild type alleles
      initByValue([1]*len(DSL), atLoci=DSL)
    ]
    # symmetric mutation model for microsatellite
    mutator = smmMutator(rate=mu, maxAllele=maxAle,
      atLoci=nonDSL)
  else: # SNP
    preOperators = [
      # initialize all loci with two haplotypes (111,222)
      initByValue(value=[[x]*sum(loci) for x in range(1,3)],
        proportions=[.5]*2), 
      # and then init DSL with all wild type alleles
      initByValue([1]*len(DSL), atLoci=DSL)
    ]      
    # k-allele model for mutation of SNP
    mutator = kamMutator(rate=mu, maxAllele=2,
      atLoci=nonDSL)
  #### burn in stage
  # mutation and recombination will always be in effective
  # so we do not have to specify them later
  operators = [
    # mutator, differ by marker type
    mutator, 
    # recombination, use intensity since loci (include DSL)
    # are not equally spaced 
    recombinator(intensity=rec),
  ]
  #### disease introduction stage
  operators.extend( [
    # need to measure allele frequency at DSL
    stat(alleleFreq=DSL, popSize=True, begin = burnin),
    pyEval( expr=r'"%d(%d): "%(gen, popSize) + " ".join(["%.3f"%(1-alleleFreq[x][1]) for x in DSL])+"\n"', 
      begin = burnin) ]
  )
  # need five conditional point mutator to
  # introduce disease. It does not matter that we introduce
  # mutant specifically to individuel i since individuals
  # are unordered
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
  # use selective pressure to let allele frequency reach 
  # [minAlleleFreq, maxAlleleFreq].
  operators.append(
    pyOperator( func=dynaAdvSelector, stage=PreMating, 
      begin=burnin, end=split) )
  # check if allele frequency has reached desired value (loosely judged by 4/5*min)
  for i in range(len(DSL)):
    operators.append(
      terminateIf("1-alleleFreq[%d][1]<%f" % (DSL[i], minAlleleFreq[i]*4/5.),
      at=[split]) )
  #### no migration stage
  if numSubPop > 1:
    operators.append( 
      # split population after burnin, to each sized subpopulations
      splitSubPop(0, proportions=[1./numSubPop]*numSubPop, at=[split]),
      )
  #
  # start selection
  if mlSelModel != SEL_None:
    operators.append( mlSelector(
      # with five multiple-allele selector as parameter
      [ maSelector(locus=DSL[x], wildtype=[1], 
        fitness=fitness[x]) for x in range(len(DSL)) ],
      mode=mlSelModel, begin=split),
    )
  #
  ### mixing stage
  # a migrator, stepping stone or island
  if numSubPop > 1 and migrModel == 'island' and mi > 0:
    operators.append( migrator(migrIslandRates(mi, numSubPop),
      mode=MigrByProbability, begin=mixing) )
  elif numSubPop > 1 and migrModel == 'stepping stone' and mi > 0:
    operators.append( migrator(migrSteppingStoneRates(mi, numSubPop, circular=True),
      mode=MigrByProbability, begin=mixing) )
  #
  # prepare for sampling:
  # use last_two numOffspringFunc, simuPOP will produce 2 offspring at 
  # the last two generations, this is when we should store parental 
  # information and trace pedigree info
  #
  operators.extend([
    # save ancestral populations starting at -2 gen
    setAncestralDepth(2, at=[-2]),
    # track pedigree
    parentsTagger(begin=-2),
    # terminate simulation is on DSL get lost.
    terminateIf("True in [alleleFreq[x][1] == 1 for x in DSL]",
      begin=split), 
    # output statistics
    pyOperator(func=outputStatistics, 
      param = (burnin, split, mixing, endGen, logFile),
      at = [burnin, split, mixing, endGen] ), 
    # Show how long the program has been running.
    pyEval(r'"Burn-in stage: generation %d\n" % gen ', step = 10, end=burnin),
    # show elapsed time
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
  # may need to run several times to
  # get a valid population (terminate if any DSL has no disease allele) 
  fixedCount = 1
  while(True):
    # create a simulator
    pop =  population(subPop=popSizeFunc(0), ploidy=2,
      loci = loci, maxAllele = maxAle, lociDist = lociDist)
    # save DSL info, some operators will use it.
    pop.dvars().DSL = DSL
    pop.dvars().numLoci = numLoci
    pop.dvars().minAlleleFreq = minAlleleFreq
    pop.dvars().maxAlleleFreq = maxAlleleFreq
    # clear log file if it exists
    try:
      os.remove(logFile)
    except:
      pass
    # start simulation.
    simu = simulator( pop, 
      randomMating(
        newSubPopSizeFunc=popSizeFunc,  # demographic model
        numOffspringFunc=last_two),     # save last two generations
      rep=1)
    # evolve! If --dryrun is set, only show info
    simu.evolve( preOps = preOperators, ops = operators, 
      end=endGen, dryrun=dryrun )
    if dryrun:
      raise exceptions.SystemError("Stop since in dryrun mode.")
    # if succeed
    if simu.gen() == endGen + 1:  # if not fixed
      pop = simu.population(0)
      # we want to save info on how this population is generated.
      # This is not required but is a good practise
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
      pop.dvars().mutationRate = mu
      pop.dvars().mutationModel = "symmetric stepwise"
      pop.dvars().migrationRate = mi
      pop.dvars().migrationModel = "circular stepping stone"
      pop.dvars().recombinationRate = rec
      # return a copy (and then simulator will be destroyed)
      #return simu.getPopulation(0)
      # simu will be destroyed after return, however, simu.getPopulation
      # will make a copy of population(0) and double the use of memory
      # at this time.... this will cause problem when population size
      # is large. 
      # I am using a trick here, namely, create a small population and 
      # swap out the population in simulator. This is undocumented
      # and should be avoid whenever possible
      tmpPop = population(1)
      tmpPop.swap(simu.population(0))
      return tmpPop      
    else:
      print "Population restarted at gen ", simu.gen()
      print "Overall fixed population ", fixedCount
      print "Allelefreq ", ( "%.3f " * numDSL + "\n\n") % \
        tuple([1-simu.dvars(0).alleleFreq[x][1] for x in DSL]) 
      fixedCount += 1

# plot the LD values for the sample.
def plotLD(pop, epsFile, jpgFile):
  ''' plot LD values in R and convert to jpg if possible '''
  # return max LD
  res = {}
  # dist: distance (location) of marker
  # ldprime: D' value
  dist = []
  ldprime = [] # D'
  ldvalue = [] # D
  for ld in pop.dvars().ctrDSLLD:
    if ld[1] == pop.dvars().ctrChromDSL:
      dist.append(pop.locusDist(ld[0]))
    else:
      dist.append(pop.locusDist(ld[1]))
    ldprime.append(pop.dvars().LD_prime[ld[0]][ld[1]])
    ldvalue.append(pop.dvars().LD[ld[0]][ld[1]])
  res['DpDSL'] = max(ldprime)
  res['DDSL'] = max(ldvalue)
  if hasRPy:
    r.postscript(file=epsFile)
    r.par(mfrow=[2,1])
    r.plot( dist, ldprime, main="D' between DSL and other markers on chrom %d" % (pop.dvars().ctrChrom+1),
      xlab="marker location", ylab="D'", type='b', ylim=[0,1])
    r.abline( v = pop.locusDist(pop.dvars().ctrChromDSL), lty=3 )
    r.axis( 1, [pop.locusDist(pop.dvars().ctrChromDSL)], ['DSL'])
  dist = [] 
  ldprime = []  # D'
  ldvalue = []  # D
  for ld in pop.dvars().noDSLLD:
    if ld[1] == pop.chromBegin(pop.dvars().noDSLChrom) + numLoci/2:
      dist.append(pop.locusDist(ld[0]))
    else:
      dist.append(pop.locusDist(ld[1]))
    ldprime.append(pop.dvars().LD_prime[ld[0]][ld[1]])    
    ldvalue.append(pop.dvars().LD[ld[0]][ld[1]])    
  res['DpNon'] = max(ldprime)
  res['DNon'] = max(ldvalue)
  if hasRPy:
    r.plot( dist, ldprime, main="D' between marker %d and other markers on chrom %d" \
      % (numLoci/2+1, pop.dvars().noDSLChrom+1),
      xlab="marker location", ylab="D'", type='b', ylim=[0,1])    
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

# penetrance generator functions. They will return a penetrance function
# with given penetrance parameter
def recessive(pen):
  ''' recessive single-locus, heterogeneity multi-locus '''
  def func(geno):
    val = 1
    for i in range(len(geno)/2):
      if geno[i*2]+geno[i*2+1] == 4:
        val *= 1 - pen
    return 1-val
  return func
  
def additive(pen):
  ''' additive single-locus, heterogeneity multi-locus '''
  def func(geno):
    val = 1
    for i in range(len(geno)/2):
      val *= 1 - (geno[i*2]+geno[i*2+1]-2)*pen/2.
    return 1-val
  return func

# if you need some specialized penetrance function, modify this
# function here.
def custom(pen):
  def func(geno):
    raise "customized penetrance function has not been defined."

def drawSamples(pop, penFun, numSample):
  ''' get samples of different type using a penetrance function '''
  # first, apply peneFunction
  PyPenetrance(pop, loci=pop.dvars().DSL, func=penFun)
  # all types of samples
  allSample = []
  for ns in range(numSample):
    print "Generating sample ", ns+1, ' of ', numSample
    allSample.append([])
    # 1. population based case control
    # get number of affected
    Stat(pop, numOfAffected=True)
    print "Number of affected individuals: ", pop.dvars().numOfAffected
    print "Number of unaffected individuals: ", pop.dvars().numOfUnaffected
    nCase = min(pop.dvars().numOfAffected , N/2)
    nControl = min(pop.dvars().numOfUnaffected, N/2)
    try:
      # if N=800, 400 case and 400 controls
      s = CaseControlSample(pop, nCase, nControl)
      # remove DSL
      s[0].removeLoci(remove=pop.dvars().DSL)
      allSample[ns].append(s[0])
    except Exception, err:
      print "Can not draw case control sample. "
      print type(err), err
      allSample[ns].append(None)
    #
    # 2. affected and unaffected sibpairs
    # this is difficult since simuPOP only has
    # methods to draw one kind of sibpairs and 
    try:
      # get number of affected/unaffected sibpairs
      # There may not be enough to be sampled
      AffectedSibpairSample(pop, countOnly=True)
      nAff = min(pop.dvars().numAffectedSibpairs, N/4)
      print "Number of (both) affected sibpairs: ", pop.dvars().numAffectedSibpairs
      AffectedSibpairSample(pop, chooseUnaffected=True, countOnly=True)
      print "Number of unaffected sibpairs: ", pop.dvars().numAffectedSibpairs
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
      affected[0].removeLoci(remove=pop.dvars().DSL)
      unaffected[0].removeLoci(remove=pop.dvars().DSL)
      allSample[ns].extend([affected[0], unaffected[0] ])
    except Exception, err:
      print type(err)
      print err
      print "Can not draw affected sibpars."
      allSample[ns].extend([None, None])
  return allSample

# apply the TDT method of GeneHunter
def TDT(DSL, cutoff, dataDir, data, epsFile, jpgFile):
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
    xlab="chromosome", ylab="-log10 p-value", type='l', axes=False, ylim=[0.01, 5])
  r.box()
  r.abline( v = [DSLafter[g]+DSLdist[g] for g in range(len(DSL))], lty=3)
  r.abline( h = cutoff )                       
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

# apply the Linkage method of GeneHunter
def Linkage(DSL, cutoff, dataDir, data, epsFile, jpgFile):
  ''' use Linkage method to analyze the results. Has to have rpy installed '''
  if not hasRPy: 
    return (0,[])
  # write a batch file and call gh
  allPvalue = []
  print "Applying Linkage (LOD) method to affected sibpairs "
  for ch in range(numChrom):
    inputfile = dataDir+data+ "_%d" % ch
    if not os.path.isfile(inputfile + ".ped"):
      print "Ped file ", inputfile+".ped does not exist. Can not apply TDT method."
      return (0,[])
    # batch file
    f=open("ghLOD.cmd","w")
    f.write("load markers " + inputfile + ".dat\n")
    f.write("single point on\n")
    f.write("scan pedigrees " + inputfile + ".ped\n")
    f.write("photo tmp.txt\n")
    f.write("total stat\n")
    f.write("q\n")
    f.close()
    # run gene hunter
    os.system(geneHunter + " < ghLOD.cmd > res.txt ")
    # get the result
    # ( I was using              
    # | grep '^loc' | tr -d 'a-zA-Z+-' > " + outputfile)
    # but this is not portable
    #
    # use re package
    # get only loc number and p-value
    # position (locxx) -- Lodscore -- NPL score -- p-value -- information
    scan = re.compile('loc(\d+)\s+[^\s]+\s+[^\s]+\s+([^\s]+)\s*.*')
    minPvalue = [1]*(numLoci-1)
    try:
      res = open("tmp.txt")  # read result
      for l in res.readlines():
        try:
          # get minimal p-value for all alleles at each locus
          (loc,pvalue) = scan.match(l).groups()
          if minPvalue[int(loc)-1] > float(pvalue):
            minPvalue[int(loc)-1] = float(pvalue)
        except:
          # does not match
          #print l
          #print err
          continue
      res.close()
    except Exception, err:
      print type(err), err
      print "Can not open result file tmp.txt. LOD failed"
      return (0,[])
    # did not get anything
    if minPvalue == [1]*(numLoci-1):
      print "Do not get any p-value, something is wrong"
    # collect -log10 pvalue
    allPvalue.extend([-math.log10(max(x,1e-6)) for x in minPvalue])
    # There will be numLoci-1 numbers, pad to the correct length
    allPvalue.append(0)
    try:
      os.unlink('res.txt')
      os.unlink('tmp.txt')
      os.unlink('ghLOD.cmd')
    except:
      pass
  # now, we need to see how Linkage works with a set of p-values around DSL
  # DSL is global
  res = []
  for d in DSL: 
    res.append( max(allPvalue[(d-2):(d+2)]))
  # use R to draw a picture
  r.postscript(file=epsFile)
  r.plot(allPvalue, main="-log10(P-value) for each marker (LOD)", ylim=[0.01,5], 
    xlab="chromosome", ylab="-log10 p-value", type='l', axes=False)
  r.box()
  r.abline( v = [DSLafter[g]+DSLdist[g] for g in range(len(DSL))], lty=3)
  r.abline( h = cutoff )                       
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
# a more friendly version of mkdir
def _mkdir(d):
  try:
    if not os.path.isdir(d):
      os.makedirs(d)
  except:
    print "Can not make output directory ", d
    sys.exit(1)

def popStat(pop, p):
  ' calculate population statistics '
  # K -- populaiton prevalance
  print "Calculating population statistics "
  Stat(pop, numOfAffected=True)
  result = {}
  result[p+'_K'] = pop.dvars().numOfAffected * 1.0 / pop.popSize()
  # P11 = [ ] = proportion of 11 | affected, 
  # P12 = [ ] = proportion of 12 | affected
  DSL = pop.dvars().DSL
  P11 = [0.]*len(DSL)
  P12 = [0.]*len(DSL)
  P22 = [0.]*len(DSL)
  for ind in range(pop.popSize()):
    if pop.individual(ind).affected():
      for x in range(len(DSL)):
        s1 = pop.individual(ind).allele(DSL[x], 0)
        s2 = pop.individual(ind).allele(DSL[x], 1)
        if s1 == 1 and s2 == 1:
          P11[x] += 1
        elif s1 == 2 and s2 == 2:
          P22[x] += 1
        else:
          P12[x] += 1
  N = pop.dvars().numOfAffected
  result[p+'_P11'] = [ x/N for x in P11 ]
  result[p+'_P12'] = [ x/N for x in P12 ]
  result[p+'_P22'] = [ x/N for x in P22 ]
  # Ks = Pr(Xs=1 | Xp=1 ) = Pr(Xs=1, Xp=1) | P(Xp=1)
  Ks = 0.
  for ind in range(pop.popSize()/2):
    s1 = pop.individual(ind*2).affected()
    s2 = pop.individual(ind*2+1).affected()
    if s1 and s2:
      Ks += 1
  result[p+'_Ks'] = 2*Ks / pop.dvars().numOfAffected
  # Lambda = Ks/K
  result[p+'_Ls'] = result[p+'_Ks'] / result[p+'_K']
  return result
   
  
def processOnePopulation(dataDir, numChrom, numLoci, markerType,
    DSLafter, DSLdist, initSize, meanInitAllele, burnin, introGen, minAlleleFreq,
    maxAlleleFreq, fitness, mlSelModel, numSubPop, finalSize, noMigrGen,
    mixingGen, popSizeFunc, migrModel, mu, mi, rec,  peneFunc, penePara, N, 
    numSample, dryrun, popIdx):
  '''
     this function organize all previous functions
     and
     1. generate (or load) a population
     2. apply different kinds of penetrance functions
     3. draw sample
     4. save samples
     5. apply TDT and/or Linkage method
     6. return a text summary and a result dictionary
  '''
  # get population
  result = {'id':popIdx, 'mu':mu, 'mi':mi, 'rec':rec}
  logFile = dataDir + "/pop_" + str(popIdx) + ".log"
  popFile = dataDir + "/pop_" + str(popIdx) + ".bin"
  genDataset = True
  if os.path.isfile(popFile) and (not overwrite):
    print "Loading a pre-existing file ", popFile
    pop = LoadPopulation(popFile)
    # check if the population is using the same parameters as requested
    if abs(pop.dvars().mutationRate - mu) + abs(pop.dvars().migrationRate - mi) \
       + abs( pop.dvars().recombinationRate - rec) < 1e-7:
      genDataset = False
      print "Dataset already exists, load it directly. (use --overwrite True to change this behavior)"
    else:
      print "A dataset is present but was generated with different parameters. Re-generating."
  if genDataset:
    print "Generating dataset ", str(popIdx)
    pop = simuComplexDisease( numChrom, numLoci, markerType, DSLafter, DSLdist,
      initSize, meanInitAllele, burnin, introGen, minAlleleFreq,
      maxAlleleFreq, fitness, mlSelModel, numSubPop, finalSize, noMigrGen,
      mixingGen, popSizeFunc, migrModel, mu, mi, rec, dryrun, logFile)
    # note that this file does not have affectedness info.
    SavePopulation(pop, popFile)
  # log file is ...
  summary = '''\
    <h2><a name="pop_%d">Dataset %d</a></h2>
    <h3>Log file (LD and other statistics): 
      <a href="pop_%d/pop_%d.log">pop_%d.log</a></h3>
    ''' % (popIdx, popIdx, popIdx, popIdx, popIdx)
  #
  # actually write the log file in the summary page
  try:
    lf = open(logFile)
    summary += "<pre>\n"+lf.read()+"\n</pre>\n"
    lf.close()
  except:
    print "Can not open log file, ignoring. "
  #
  # save Fst, Het in res
  result['Fst'] = pop.dvars().AvgFst
  result['AvgHet'] = pop.dvars().AvgHetero
  result['alleleFreq'] = [1- pop.dvars().alleleFreq[i][1] for i in pop.dvars().DSL]
  #
  # plot LD, res = 0, fail, 1: eps file, 2: converted to jpg
  epsFile = dataDir + "/LD_" + str(popIdx) + ".eps"
  jpgFile = dataDir + "/LD_" + str(popIdx) + ".jpg"
  
  # ldres has max D' on a chrom with DSL, and a chrom without DSL
  (suc,ldres) = plotLD(pop, epsFile, jpgFile)
  if suc > 0 : # eps file successfully generated
    summary += """<p>D' measures on two chromosomes with/without DSL at the last gen: 
    <a href="pop_%d/LD_%d.eps">LD.eps</a></p>\n""" % (popIdx, popIdx)
  if suc > 1 : # jpg file is generated
    summary += '''<img src="pop_%d/LD_%d.jpg"
      width=800, height=600>'''  % (popIdx, popIdx)
  result['DpDSL'] = ldres['DpDSL']
  result['DpNon'] = ldres['DpNon']
  result['DDSL'] = ldres['DDSL']
  result['DNon'] = ldres['DNon']
  #
  # apply penetrance and get numSample for each sample
  summary += "<h3>Samples using different penetrance function</h3>\n"
  for p in range(len(peneFunc)):
    if peneFunc[p] == 'recessive':
      print "Using recessive penetrance function"
      summary += "<h4>Recessive single-locus, heterogeneity multi-locus</h4>\n"
      s = drawSamples(pop, recessive( penePara[p]), numSample)
    elif peneFunc[p] == 'additive':
      print "Using additive penetrance function"
      summary += "<h4>Additive single-locus, heterogeneity multi-locus</h4>\n"
      s = drawSamples(pop, additive(penePara[p]), numSample)
    elif peneFunc[p] == 'custom':
      print "Using customized penetrance function"
      summary += "<h4>Customized penetrance function</h4>\n"
      s = drawSamples(pop, customPene(penePara[p]), numSample)
    # calculate population statistics like prevalence
    result.update( popStat(pop, peneFunc[p]) )
    # for each sample
    for sn in range(numSample):
      print "Processing sample %s%d" % ( peneFunc[p], sn)
      # save these samples
      penDir = dataDir + "/" + peneFunc[p] + str(sn)
      relDir = 'pop_%d/%s%d/' % (popIdx, peneFunc[p], sn)
      _mkdir(penDir)
      summary += "<p>Case-control, affected and unaffected sibpairs saved in different formats. Sample sizes are "
      if s[sn][0] != None:
        summary += str(s[sn][0].popSize()) + " (case-control), "
      if s[sn][1] != None:
        summary += str(s[sn][1].popSize()) + " (affected sibs), "
      if s[sn][2] != None:
        summary += str(s[sn][2].popSize()) + " (unaffected sibs) "
      summary += '<ul>'
      # write samples to different format
      if 'simuPOP' in saveFormat:
        print "Write to simuPOP binary format"
        summary += '<li>simuPOP binary format:'
        if s[sn][0] != None: # has case control
          binFile = penDir + "/caseControl.bin"
          s[sn][0].savePopulation(binFile)
          summary += '<a href="%scaseControl.bin"> caseControl.bin</a>, ' % relDir
        if s[sn][1] != None: # has affected sibpairs
          binFile = penDir + "/affectedSibpairs.bin"
          s[sn][1].savePopulation(binFile)
          summary += '<a href="%saffectedSibpairs.bin"> affectedSibpairs.bin</a>, ' % relDir
        if s[sn][2] != None:
          binFile = penDir + "/unaffectedSibpairs.bin"
          s[sn][2].savePopulation(binFile)
          summary += '<a href="%sunaffectedSibpairs.bin"> unaffectedSibpirs.bin</a>, ' % relDir
        summary += '</li>\n'
      if 'Linkage' in saveFormat:
        print "Write to linkage format"
        summary += '<li>Linkage format by chromosome:'
        linDir = penDir + "/Linkage/"
        _mkdir(linDir)
        # here we are facing a problem of using which allele frequency for the sample
        # In reality, if it is difficult to estimate population allele frequency,
        # sample frequency has to be used. Otherwise, population frequency should 
        # be used whenever possible. Here, we use population allele frequency, with only
        # one problem in that we need to remove frequencies at DSL (since samples do not
        # have DSL).
        af = []
        Stat(pop, alleleFreq=range(pop.totNumLoci()))
        for x in range( pop.totNumLoci() ):
          if x not in pop.dvars().DSL:
            af.append( pop.dvars().alleleFreq[x] )
        if s[sn][1] != None: # has case control
          for ch in range(0, pop.numChrom() ):
            SaveLinkage(pop=s[sn][1], popType='sibpair', output = linDir+"/Aff_%d" % ch,
              chrom=ch, recombination=pop.dvars().recombinationRate,
              alleleFreq=af, daf=0.1)        
          summary +=  '<a href="%sLinkage">affected</a>, ' % relDir
        if s[sn][2] != None:
          for ch in range(0,pop.numChrom() ):
            SaveLinkage(pop=s[sn][2], popType='sibpair', output = linDir+"/Unaff_%d" % ch,
              chrom=ch, recombination=pop.dvars().recombinationRate,                            
              alleleFreq=af, daf=0.1)        
          summary += '<a href="%sLinkage">unaffected</a>' % relDir
        summary += '</li>\n'
      summary += '</ul>\n'
      # if there is a valid gene hunter program, run it
      (suc,res) = TDT(pop.dvars().DSL, -math.log10(0.05/pop.totNumLoci()), 
        penDir, "/Linkage/Aff", penDir + "/TDT.eps", penDir + "/TDT.jpg")
      #  if suc > 0 : # eps file succe
      if suc > 0 : # eps file successfully generated
        summary += """<h4>TDT analysis for affected sibpair data:  <a href="%s/TDT.eps">TDT.eps</a>""" % relDir
      if suc > 1 : # jpg file is generated
        summary += '''<p><img src="%s/TDT.jpg" width=800, height=600></p>''' % relDir
      # keep some numbers depending on the penetrance model
      result['TDT_%s_%d' % (peneFunc[p], sn)] = res
      # then the Linkage method
      (suc,res) = Linkage(pop.dvars().DSL, -math.log10(0.05/pop.totNumLoci()), 
        penDir, "/Linkage/Aff", penDir + "/LOD.eps", penDir + "/LOD.jpg")
      #  if suc > 0 : # eps file succe
      if suc > 0 : # eps file successfully generated
        summary += """<h4>LOD analysis for affected sibpair data:  <a href="%s/LOD.eps">LOD.eps</a>""" % relDir
      if suc > 1 : # jpg file is generated
        summary += '''<p><img src="%s/LOD.jpg" width=800, height=600></p>''' % relDir
      # keep some numbers depending on the penetrance model
      result['LOD_%s_%d' % (peneFunc[p], sn)] = res
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
  <TITLE>Summary of simulations</TITLE>
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
  rate, recombination rate, Fst, average heterozygosity, highest D'/D between a 
  DSL and all surrounding markers (not necessarily its cloest marker), highest
  D'/D between a marker with its surrounding markers on a chromsome without
  DSL, allele frequency at DSL, -log10 p-values (TDT method and Linkage method, + for
  exceeds and - for less than cutoff value -log10(pvalue/total number of loci) ) at
  all relevant DSL. Other statistics include K (population prevalence), Ks (sibling 
  recurrance risk), Ls (lambda_s, sibling recurrance ratio), P11 (P(NN | affected)),
  P12 (P(NS | affected)), P13 (P(SS|affected))</p>
  <table border="1">
  <tr><th>id </th>
  <th>mu</th>  <th>mi</th>   <th>rec</th>
  <th>Fst</th> <th>Het</th> 
  <th>D'(dsl)</th> <th>D (dsl)</th>
  <th>D'(non)</th> <th>D (non)</th>
  ''')
  for i in range(len(allParam[3])):
    summary.write('<th>allele Frq%d</th>'%(i+1))
  #
  # has TDT and some penetrance function
  if len(allParam[-9]) > 0:
    summary.write('<th>K</th><th>Ks</th><th>Ks/K</th>')
    summary.write('<th>P11</th><th>P12</th><th>P22</th>')
    for p in allParam[-9]:
      summary.write('<th>%s:TDT</th><th>%s:LOD</th>'%(p,p))
  summary.write('</tr>')
  #
  # end of headers, now result
  for res in results:
    summary.write('''<tr><td><a href="#pop_%d">%d</a></td> ''' \
      % (res['id'], res['id']))    
    summary.write('<td>%.5g</td>' % res['mu'])
    summary.write('<td>%.5g</td>' % res['mi'])
    summary.write('<td>%.5g</td>' % res['rec'])
    summary.write('<td>%.2g</td>' % res['Fst'])
    summary.write('<td>%.2g</td>' % res['AvgHet'])
    summary.write('<td>%.2g</td>' % res['DpDSL'])
    summary.write('<td>%.2g</td>' % res['DDSL'])
    summary.write('<td>%.2g</td>' % res['DpNon'])
    summary.write('<td>%.2g</td>' % res['DNon'])
    for i in range(len(allParam[3])):  # len(DSLafter)
      summary.write('<td>%.3f</td>'% res['alleleFreq'][i] )
    # for each penetrance function
    if len(allParam[-9]) > 0:
      for p in allParam[-9]: # penetrance function
        summary.write('<td>%.3g</td>' % res[p+'_K'])
        summary.write('<td>%.3g</td>' % res[p+'_Ks'])
        summary.write('<td>%.3g</td>' % res[p+'_Ls'])
        summary.write('<td>' + ','.join( ['%.2g'%x for x in res[p+'_P11'] ]) + '</td>')
        summary.write('<td>' + ','.join( ['%.2g'%x for x in res[p+'_P12'] ]) + '</td>')
        summary.write('<td>' + ','.join( ['%.2g'%x for x in res[p+'_P22'] ]) + '</td>')
        for met in ['TDT', 'LOD']:
          for num in range(int(allParam[-6])): # samples
            plusMinus = ''
            for pvalue in res[met+'_'+p+'_'+str(num)]:
              if pvalue > -math.log10(0.05/(allParam[0]*allParam[1])):
                plusMinus += '+'
              else:
                plusMinus += '-'
            summary.write('<td>'+plusMinus+'</td>')
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
  (numChrom, numLoci, markerType, DSLafter, DSLdistTmp, initSize, meanInitAllele, 
    burnin, introGen, minAlleleFreqTmp, maxAlleleFreqTmp, fitnessTmp, 
    mlSelModelTmp, numSubPop, finalSize, noMigrGen,
    mixingGen, growth, migrModel, migrRate, mutaRate, recRate,
    numRep, saveFormat, peneFunc, peneParaTmp, N, numSample, outputDir,
    overwrite, geneHunter, dryrun, saveConfig) = allParam
  #
  # this may not happen at all but we do need to be careful
  if initSize < len(DSLafter):
    raise exceptions.ValueError("Initial population size is too small. (Less than number of DSL)")
  #
  # expand .5 => .5,.5,.5 etc
  if len(DSLdistTmp) == 1:
    DSLdist = DSLdistTmp * len(DSLafter)
  else:
    DSLdist = DSLdistTmp
  if len( DSLafter ) != len(DSLdist):
    raise exceptions.ValueError("Please specify DSL distance for each DSL.")
  #  
  # handle minAlleleFreq and maxAlleleFreq
  if len(minAlleleFreqTmp) == 1:
    minAlleleFreq = minAlleleFreqTmp * len(DSLafter)
  elif len(minAlleleFreqTmp) == len(DSLafter):
    minAlleleFreq = minAlleleFreqTmp
  else:
    raise exceptions.ValueError("min allele frequency should be either a number or a list\n" +
      " of the same length as DSLafter")
  if len(maxAlleleFreqTmp) == 1:
    maxAlleleFreq = maxAlleleFreqTmp * len(DSLafter)
  elif len(maxAlleleFreqTmp) == len(DSLafter):
    maxAlleleFreq = maxAlleleFreqTmp
  else:
    raise exceptions.ValueError("max allele frequency should be either a number or a list\n" +
      " of the same length as DSLafter")
  #  
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
      
  # check penePara
  if len(peneParaTmp) == 1:
    penePara = peneParaTmp*len(peneFunc)
  else:
    penePara = peneParaTmp
  if len(peneFunc) != len(penePara):
    raise exceptions.ValueError("Please give penetrance parameter to each chosen penetrance function")
  #
  # everything is ready
  popIdx = 1
  summary = ''
  results = []
  for mu in mutaRate:
    for mi in migrRate:
      for rec in recRate:
        for rep in range(numRep):
          # outputDir should already exist
          dataDir = outputDir + "/pop_" + str(popIdx)
          _mkdir(dataDir)      
          (text, result) =  processOnePopulation(dataDir,
            numChrom, numLoci, markerType, 
            DSLafter, DSLdist, initSize, meanInitAllele, burnin, introGen, 
            minAlleleFreq, maxAlleleFreq, fitness, mlSelModel, numSubPop, 
            finalSize, noMigrGen, mixingGen, popSizeFunc, migrModel, 
            mu, mi, rec, peneFunc, penePara, N, numSample, dryrun, popIdx)
          summary += text
          results.append( result)
          popIdx += 1
  writeReport(summary, allParam, results)
  print "Done!"
