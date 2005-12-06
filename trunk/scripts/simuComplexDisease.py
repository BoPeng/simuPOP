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

The result of the simulation is a large multi-generation population. To analyze 
the population, you will typically need to 

  1. Apply a penetrance function to the population and
     determine the affectedness for each individual

  2. Draw Population and/or pedigree based samples and save in popular 
     formats so that the samples can be analyzed by other software like
     genehunter.
     
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


Statistics Monitored
====================

A number of statistics will be measured and saved. They are:

  1. Fst before and after mixing
  2. Observed heterogeneity before and after mixing
  3. LD (D prime) between markers and between a DSL and surrounding markers
     at the end of each stage.
  4. Disease allele frequency at the end of each stage.

"""

import simuOpt, simuUtil
import os, sys, exceptions, types

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
   'configName': 'Number of chromosomes',
   'prompt': 'Number of chromosomes (20):  ',
   'description': 'Number of chromosomes.',
   'allowedTypes': [types.IntType],
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'numLoci=',
   'default': 20,
   'configName': 'Number of loci on each chromosome',
   'prompt': 'Number of loci on each chromosome (20):  ',
   'description': '''Number of loci on each chromosome, current there 
       only equal number of markers on each chromosome is supported.''',
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
   'configName': 'DSL after marker (0-indexed)',
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
   'configName': 'DSL location between markers',
   'prompt': 'Enter the position of each DSL between two markers ([.01,.5,.75]):  ',
   'description': '''A list of loci location between two markers.
        Since all disease loci will be *between* equal spaced markers,
        the location should be between 0 and 1. A single value is acceptable
        as the location of all DSL.''',
   'allowedTypes': [types.TupleType, types.ListType],
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },  
  {'longarg': 'initSize=',
   'default': 10000,
   'configName': 'Initial population size',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Initial Population size (10000):  ',
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
   'default': 500,
   'configName': 'Length of burn-in stage',
   'allowedTypes': [types.IntType],
   'prompt': 'Length of burn in stage (500):  ',
   'description': 'Number of generations of the burn in stage.',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'introGen=',
   'default': 50,
   'configName': 'Length of disease-intro stage',
   'allowedTypes': [types.IntType],
   'prompt': 'Length of disease-introduction stage (50):  ',
   'description': '''Number of generations to introduce the disease
        Since the disease will be under positive selection during this
        stage, large introGen will lead to over-common diseases. Note 
        that the larger the initial population size is, the longer 
        introGen is required.''',
   'validate':  simuOpt.valueGT(0),
  },
  {'longarg': 'minAlleleFreq=',
   'default': [0.045],
   'configName': 'Minimal allele frequency',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Minimum allele frequency required for each DSL (0.045):  ',
   'description': '''Mininal allele frequencies required for each DSL.
        If a number is given, it is assumed to be the lower bound for all
        DSL.''',
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },
  {'longarg': 'maxAlleleFreq=',
   'default': [0.055],
   'configName': 'Maximum allele frequency',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Maximum allele frequency required for each DSL (0.055):  ',
   'description': '''Maximum allele frequencies required for each DSL.
        If a number is given, it is assumed to be the upper bound for all
        DSL.''',
   'validate':  simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },
  {'longarg': 'fitness=',
   'default': [1,1,1],
   'configName': 'Fitness of genotype AA,Aa,aa',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': ':  Fitness of genotype AA,Aa,aa ([1,1,1]): ',
   'description': '''Fitness of genotype AA,Aa,aa for every DSL after disease introduction
        stage. This can be an array (default) of fitness values for genotype AA,Aa,aa ( A
        is wild type) or an array of such array for each DSL.''',
   'validate':  simuOpt.valueGE(0.),
  },
  {'longarg': 'selMultiLocusModel=',
  'default': 'none',
  'configName': 'Multi-locus selection model',
  'prompt': 'Multi-locus selection model (additive): ',
  'description': '''Model of overall fitness value given fitness values for each DSL.
        fitness values are Prod(f_i) for multiplicative model and
        1-Sum(1-f_i) for additive model. ''',
  'allowedTypes': [types.StringType],
  'chooseOneOf': [ 'additive', 'multiplicative', 'none']
  }, 
  {'longarg': 'numSubPop=',
   'default': 10,
   'configName': 'Number of subpopulations to split',
   'allowedTypes': [types.IntType],
   'prompt': 'Number of subpopulations to split (10):  ',
   'description': 'Number of subpopulations to be split into after burnin stage.',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'finalSize=',
   'default': 200000,
   'configName': 'Final population size',
   'prompt': 'Final population size (sum of all subpopulations) (200000):  ',
   'allowedTypes': [types.IntType, types.LongType],
   'description': 'Final population size after population expansion.',
   'validate':  simuOpt.valueGT(0)
  }, 
  {'longarg': 'noMigrGen=',
   'default': 400,
   'configName': 'Length of split-and-grow stage',
   'prompt': 'Length of split-and-grow stage  (400):  ',
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
   'configName': 'Population growth model',
   'prompt': 'Population growth style, linear or exponential. (exponential):  ',
   'description': '''How population is grown from initSize to finalSize.
        Choose between linear and exponential''',
   'chooseOneOf': ['exponential', 'linear'],
  },
  {'longarg': 'migrModel=',
   'default': 'stepping stone',
   'configName': 'Migration model',
   'prompt': 'Migration model. (stepping stone):  ',
   'allowedTypes': [types.StringType],
   'description': '''Migration model. Choose between stepping stone (circular),
        island and none. ''',
   'validate':  simuOpt.valueOneOf(['island', 'stepping stone', 'none']),
   'chooseOneOf': ['stepping stone', 'island', 'none']
  }, 
  {'longarg': 'migrRate=',
   'default': 0.0001,
   'configName': 'Migration rate',
   'prompt': 'Migration rate during mixing stage. (0.0001):  ',
   'description': '''Migration rate during mixing stage. 
        A circular stepping stone migration model will be used. ''',
   'allowedTypes': [types.IntType, types.FloatType],
   'validate':  simuOpt.valueBetween(0,1)
  },
  {'longarg': 'mutaRate=',
   'default': 0.0001,
   'configName': 'Mutation rate',
   'prompt': 'Mutation rate at non-DSL markers. (0.0001):  ',
   'allowedTypes': [types.IntType, types.FloatType],
   'description': '''Microsatellite markers are mutated using  
        symmetric stepwise mutation wile SNP markers are mutaed
        using a 2-allele model (kam) DSL are not mutated unless in disease
        introduction stage.''',
   'validate':  simuOpt.valueBetween(0,1)
  },
  {'longarg': 'recRate=',
   'default': 0.0005,
   'configName': 'Recombination rate',
   'allowedTypes': [types.IntType, types.FloatType],
   'prompt': 'Recombination rate between adjacent markers. (0.0005): ',
   'description': '''Recombination rate between adjacent markers. Note 
         that this value will be different between marker and DSL depending
         on the location of DSL. ''',
   'validate':  simuOpt.valueBetween(0,1)
  },
  {'longarg': 'dryrun',
   'default': False,
   'allowedTypes': [types.IntType],
   'validate':  simuOpt.valueOneOf([True, False]),
   'description':  'Only display how simulation will perform.'
   # do not save to config, do not prompt, so this appeared to be an undocumented option.
  },
  {'longarg': 'simuName=',
   'default': 'simu',
   'allowedTypes': [types.StringType],
   'configName': 'Simulation name',
   'prompt': 'Name of the simulation (simu):  ',
   'description': '''Name of simulation, files saved will be 
          name + '.log': statistics output
          name + '.cfg': configuration
          name + '.bin': saved popuation''',
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
    '''  This program simulates the evolution of a complex common disease, subject 
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
    simuOpt.saveConfig(options, allParam[-2]+'.cfg', allParam)
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

def outputStatistics(pop, args): 
  ''' This function will be working with a pyOperator
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
  toCtrDist = [abs(pop.locusPos(x)-numLoci/2) for x in DSL]
  ctrChromDSL = DSL[ [x==min(toCtrDist) for x in toCtrDist].index(True)]
  ctrChrom = pop.chromLocusPair(ctrChromDSL)[0]
  # 2. first chromosome without DSL
  try:
    noDSLChrom = [pop.numLoci(x)==numLoci for x in range(pop.numChrom())].index(True)
  except:  # there is no such chromosome!
    noDSLChrom = -1
  #
  # loci pairs on the selected chromosomes when calculating LD
  i = pop.chromBegin( ctrChrom)
  ctrDSLLD = [ [x, ctrChromDSL] for x in range(i, ctrChromDSL)] + \
    [ [ctrChromDSL, x] for x in range(ctrChromDSL+1, i+numLoci+1)]
  if noDSLChrom > -1:
    i = pop.chromBegin( noDSLChrom)
    noDSLLD = [ [i+x, i+numLoci/2] for x in range(numLoci/2)] + \
      [ [i + numLoci/2, i+x] for x in range(numLoci/2+1, numLoci)]
  else:
    noDSLLD = []
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
    if noDSLChrom > -1 :
      print >> output, "\n\nD between a center marker %d (chrom %d) and surrounding markers at gen %d" \
        % (pop.chromBegin(noDSLChrom)+numLoci/2, noDSLChrom, gen)
      for ld in noDSLLD:
        print >> output, '%.4f ' % pop.dvars().LD[ld[0]][ld[1]],
    print >> output, "\n\nD' between DSL %d (chrom %d) and surrounding markers at gen %d" \
      % (ctrChromDSL, ctrChrom, gen)
    for ld in ctrDSLLD:
      print >> output, '%.4f ' % pop.dvars().LD_prime[ld[0]][ld[1]],
    if noDSLChrom > -1:
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

def dynaAdvSelector(pop):
  ''' This selector will apply advantage/purifying selection to DSl according
    to allele frequency at each DSL. minAlleleFreq and maxAlleleFreq
    are stored in pop.dvars().
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
      print ' ',
      sel.append( maSelector(locus=DSL[i], wildtype=[1], fitness=[1,1.02,1.04]) )
    # apply multi-locus selector, note that this operator will only
    # set a variable fitness in pop, actual selection happens during mating.
    if len(sel ) > 0:  # need adjustment (needed if 'else' part is empty)
      MlSelect( pop, sel, mode=SEL_Multiplicative)
  return True

# simulate function, using a single value of mutation, migration,
# recombination rate
def simuComplexDisease( numChrom, numLoci, markerType, DSLafter, DSLdistTmp,
    initSize, meanInitAllele, burnin, introGen, minAlleleFreqTmp,
    maxAlleleFreqTmp, fitnessTmp, mlSelModelTmp, numSubPop, finalSize, noMigrGen,
    mixingGen, growth, migrModel, mu, mi, rec, dryrun, filename):
  ''' run a simulation of complex disease with given parameters. 
  '''
  ### CHECK PARAMETERS ###
  if initSize < len(DSLafter):
    raise exceptions.ValueError("Initial population size is too small. (Less than number of DSL)")
  # expand .5 -> [0.5, 0.5...]
  if len(DSLdistTmp) == 1:
    DSLdist = DSLdistTmp * len(DSLafter)
  else:
    DSLdist = DSLdistTmp
  if len( DSLafter ) != len(DSLdist):
    raise exceptions.ValueError("Please specify DSL distance for each DSL.")
  # handle minAlleleFreq and maxAlleleFreq
  # expand 0.5 -> [0.5,0.5,...]
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
  if markerType == 'microsatellite':
    maxAle = 99                   # max allele
  else:
    maxAle = 2                    # SNP 
  #### translate numChrom, numLoci, DSLafterLoci to
  #### loci, lociPos, DSL, nonDSL in the usual index
  #### (count markers along with DSL)
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
  i = 0  # current absolute locus 
  j = 0  # current DSL
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
  # non DSL loci, for convenience
  nonDSL = range(0, reduce(operator.add, loci))
  for loc in DSL:
    nonDSL.remove(loc)
  #### Change generation duration to 'generation-points'
  split  = burnin + introGen
  mixing = split  + noMigrGen
  endGen = mixing + mixingGen  
  #### initialization and mutation
  if maxAle > 2:  # Not SNP
    preOperators = [
      # initialize all loci with 5 haplotypes
      initByValue(value=[[x]*reduce(operator.add, loci) \
        for x in range(meanInitAllele-2, meanInitAllele+3)],
        proportions=[.2]*5), 
      # and then init DSL with all wild type alleles
      initByValue([1]*len(DSL), atLoci=DSL)
    ]
    # symmetric mutation model for microsatellite
    mutator = smmMutator(rate=mu, maxAllele=maxAle, atLoci=nonDSL)
  else: # SNP
    preOperators = [
      # initialize all loci with two haplotypes (111,222)
      initByValue(value=[[x]*reduce(operator.add, loci) for x in range(1,3)],
        proportions=[.5]*2), 
      # and then init DSL with all wild type alleles
      initByValue([1]*len(DSL), atLoci=DSL)
    ]      
    # k-allele model for mutation of SNP
    mutator = kamMutator(rate=mu, maxAllele=2, atLoci=nonDSL)
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
    # output progress 
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
  # if not, terminate the simulation.
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
  # start selection
  try:
    mlSelModel = {'additive':SEL_Additive, 
      'multiplicative':SEL_Multiplicative,
      'none':SEL_None}[mlSelModelTmp]
  except:
    raise exceptions.ValueError("Wrong multi-locus seleciton model. " + mlSelModelTmp)
  
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
      param = (burnin, split, mixing, endGen, filename+'.log'),
      at = [burnin, split, mixing, endGen] ), 
    # Show how long the program has been running.
    pyEval(r'"Burn-in stage: generation %d\n" % gen ', step = 10, end=burnin),
    # show elapsed time
    ticToc(at=[burnin, split, mixing, endGen]),
    ]
  )
  # with all operators, we can set up a simulator
  # produce two offsprings only at the last generations.
  def last_two(gen):
    if gen >= endGen -2:
      return 2
    else:
      return 1
  # population growth model
  if growth == 'linear':
    popSizeFunc = LinearExpansion(initSize, finalSize, endGen,
      burnin, split, numSubPop)
  elif growth == 'exponential':
    popSizeFunc = ExponentialExpansion(initSize, finalSize, endGen,
      burnin, split, numSubPop)
  else:
    raise exceptions.ValueError("Growth model can be one of linear and exponential. Given " + growth)
  # may need to run several times to
  # get a valid population (terminate if desired allele frequencyany can not be reached.)
  fixedCount = 1
  while(True):
    # create a simulator
    pop =  population(subPop=popSizeFunc(0), ploidy=2,
      loci = loci, maxAllele = maxAle, lociPos = lociPos)
    # save DSL info, some operators will use it.
    pop.dvars().DSL = DSL
    pop.dvars().numLoci = numLoci
    pop.dvars().minAlleleFreq = minAlleleFreq
    pop.dvars().maxAlleleFreq = maxAlleleFreq
    # clear log file if it exists
    try:
      os.remove(filename+'.log')
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
      print "Saving population to " + filename + ".bin\n"
      simu.population(0).savePopulation(filename+'.bin')
      return True
    else:
      print "Population restarted at gen ", simu.gen()
      print "Overall fixed population ", fixedCount
      print "Allelefreq ", ( "%.3f " * numDSL + "\n\n") % \
        tuple([1-simu.dvars(0).alleleFreq[x][1] for x in DSL]) 
      fixedCount += 1

if __name__ == '__main__':
  ############## GET OPTIONS ####################################
  allParam = getOptions()
  # unpack options
  (numChrom, numLoci, markerType, DSLafter, DSLdist, initSize, 
    meanInitAllele, burnin, introGen, minAlleleFreq, maxAlleleFreq, 
    fitness, mlSelModel, numSubPop, finalSize, noMigrGen,
    mixingGen, growth, migrModel, migrRate, mutaRate, recRate,
    dryrun, filename) = allParam
  #
  # load simuPOP libraries
  from simuPOP import *
  from simuUtil import *
  #
  # check minimal requirement of simuPOP version
  if simuRev() < 47  :
    raise exceptions.SystemError('''This scripts requires simuPOP revision >=47. 
      Please upgrade your simuPOP distribution.''' )
  #
  ################## RUN THE SIMULATION ###############
  simuComplexDisease( numChrom, numLoci, markerType, 
    DSLafter, DSLdist, initSize, meanInitAllele, burnin, introGen, 
    minAlleleFreq, maxAlleleFreq, fitness, mlSelModel, numSubPop, 
    finalSize, noMigrGen, mixingGen, growth, migrModel, 
    mutaRate, migrRate, recRate, dryrun, filename)
  print "Done!"
