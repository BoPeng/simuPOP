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
# Known limitations/bugs:
# 
#  * Currently, this script only handles constant selection pressure
#    and independent multi-locus selection model. 
#    

"""

Introduction
=============

NOTE: this version (> rev193) of simuComplexDisease.py uses a significantly
different mechanism to control the allele frequency of disease susceptibility
loci. I will describe the method briefly here. More details please see 
an unpublished paper. :-) If you prefer the previous method, you should check
out revision 193 of this file from the svn/web interface of simupop.sf.net.

This program simulates the evolution of a complex common disease under the 
influence of mutation, migration, recombination and population size change. 
Starting from a small founder population, each simulation will go through
the following steps:

  1. Simulate the trajectory of allele frequency using specified disease model.
  2. After determining mutant age, start simulation several thousands generations
     before the disease mutants are introduced. 
  3. Burn-in the population with mutation and recombination
  4. Introduce disease alleles and evolve the population with pre-simulated 
     allele frequency.
  5. Population structure and migration are specified along with demographic
     models. The population can be split into several equally-sized subpopulations
     and then evolve independently, or with migration. 

The result of the simulation is a large multi-generation population. To analyze 
the population, you will typically need to 

  1. Apply a penetrance function to the population and determine the affectedness
    for each individual

  2. Draw Population and/or pedigree based samples and save in popular 
     formats so that the samples can be analyzed by other software like
     genehunter.
     
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


Statistics Monitored
====================

A number of statistics will be measured and saved. They are:

  1. Fst before and after mixing
  2. Observed heterogeneity before and after mixing
  3. LD (D prime) between markers and between a DSL and surrounding markers
     at the end of each stage.
  4. Disease allele frequency trajectory.

"""

import simuOpt
import os, sys, exceptions, types, math

try:
  from rpy import *
  hasRPy = True
except:
  hasRPy = False

#
# declare all options. getParam will use these information to get parameters
# from a tk/wxPython-based dialog, command line, config file or user input
#
options = [
  {'arg': 'h',
   'longarg': 'help',
   'default': False, 
   'description': 'Print this usage message.',
   'allowedTypes': [types.NoneType, type(True)],
   'jump': -1          # if -h is specified, ignore any other parameters.
  },
  #
  # 
  {'separator': 'Genotype structure:'},  
  {'longarg': 'numChrom=',
   'default': 10,
   'configName': 'Number of chromosomes',
   'prompt': 'Number of chromosomes (10):  ',
   'description': 'Number of chromosomes.',
   'allowedTypes': [types.IntType],
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'numLoci=',
   'default': 20,
   'configName': 'Number of loci on each chrom',
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
   'default': [45, 85, 125],
   'configName': 'DSL after marker (0-indexed)',
   'prompt': 'Enter a list of 0-indexed disease loci ([45,85,125]):  ',
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
  #
  #
  {'separator': 'Demographic model:'},
  {'longarg': 'initSize=',
   'default': 10000,
   'configName': 'Initial population size',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Initial Population size (10000):  ',
   'description': '''Initial population size. This size will be maintained
        till the end of burnin stage''',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'endingSize=',
   'default': 200000,
   'configName': 'Final population size',
   'prompt': 'Final population size (sum of all subpopulations) (200000):  ',
   'allowedTypes': [types.IntType, types.LongType],
   'description': 'Final population size after population expansion.',
   'validate':  simuOpt.valueGT(0)
  }, 
  {'longarg': 'growthModel=',
   'default': 'exponential',
   'configName': 'Population growth model',
   'prompt': 'Population growth style, linear or exponential. (exponential):  ',
   'description': '''How population is grown from initSize to endingSize.
        Choose between linear and exponential''',
   'chooseOneOf': ['exponential', 'linear'],
  },
  {'longarg': 'burninGen=',
   'default': 3000,
   'configName': 'Length of burn-in stage',
   'allowedTypes': [types.IntType],
   'prompt': 'Length of burn in stage (3000):  ',
   'description': 'Number of generations of the burn in stage.',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'splitGen=',
   'default': 5000,
   'configName': 'When to split population',
   'prompt': 'When to split population (5000):  ',
   'allowedTypes': [types.IntType, types.LongType],
   'description': '''At which generation to split the population. 
        The population will start to grow after burnin stage but will
        not split till this generation. Note that if the disease is 
        introduced after this stage, it will be in one of subpopulations.
   ''',
   'validate':  simuOpt.valueGT(0)
  },
  {'longarg': 'mixingGen=',
   'default': 8000,
   'configName': 'When to mix population',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'When to mix population: (8000):  ',
   'description': '''At which generation to start mixing (allow migration.
        This number should be greater than or equal to split gen.''',
   'validate':  simuOpt.valueGE(0)
  },
  {'longarg': 'endingGen=',
   'default': 10000,
   'configName': 'Ending generation number',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'When to stop the simulation (10000):  ',
   'description': '''At which generation to stop the simulation.
        This is the total generation number.''',
   'validate':  simuOpt.valueGE(0)
  },
  #
  #
  {'separator': 'Migration parameters:'},  
  {'longarg': 'numSubPop=',
   'default': 5,
   'configName': 'Number of subpopulations to split',
   'allowedTypes': [types.IntType],
   'prompt': 'Number of subpopulations to split (5):  ',
   'description': 'Number of subpopulations to be split into after burnin stage.',
   'validate':  simuOpt.valueGT(0)
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
  {'longarg': 'alleleDistInSubPop=',
   'default': 'uneven',
   'configName': 'Allele distribution in subpopulations',
   'prompt': 'Allele distribution in subpopulations (exponential):  ',
   'description': '''If disease allele frequencies in each subpopulation are
        not specified, how to distribute disease alleles in subpopulations.
        If 'even' is chose, number of disease alleles will be distributed
        according to a multinomial distribution so is largely even among 
        subpopulations. If 'exponential' is chosen, the number of disease alleles
        will be proportional to the interval lengths of 0 x x x 1 while x are uniform
        [0,1]. The distribution of interval lengths, are roughly exponential 
        (conditional on overall length 1). ''',
   'validate':  simuOpt.valueOneOf(['even', 'uneven']),
   'chooseOneOf': ['even', 'uneven']
  },
  #
  #
  {'separator': 'Disease model:'},
  {'longarg': 'curAlleleFreq=',
   'default': [0.05, 0.05, 0.05],
   'configName': 'Final allele frequencies',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Allele frequency at the ending generation (0.05,0.05,0.05):  ',
   'description': '''Current allele frequencies for each DSL.
        If a number is given, it is assumed to be the frequency
        for all DSL.''',
   'validate': simuOpt.valueListOf( simuOpt.valueBetween(0,1))
  },
  {'longarg': 'minMutAge=',
   'default': 0,
   'configName': 'Minimum mutant age',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Minimum mutation age (0)',
   'description': '''Minimum mutation age. Default to 0. Because the population
        may be split into subpopulations at splitGen, if the mutation age is too
        short, it may fall in one of the subpopulations. To void this, you can 
        set minMutAge to be endingGen - splitGen. Note that all mutants have the same
        minimum mutant age regardless final allele frequency.'''
  },
  {'longarg': 'maxMutAge=',
   'default': 0,
   'configName': 'Maximum mutant age',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Maximum mutation age (0)',
   'description': '''Maximum mutant age. Default to 0, means no max and the real
        maximum mutant age will be endingGen. However, if you do not want mutant to
        be generated in the burnin stage. You can set maxMutAge to endingGen-burnin.
        Note that all mutants have the same maximum mutant age regardless final 
        allele frequency.'''
  },
  {'longarg': 'fitness=',
   'default': [1, 1.0001, 1.0002],
   'configName': 'Fitness of genotype AA,Aa,aa',
   'allowedTypes': [types.ListType, types.TupleType],
   'prompt': 'Fitness of genotype AA,Aa,aa ([1, 1.0001, 1.0002]): ',
   'description': '''Fitness of genotype AA,Aa,aa in the monogenic case, and 
        [AA Aa aa BB Bb bb ...] in the polygenic case. ''',
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
  #
  #
  {'separator': 'Evolutionary parameters:'},
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
   'default': [0.0005],
   'configName': 'Recombination rate',
   'allowedTypes': [types.TupleType, types.ListType],
   'prompt': 'Recombination rate(s). (0.0005): ',
   'description': '''If a number is given, it is the recombination rate between
         adjacent markers. If a list is given, it should be the recombination
         rate between all markers and DSL. For example, if you have two chromosome
         with 5 markers, and with DSL after 3. The marker/DSL layout is
           0 1 2 3 x 5    6 7 8 9 10  
         The specified recombination rate should be those between 
           0-1, 1-2, 2-3, 3-x, x-5, 6-7, 7-8, 8-9, 9-10.
         Note that the distance between 3-x, x-5 is smaller than distance
         between markers.
   ''',
   'validate':  simuOpt.valueListOf(simuOpt.valueBetween(0,1))
  },
  #
  # 
  {'separator': 'Miscellaneous:'},
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
          name + .bin/txt/xml: saved popuation''',
  },
  {'longarg': 'saveFormat=',
   'default': 'txt',
   'allowedTypes': [types.StringType],
   'configName': 'Format to save population',
   'prompt': 'Format to save population):  ',
   'description': '''Format to save population, can be 
         text, bin or xml. Note that the binary format, although
         smallest, may not be portable between different machines.''',  
   'chooseOneOf': [ 'txt', 'bin', 'xml']
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
  # change current directory, all files wil be saved here.
  try:
    os.mkdir(allParam[-3])
  except:
    pass
  # --saveConfig
  if allParam[-3] != None: # saveConfig
    simuOpt.saveConfig(options, os.path.join(allParam[-3], allParam[-3]+'.cfg'), allParam)
  # --verbose or -v (these is no beautifying of [floats]
  if allParam[-1]:         # verbose
    simuOpt.printConfig(options, allParam)
  # return the rest of the parameters
  return allParam[1:-1]

def plotScenario(name, burninGen, splitGen, mixingGen, 
  endingGen, popSizeFunc, trajectory, legend):
  ''' plot the simulation scenario '''
  if not hasRPy:
    print "Module rpy is not present. Can not draw figure."
    return
  r.postscript(name+'.eps')
  r.par(mar=[5.1, 4.1, 4.1, 4.1])
  sg = max(min( endingGen - max([len(x) for x in trajectory]) - 1000, burninGen-1000 ), 0)
  r.plot( x=[0, burninGen], y=[sum(popSizeFunc(0)), sum(popSizeFunc(burninGen))],
    xlim=[sg, endingGen], ylim=[0, sum(popSizeFunc(endingGen))],
    xlab='Generations', ylab='Population size', type='l',
    main='Demographic model and trajectory of disease allele frequency of simulation ' + name, 
    axes=False)
  r.box()
  r.axis(1, [burninGen, splitGen, mixingGen, endingGen],
    [str(x) for x in  [burninGen, splitGen, mixingGen, endingGen]])
  r.axis(2)
  r.abline(v=[burninGen, splitGen, mixingGen, endingGen], lty=3)
  r.lines(x=[0,endingGen], y=[0,0], lty=1)
  r.lines( x=range(burninGen, splitGen), y=[sum(popSizeFunc(x)) for x in range(burninGen, splitGen)],
    type='l')
  # split
  if mixingGen > splitGen:
    for sp in range( len( popSizeFunc(endingGen) )):
      r.lines( x= range(splitGen, mixingGen),
        y=[ sum(popSizeFunc(x)[:(sp+1)]) for x in range(splitGen, mixingGen) ],
        lty=1
    )
  if endingGen > mixingGen:
    for sp in range( len( popSizeFunc(endingGen) )):
      if sp != len( popSizeFunc(endingGen) ) -1:
        lineType = 2
      else:
        lineType = 1
      r.lines( x= range(mixingGen, endingGen),
        y=[ sum(popSizeFunc(x)[:(sp+1)]) for x in range(mixingGen, endingGen) ],
        lty=lineType
    )
  # ploting trajectories
  r.par(new=True)
  maxAF = max( [ max(x) for x in trajectory] )
  r.plot(0, 0, type='n', xlim=[sg, endingGen], ylim=[0,maxAF],
    axes=False, xlab='', ylab='')
  r.axis(4)
  r.mtext('Allele frequency', 4, line=2)
  col = 1
  for traj in trajectory:
    r.lines( x = range( endingGen - len(traj), endingGen ),
      y=traj, col=col)
    col += 1
  r.text( x=sg, y=maxAF, labels=legend, adj=[0,1])
  r.dev_off()
  
  
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
    mutaRate, recRate, 
    dryrun, filename, format):
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
  else:
    if fitnessTmp == []:  # neutral process
      fitness = [1,1,1]*numDSL
    else:
      # for a single DSL
      if len(fitnessTmp) == 3:
        fitness = fitnessTmp*numDSL
      elif len(fitnessTmp) != numDSL*3:
        raise exceptions.ValueError("Please specify fitness for each DSL")
      else:
        fitness = fitnessTmp
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
    else:  # exponential
      rate =  (math.log(endingSize)-math.log(initSize))/(endingGen-burninGen)
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
  plotScenario(filename, burninGen, splitGen, mixingGen,
    endingGen, popSizeFunc, traj, '''Initial pop size: %d
Final pop size: %d
Burnin: %d
Split at gen: %d
Mix at gen: %d
End at gen: %d
Number of subPop: %d

Current allele freq: %s
Fitness: %s
Min mutant age: %d
Max mutant age: %d ''' % \
(initSize, endingSize, burninGen, splitGen, mixingGen, endingGen, numSubPop, \
','.join([str(x) for x in curAlleleFreqTmp]), \
','.join([str(x) for x in fitnessTmp]), minMutAge, maxMutAge) )
  print '\n\nTrajectories simulated successfully. Their lengths are %s.\n Check %s.eps for details. \n\n' \
    % ( ', '.join([ str(len(x)) for x in traj]), filename)
  ###
  ### translate numChrom, numLoci, DSLafterLoci to
  ### loci, lociPos, DSL, nonDSL in the usual index
  ### (count DSL along with markers)
  ###
  if markerType == 'microsatellite':
    maxAle = 99                   # max allele
  else:
    maxAle = 1                    # SNP (0 and 1)
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
  # non-DSL loci, for convenience
  nonDSL = range(0, reduce(operator.add, loci))
  for loc in DSL:
    nonDSL.remove(loc)
  ###
  ### initialization 
  ###
  if maxAle > 1:  # Not SNP
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
  if maxAle > 1:  # Not SNP
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
      step=10),
    # output to file (append)
    pyEval( expr=r'"%d %d " %(gen, popSize) + " ".join(["%.5f"%(1-alleleFreq[x][0]) for x in DSL])+"\n"',
      output='>>>'+filename+'.traj')
  ]
  ###
  ### introduction of disease mutants
  ###
  ###                 endingGen
  ###      0 1 ...... i_T
  for i in range( numDSL ):
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
  try:
    mlSelModel = {'additive':SEL_Additive, 
      'multiplicative':SEL_Multiplicative,
      'none':SEL_None}[mlSelModelTmp]
  except:
    raise exceptions.ValueError("Wrong multi-locus seleciton model. " + mlSelModelTmp)
  
  if mlSelModel != SEL_None:
    operators.append( mlSelector(
      # with five multiple-allele selector as parameter
      [ maSelector(locus=DSL[x], wildtype=[0], 
        fitness=[fitness[3*x],fitness[3*x+1],fitness[3*x+2]]) for x in range(len(DSL)) ],
      mode=mlSelModel, begin=splitGen),
    )
  ###
  ### migration
  ###
  if numSubPop > 1 and migrModel == 'island' and migrRate > 0:
    operators.append( migrator(migrIslandRates(migrRate, numSubPop),
      mode=MigrByProbability, begin=mixingGen) )
  elif numSubPop > 1 and migrModel == 'stepping stone' and migrRate > 0:
    operators.append( migrator(migrSteppingStoneRates(migrRate, numSubPop, 
      circular=True),  mode=MigrByProbability, begin=mixingGen) )
  ###
  ### prepare for sampling:
  ###
  operators.extend([
    # save ancestral populations starting at -2 gen
    setAncestralDepth(1, at=[-2]),
    # track pedigree
    parentsTagger(begin=-2),
  ])
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
  ### 
  ### prepare mating scheme
  ###
  ### use last_two numOffspringFunc, simuPOP will produce 2 offspring at 
  ### the last two generations
  ###
  def last_two(gen):
    if gen >= endingGen -2:
      return 2
    else:
      return 1
  # create a simulator
  pop =  population(subPop=popSizeFunc(0), ploidy=2,
    loci = loci, maxAllele = maxAle, lociPos = lociPos,
    infoFields=['fitness'])
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
      newSubPopSizeFunc=popSizeFunc,      # demographic model
      numOffspringFunc=last_two,          # save last two generations
      loci=DSL,                           # which loci to control
      alleles=[1]*numDSL,                 # which allele to control
      freqFunc=trajFunc                   # frequency control function
    ),
    rep=1)
  # evolve! If --dryrun is set, only show info
  simu.evolve( preOps = preOperators, ops = operators, 
    end=endingGen, dryrun=dryrun )
  if dryrun:
    raise exceptions.SystemError("Stop since in dryrun mode.")
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
    mutaRate, recRate, 
    dryrun, filename, format) = allParam
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
  if simuRev() < 199  :
    raise exceptions.SystemError('''This scripts requires simuPOP revision >=199. 
      Please upgrade your simuPOP distribution.''' )
  #
  ################## RUN THE SIMULATION ###############
  simuComplexDisease(numChrom, numLoci, markerType, DSLafter, DSLdist, 
    initSize, endingSize, growthModel, 
    burninGen, splitGen, mixingGen, endingGen, 
    numSubPop, migrModel, migrRate, alleleDistInSubPop, 
    curAlleleFreq, minMutAge, maxMutAge, fitness, selMultiLocusModel, 
    mutaRate, recRate, 
    dryrun, os.path.join(filename, filename), format)
  
  print "Done!"
