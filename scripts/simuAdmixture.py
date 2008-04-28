#!/usr/bin/env python
'''
This script simulates an admixed population based on the HapMap dataset. Please
read this help message carefully, making sure your know how this script works,
then run a few test commands, before you explore the capacity of this script.


Step 0: Prepare HapMap dataset (call scripts/loadHapMap.py if available)
=========================================================================

This script makes use of the HapMap dataset. The dataset is downloaded, imported
and saved in simuPOP format automatically, using script scripts/loadHapMap.py.
If loadHapMap.py can not be imported (not in the working directory or in $PYTHONPATH),
please try to run loadHapMap.py manually and provide a path with files hapmap_XX.bin
to parameter --HapMap_dir.


Step 1: Generate a seed population (scripts/simuAdmixture.py)
==============================================================

Determine which populations and markers to use.
-----------------------------------------------

The initial population consists of one, two or three hapmap populations
(param --pops), and a selected set of markers. You can start from 
all hapmap markers or a marker list in the format of 
'chromosome position name', such as (no header is required).

1 2103664 rs1496555
1 2708522 rs1338382
1 2719853 rs10492936
1 2734227 rs10489589
...

If more fields are given, they will be ignored.

You can refine your selection using 
  1. which chromosome(s) to use (param --chrom)
  2. number of markers to use on each chromosome (param --numMarkers)
  3. starting and ending positions (param --startPos and --endPos)
  4. minimal allele frequency (param --minAF)
  5. minimal allele frequency differences between hapmap populations
    (param --minDiffAF). This criteria can be used to find markers
    of high ancestray information content.
  6. minimal distance between markers (param --minDist)

Note that 'position' in the marker list is assumed to be in base pair
(as in Affymetrix or Illumina annotation data). Marker positions 
are converted to centiMorgan by dividing the positions by 10^6 (i.e. 
1 centiMorgan = 1 M basepairs). Starting and ending position should be
inputted in centiMorgan. Not all requirements can be met at the same time.

If an initial population is given, it will be used directly (without
checking parameters), and to save the generated initial population.

Each simulation will be given a name and all files will be saved to a directory
named after it. A file markers.lst will be saved to directory $name.


Generate a seed population. 
----------------------------

The population is expanded instantly by copying individual 10 (param --initCopy)
time to avoid quick loss of heterogenity due to small population sizes. It
then evolves for a relatively short period of time (param --initGen )
without migration. This stage is designed to mix these copied individuals,
and make HapMap populations more distinct, because the HapMap populations are
already admixed.

During the evolution, recombination at a rate of 0.01 per cM (param --recIntensity),
and mutation at a rate of 1e-7 per nucleotide per generation (param --mutaRate) 
is applied. Recombination and mutation will also be applied at the subsequent stages.

This population will be saved as 'seed.bin' (param --seed) that will be used
for the subsequent replicate simulations. We consider the subpopulations of
the seed population as the populations around 2000 years ago, that have 
been evolved separately from the Out-Of-Africa population for abut 40000
generations



Step 2. Evolve the seed population.
=====================================

After a seed population is generated (or loaded if already exists and parameter
--useSavedSeed is given), this script will evolve it as follows:

Evolve the seed population subject to rapid population expansion for
100 (param --expandGen) generations. This is to mimic the rapid
population expansion in the past 2000 years. Exponential population 
growth model is used. No migration is allowed.

Selection on a number of loci can be applied to selected loci at this, and
the next mixing stage. These loci are supposed to have disease predisposing
alleles (though selection does not have to work against them), and there
frequencies are controlled.

A special controlled mating scheme can be used if the allele frequencies
of some markers at the end of this stage is specified. This can be done
backward in time, assuming no disease allele at these loci in 
the seed population (c.f. Peng 2007 PLoS genetics), or forward in time,
by starting from existing allele frequencies in the seed population.
These two kinds of loci can not yet be mixed.

The resulting population of this generation will be saved as $name/expanded.bin
(--param expandedName).


Step 3. Mix the subpopulations (optional)
==========================================

Mix the populations using given migration model. Two basic models,
namely 'Hybrid Isolations' and 'Continuous Gene Flow' are provided. 
The first model simply mix the subpopulations and form a single population.
The second model allows stable migration between subpopulations during
the migration period. More advanced migration model, with changing
migration rates can be specified using a 'Customized' model. To use 
the last model, you will have to modify the 'migrFunc' function in
this script. Note that it is possible to create a separate population,
with individuals migrated from existing populations.

e.g.

'Hybrid Isolations':  Merge all HapMap populations instantly.

'Continuous Gene Flow', two HapMap populations, with migration rate
    [[1, 0],
     [0.1, 0.9]]
  10% of individuals will migrate from population 1 to 0 at each 
  generation. Note that the population sizes will change as a result
  of migration so the exact number of migrants varies from generation
  to generation.
  
'Continuous Gene Flow', three HapMap populations, with migration rate
    [[1, 0, 0],
     [0.1, 0.9, 0],
     [0.1, 0, 0.9]]
   10% of individuals from population 2 and 3 migrate to population 1.

'Continuous Gene Flow', two HapMap populations, with migration rate
    [[0.8, 0, 0.2],
     [0, 0.8, 0.2]]
    A new population is created, and gets 20% of individuals from
    population 1 and 2 at each generation.


The result of this stage will be saved to $name/admixed.bin (--param admixedName)

The script will skip population expansion and/or admixture if the corresponding
file is found and parameter useSavedExpanded and/or useSavedAdmixed is
specified. In the useSavedExpanded case, the simulation starts mixing population
$name/expanded.bin directly. In the useSavedAdmixed case, the simulation draw
sample from the admixed population directly.

The result of this script is expanded.bin (if no migration) or admixed.bin.
Users can apply arbitrary penetrance or quantitative trait model and draw
samples from this population. Script analAdmixture.py implements the examples
used in paper Peng 2008.


Test scripts
==============

The following test scripts demonstrate the use of this script using a small
number of loci. Note that some parameters can be ignored if their
default values are used.

simuAdmixture.py --noDialog  --HapMap_dir='../../HapMap' \
    --chrom="range(1,3)"  --markerList='../../Affy/mapAffySNPs.txt' \
    --startPos="[0]" --endingPos='[0]' --numMarkers="[100,100]" --minAF='0' --minDist='0'  \
    --pops="['CEU', 'YRI', 'JPT+CHB']" --initCopy='10' --initGen='20' --seedSize='3600'  \
    --seedName='test_seed.bin'

simuAdmixture.py --noDialog  --useSavedSeed --seedName=test_seed.bin --migrGen='5' \
    --migrRate="([1, 0, 0], [1, 0, 0], [1, 0, 0])" 
    

simulation for XJ Gu et al (2008)
====================================

# 1: Hybrid isolation.
simuAdmixture.py --noDialog  --name='IH' \
    --HapMap_dir='../../HapMap' --pops="['CEU']" --chrom="range(1, 23)" \
    --markerList='../../Affy/mapAffySNPs.txt'  --initCopy='10' \
    --gen='200' --size='4800' --expandGen='100' --expandSize='24000' \
    --useSavedAdmixed --migrModel='None' --migrGen='1' \
    --migrRate="()" --chromWithDSL="(1, 2, 3, 4)" \
    --freqDSL='0.2' --freqDev='0.02' --dslVar="(0.005, 0.01, 0.03, 0.05)" \
    --cutoff='-0.5' --DSLpene='[]' --peneFunc='None' --parameter='[0.5]' --ccSampleSize="(600, 600)" \
    --ccSampleName='case-control' --randomSampleSize='800' --randomSampleName='random'

simuAdmixture.py --noDialog  --name='admix' --useSavedSeed --initPop='' --HapMap_dir='../../HapMap' \
    --pops="['CEU', 'YRI', 'JPT+CHB']" --chrom="(2, 3)" --markerList='' --numMarkers="(1000, 1000)" \
    --startPos='[0]' --endingPos='[0]' --minAF='0.05' --minDiffAF='0' --minDist='0' --initCopy='10' \
    --gen='200' --size='4800' --useSavedExpanded --expandGen='100' --expandSize='24000' \
    --useSavedAdmixed --migrModel='Continuous Gene Flow' --migrGen='5' \
    --migrRate="([0.90000000000000002, 0.10000000000000001], [0.0, 1.0])" --chromWithDSL="(1, 2, 3, 4)" \
    --freqDSL='0.2' --freqDev='0.02' --dslVar="(0.050000000000000003, 0.10000000000000001, 0.29999999999999999, 0.5)" \
    --cutoff='-0.5' --DSLpene='[]' --peneFunc='None' --parameter='[0.5]' --ccSampleSize="(600, 600)" \
    --ccSampleName='case-control' --randomSampleSize='800' --randomSampleName='random'


Evolve seed population:

simuAdmixture.py --noDialog  --expandGen='100' --expandSize='24000' \
    --migrModel='Continuous Gene Flow' --migrGen='5' \
    --migrRate="([0, 0.10000000000000001], [0.0, 1.0])" \
    --DSLpene='[0]' --pene="(0.10000000000000001, 0.25, 0.5)"  \
    --resample='False' --sampleSize="(500, 500)"  \
    --name='admix'

Or use some default values

simuAdmixture.py --noDialog  --seed=admix_seed.bin --migrGen='5' \
    --migrRate="([1, 0, 0], [1, 0, 0], [1, 0, 0])" --chromWithDSL="[1,2,3,4]" \
    --sampleSize="(600, 600)" --freqDSL=0.15 --freqDev=0.01 --cutoff=1 \
    --dslVar="[0.005, 0.01, 0.03, 0.05]"  --name='admix'

# test
simuAdmixture.py --noDialog  --seed=test_seed.bin --migrGen='5' \
    --migrRate="([1, 0, 0], [1, 0, 0], [1, 0, 0])" --chromWithDSL="[1,2]" \
    --sampleSize="(5, 5)" --freqDSL=0.1 --cutoff=1 \
    --dslVar="[0.05, 0.1]"   --name='test'

Evolve a different seed population

simuAdmixture.py --noDialog  -s=seed1.bin --migrGen='5' \
    --migrRate="([0, 0.1], [0, 1])" --DSLpene=100 \
    --pene="(0.1, 0.25, 0.5)"  --sampleSize="(500, 500)"  \
    --name='admix'

Resample from the saved population:

simuAdmixture.py --noDialog --resample=True --DSLpene=100 \
    --pene="(0.1, 0.25, 0.5)"  --sampleSize="(500, 500)"  \
    --name='admix'

Two disease susceptibility loci:

simuAdmixture.py --noDialog --resample=True --DSLpene='[10, 50]' \
    --pene="(0.1, 0.25, 0.5, 0.1, 0.25, .5, 0.2, 0.4, 0.6)" \
    --sampleSize="(500, 500)"   --name='admix'

Create a new subpopulation (the third population is created):

simuAdmixture.py --noDialog --migrRate='[[0, 0, 0.1], [0, 0, 0.1]]' \
    --DSLpene='[10, 50]' \
    --pene="(0.1, 0.25, 0.5, 0.1, 0.25, .5, 0.2, 0.4, 0.6)" \
    --sampleSize="(500, 500)"   --name='admix'

Sample only from the last (newly created population)

simuAdmixture.py --noDialog --migrRate='[[0, 0, 0.1], [0, 0, 0.1]]' \
    --DSLpene='[10, 50]' \
    --pene="(0.1, 0.25, 0.5, 0.1, 0.25, .5, 0.2, 0.4, 0.6)" \
    --sampleSize="([0,0,500], [0,0,500])"   --name='admix'

Use a varying migration rate model:

Modify migrFunc to fit your need. Then

simuAdmixture.py --noDialog --migrModel='Customized' --DSLpene='[10, 50]' \
    --pene="(0.1, 0.25, 0.5, 0.1, 0.25, .5, 0.2, 0.4, 0.6)" \
    --sampleSize="([0,0,500], [0,0,500])" --name='admix'


Dr. Reddon simulation one:

Seed population: 

500 markers from chromosome 2, initial allele frequency > 0.1, initial 
allele frequency difference between CEU and YRI populations > 0.2, 
minimal distance between adjacent markers 0.05cM

simuAdmixture.py --noDialog  --HapMap_dir='../../HapMap' --chrom='[2]' \
   --numMarkers=500 --startPos=100  --minAF=0.1 --minDiffAF=0.2 \
   --minDist=0.05 --pops="['CEU', 'YRI']"

250 sample from CEU, 250 sample from YRI, 500 sample from admixed population
when 10% of the CEU popopulation migrate to YRI for 5 generations.

Round 1: expand and get sample from CEU

simuAdmixture.py  --noDialog --useSavedExpanded --migrModel='None' --migrGen='0'

Round 2: load expanded population and get sample from YRI

simuAdmixture.py  --noDialog --remix=True --migrModel='None' --migrGen='0' --sampleType='random'\
  --sampleSize="(0, 250)" --sampleName='YRI' --name='simu2'
  
'''

from simuOpt import *
setOptions(alleleType='binary')
from simuPOP import *
from hapMapUtil import getMarkersFromName, getMarkersFromRange

import os, sys, math 
from types import *
from exceptions import ValueError, SystemError
from simuUtil import SaveQTDT

try:
    from rpy import *
    hasRPy = True
except:
    hasRPy = False

HapMap_pops = ['CEU', 'YRI', 'JPT+CHB']

options = [
    {'arg': 'h',
     'longarg': 'help',
     'default': False, 
     'description': 'Print this usage message.',
     'allowedTypes': [NoneType, type(True)],
     'jump': -1                    # if -h is specified, ignore any other parameters.
    },
    {'longarg': 'name=',
     'default': 'simu',
     'useDefault': True,
     'allowedTypes': [StringType],
     'label': 'Name of the simulation',
     'description': '''Name of this simulation. A directory with this name
                will be created. Configuration file (.cfg), marker list and
                various populations will be saved to this directory''',
    },
    {'longarg': 'useSavedSeed',
     'default': False,
     'label': 'Use saved seed population',
     'useDefault': True,
     'allowedTypes': [BooleanType],
     'jump': 12,
     'description': '''Use specified or a default seed population, if available.
        The default seed population is seed.bin under the simulation directory.
        '''
    },
    {'longarg': 'useSavedExpanded',
     'default': False,
     'useDefault': True,
     'allowedTypes': [BooleanType],
     'jump': 16,
     'label': 'Use saved expanded population',
     'description': '''If set to true, load specified or saved $name/expanded.bin and 
                skip population expansion'''
    },
    #
    {'separator': 'Progress report and visualization'},
    {'longarg': 'step=',
     'default': 10,
     'label': 'Progress report interval',
     'useDefault': True,
     'allowedTypes': [IntType, LongType],
     'description': '''Gap between generations at which population statistics are 
                calculated and reported.'''
    },
    {'longarg': 'showAlleleFreq',
     'default': True,
     'useDefault': True,
     'allowedTypes': [BooleanType],
     'label': 'Show allele frequency at specified loci',
     'description': '''If set, display allele frequency of loci specified in parameters
                --controlledLoci or --backwardControlledLoci''',
    },
    {'longarg': 'figureStep=',
     'default': 20,
     'label': 'Figure update interval',
     'useDefault': True,
     'allowedTypes': [IntType, LongType],
     'description': '''Gap between generations at which LD plots are 
                draw. Default to 20.'''
    },
    {'longarg': 'showLDPlot',
     'default': False,
     'useDefault': True,
     'allowedTypes': [BooleanType],
     'label': 'Show allele frequency at specified loci',
     'description': '''If set, draw and display LD structure (using haploview, available)
                on specified regions (--ldRegions). The figures will be saved in names
                such as CEU_stage_start-end_gen.PNG''',
    },
    {'longarg': 'haploview=',
     'default': 'haploview',
     'useDefault': True,
     'allowedTypes': [StringType],
     'label': 'Command to start haploview',
     'description': '''Path to haploview or command to start haploview. It can simply
                be haploview, but can be something like '/path/to/jave /path/to/haploview.jar'.
                If haploview is not found, no LD plot will be displayed.'''
    },
    {'longarg': 'ldRegions=',
     'default': [],
     'useDefault': True,
     'allowedTypes': [TupleType, ListType],
     'label': 'Regions to plot LD structure',
     'description': '''A list of regions, in terms of marker position that LD structure
                is plotted. If there are 1000 markers on each chromosome, viewRegions
                can be [[0, 500], [1200, 1700]]. A single region such as [0, 500] is
                allowed.'''
    },
    {'longarg': 'ldSampleSize=',
     'default': 200,
     'label': 'Sample size for ld plotting',
     'useDefault': True,
     'allowedTypes': [IntType, LongType],
     'description': '''A random sample of specified size will be draw from each subpopulation
                for haploview plot.'''
    },
    #
    {'separator': 'Populations and markers to use'},
    {'longarg': 'HapMap_dir=',
     'default': '../HapMap',
     'useDefault': True,
     'label': 'HapMap data directory',
     'description': '''Directory to store HapMap data in simuPOP format. Hapmap
                data file hapmap_??.bin in this directory, if exits, will be
                loaded directly. Otherwise, module loadHapMap.py (usually under
                /path/to/simuPOP/scripts) will be used to download, import, and
                save HapMap data in simuPOP formats. If this module can not be
                imported, you can either add its path to environmental variable
                $PYTHONPATH or run this script manually.''',
     'allowedTypes': [StringType],
     #'validate': valueValidDir(),
    },
    {'longarg': 'pops=',
     'default' : ['CEU', 'YRI', 'JPT+CHB'],
     'useDefault': True,
     'label' : 'HapMap populations',
     'description': '''Which HapMap populations to use?''',
     'allowedTypes': [TupleType, ListType],
     'chooseFrom': HapMap_pops,
     'validate': valueListOf(valueOneOf(HapMap_pops)),
    },
    {'longarg': 'markerList=',
     'default': '',
     'useDefault': True,
     'label': 'Marker list file',
     'description': '''A file with a list of marker names, in the form of
                "chrom_number marler_pos marker_name". Markers that on a chromosome not
                in the chromosome list (parameter chrom) are ignored. lines start with 
                # is ignored. If numMarkers, startPos, endingPos, minDist are also specified, 
                the first numMarkers between startPos and endingPos will be used.
                This script assumes that the marker position in the 
                list file is in base pair, and will use pos/1000000 as cM to
                compare marker location. If more fields are given, others are ignored.''',
     'allowedTypes': [StringType],
     'validate': valueOr(valueEqual(''), valueValidFile()),
    },
    {'longarg': 'chrom=',
     'default': [2, 3],
     'label': 'Chromosomes to use',
     'description': '''A list of chromosomes to use from the HapMap data. When multiple
                chromosomes are involves, numMarkers, if used, should be a list that specicy
                number of markers on each chromosome. The same rule applies to startPos
                and endingPos as well.''',
     'allowedTypes': [TupleType, ListType],
     'validate': valueListOf(valueBetween(1, 22)),
    },
    {'longarg': 'numMarkers=',
     'default': [1000, 1000],
     'label': 'Number of markers to use',
     'description': '''Number of markers to use from the marker list file. If 0 is used,
                all markers that satisfy conditions startPos, endingPos, minDist will
                be used.''',
     'allowedTypes': [TupleType, ListType],
     'validate': valueOr(valueGT(0), valueListOf(valueGE(0)))
    },
    {'longarg': 'startPos=',
     'default': 0,
     'useDefault': True,
     'label': 'staring position',
     'description': '''Starting position of the markers. If multiple
                chromosomes are used, the positions for each 
                chromosome can be specified as a list.''',
     'allowedTypes': [TupleType, ListType],
     'validate': valueOr(valueGE(0), valueListOf(valueGE(0)))
    },
    {'longarg': 'endingPos=',
     'default': 0,
     'useDefault': True,
     'label': 'Ending position',
     'description': '''Ending position of the markers. Ignored if its value 
                is 0.  If multiple chromosomes are used, the positions for each 
                chromosome can be specified as a list. ''',
     'allowedTypes': [TupleType, ListType],
     'validate': valueOr(valueGE(0), valueListOf(valueGE(0)))
    },
    {'longarg': 'minAF=',
     'default': 0.05,
     'useDefault': True,
     'label': 'Minimal allele frequency',
     'description': '''Minimal allele frequency, only used for picking markers
                from the HapMap dataset''',
     'allowedTypes': [IntType, LongType, FloatType],
     'validate': valueGE(0)
    },
    {'longarg': 'minDiffAF=',
     'default': 0,
     'useDefault': True,
     'label': 'Minimal allele frequency difference',
     'description': '''Minimal allele frequency difference between two HapMap population,
                , can only be used when two HapMap populations are used. This options can be
                used to choose markers with high ancestry information content.''',
     'allowedTypes': [IntType, LongType, FloatType],
     'validate': valueGE(0)
    },
    {'longarg': 'minDist=',
     'default': 0,
     'useDefault': True,
     'label': 'Minimal distance between markers (cM)',
     'allowedTypes': [IntType, LongType, FloatType],
     'description': '''Minimal distance between markers (in the unit of cM).
                Can be used for both methods.''',
    },
    #
    {'separator': 'Generate seed population'},
    {'arg': 's:',
     'longarg': 'seedName=',
     'default': 'seed.bin',
     'useDefault': True,
     'description': '''Name of the seed population''',
     'allowedTypes': [StringType],
    },
    {'longarg': 'initCopy=',
     'default': 10,
     'useDefault': True,
     'label': 'Initial propagation',
     'description': '''How to expand the initial small HapMap sample to
                 avoid quick loss of heterogenity. By default, each individual
                 is copied 10 times.''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'initGen=',
     'default': 100,
     'useDefault': True,
     'label': 'Generations to evolve',
     'description': '''Number of generations to evolve to get the seed
                population.''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'seedSize=',
     'default': 4800,
     'useDefault': True,
     'label': 'Size of the seed population',
     'description': '''Size of the seed population. The default value is the recommended
                value when all hapmap populations are used (60+60+90)*20. You may want
                to reduce it according to the populations used.''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGE(100)
    },
    #
    {'separator': 'Mutation, selection and recombination'},
    {'longarg': 'mutaRate=',
     'default': 1e-6,
     'label': 'Mutation rate',
     'allowedTypes': [IntType, FloatType],
     'description': '''Mutation rate using a 2-allele model (kam). Note that mutation
                can generally be ignored during short period of time unless you
                intentionally set a higher mutation rate.''',
     'validate': valueBetween(0,1),
    },
    {'longarg': 'recIntensity=',
     'default': 0.01,
     'label': 'Recombination intensity',
     'allowedTypes': [FloatType],
     'description': '''Recombination intensity. The actually recombination rate between
                two adjacent markers depends on distance (in cM) between them. For example,
                two markers that are 10kb apart (0.00001 cM apart) will have recombination
                rate 10^-5*0.01 (the default value) = 10^-6.
     ''',
     'validate': valueBetween(0,1),
    },
    {'longarg': 'controlledLoci=',
     'label': 'Loci with controlled allele frequency',
     'default': [],
     'useDefault': True,
     'allowedTypes': [TupleType, ListType],
     'description': '''A list of markers (by name) whose allele frequency will be
                controlled during this stage of evolution. A forward-time trajectory
                simulation algorithm will be used. Currently, only one of
                --controlledLoci and --backwardControlledLoci is allowed.'''
    },
    {'longarg': 'controlledFreq=',
     'label': 'Ending allele frequency at controlled loci',
     'default': [],
     'useDefault': True,
     'allowedTypes': [TupleType, ListType],
     'description': '''A list of allele frequency ranges for each controlled locus.
                If a single range is given, it is assumed for all markers. An example
                of the parameter is [[0.18, 0.20], [0.09, 0.11]].'''
    },
    {'longarg': 'backwardControlledLoci=',
     'label': 'Backward controlled loci',
     'default': [],
     'useDefault': True,
     'allowedTypes': [TupleType, ListType],
     'description': '''A list of markers (by name) whose mutants, if any, will be removed
                at the beginning of population expansion stage. A mutant will be introduced
                as the result of mutation. The frequency trajectory will be simulated
                using a backward approach (see Peng 2007, PLoS Genetics). Currently,
                only one of --controlledLoci and --backwardControlledLoci is allowed.''',
    },
    {'longarg': 'backControlledFreq=',
     'label': 'Ending allele frequency at backward controlled loci',
     'default': [],
     'useDefault': True,
     'allowedTypes': [TupleType, ListType],
     'description': '''A list of allele frequency (not a list of ranges as parameter controlledFreq)''',
    },
    {'longarg': 'fitness=',
     'default': [1, 1.0001, 1.0002],
     'label': 'Fitness of genotype AA,Aa,aa',
     'allowedTypes': [ListType, TupleType],
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
     'validate': valueListOf(valueGE(0.)),
    },
    {'longarg': 'selMultiLocusModel=',
    'default': 'none',
    'label': 'Multi-locus selection model',
    'description': '''Model of overall fitness value given fitness values for each DSL.
                multiplicative: f =  Prod(f_i)
                additive: f = 1-Sum(1-f_i)
                interaction: the intepretation of fitness parameter is different.
                    see fitness.
                Note that selection will be applied to all generations, but backControlledLoci
                will only have wild-type allele before a mutant is introduced.
                ''',
    'allowedTypes': [StringType],
    'chooseOneOf': ['additive', 'multiplicative', 'interaction', 'none']
    },
    #
    {'separator': 'Population expansion'},
    {'longarg': 'expandedName=',
     'default': 'expanded.bin',
     'useDefault': True,
     'description': '''Name of the expanded population, relative to simulation path''',
     'allowedTypes': [StringType],
    },
    {'longarg': 'expandGen=',
     'default': 100,
     'useDefault': True,
     'label': 'Generations to expand',
     'description': '''Number of generations to evolve during the population
                expansion stage''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'expandSize=',
     'default': 42000,
     'useDefault': True,
     'label': 'Expanded population size',
     'description': '''Size of the expanded population. The default value if the recommended
                value when all hapmap populations are used (60+60+90)*200. You may want to
                reduce it according to the population used, or increase it if disease
                prevalence if low and insufficient cases are generated.''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGE(100)
    },
     #
    {'separator': 'Population admixture'},
    {'longarg': 'admixedName=',
     'default': 'admixed.bin',
     'useDefault': True,
     'description': '''Name of the admixed, relative to simulation path''',
     'allowedTypes': [StringType],
    },
    {'longarg': 'migrModel=',
     'default': 'Continuous Gene Flow',
     'useDefault': True,
     'label': 'Migration model',
     'description': '''Migration model used in the admixture stage, It can be 
                'Hybrid Isoloation', 'Continuous Gene Flow', 'Customized' and 'None'. 
                In the HI model, admixture occurs in a single generation and
                is followed by recombination and drift, with no further genetic
                contribution from either parental populations. In the CGF
                model, migration happens with given migration matrix during
                the migration stage. For example, [[0.9, 0.1], [0, 1]] means
                moving 10% from population 1 to 2. 'Custimized' migration model
                allows you to define your own migration model. A function
                migrModel needs to be defined in this script, which returns
                a migration rate at each generation. See the 'migrFunc' function 
                in this script for details. If 'None' is chose, there will be
                no migration. Note that the merge of two populations can be 
                mimiced by a Hybrid Isolation migration of rate [[1, 0], [1, 0]].
                That is to say, everyone from the second subpopulationmoves to the
                first. The three subpopulation case is similar.''',
     'chooseOneOf': ['Hybrid Isolation', 'Continuous Gene Flow', 'Customized', 'None'],
     'allowedTypes': [StringType],
     'validate': valueOneOf(['Hybrid Isolation', 'Continuous Gene Flow', 'Customized', 'None'])
    },
    {'longarg': 'migrGen=', 
     'default': 5,
     'useDefault': True,
     'label': 'Migration generations',
     'description': '''Length of migration stage. If set to zero, the migration stage
                is ignored''',
     'allowedTypes': [IntType, LongType],
     'validate': valueGE(0),
    },
    {'longarg': 'migrRate=',
     'default': [[0.9, 0.1], [0., 1.]],
     'label': 'Migration rate matrix',
     'useDefault': True,
     'description': '''Migration rate matrix. Use only for the continuous gene flow model.
                A_ij of this matrix represents the probability of moving from population i 
                to j, and A_ii is the probability of staying in the same population, and 
                is calculated as 1-sum_(j \\ne i) A_ij. It is possible to create another 
                subpopulation in this way, like sending some individuals from both parental
                populations to a new subpopulation. ''',
     'allowedTypes': [TupleType, ListType],
     'validate': valueListOf(valueListOf(valueBetween(0,1))),    
    },
]


class admixtureParams:
    ''' This class is used to wrap all parameters to a single object so that
    I do not have to pass a bunch of parameters here and there. 
    This class also clean up/validate parameters and calcualtes some derived
    parameters for later uses.
    '''
    def __init__(self, options, allParam):
        # when user click cancel ...
        if len(allParam) == 0:
            sys.exit(1)
        # -h or --help
        if allParam[0]:    
            print usage(options, __doc__)
            sys.exit(0)
        # expand all params to different options
        (self.name, self.useSavedSeed, self.useSavedExpanded,
            self.step, self.showAlleleFreq, self.figureStep,
            self.showLDPlot, self.haploview, self.ldRegionsTmp,
            self.ldSampleSize,
        self.HapMap_dir, self.pops, self.markerList, self.chrom,
        self.numMarkers, self.startPos, self.endingPos, 
        self.minAF, self.minDiffAF, self.minDist, 
            self.seedName, self.initCopy, self.initGen, self.seedSize,
        self.mutaRate, self.recIntensity, self.forCtrlLoci, self.forCtrlFreq,
        self.backCtrlLoci, self.backCtrlFreq, self.fitness, self.mlSelModel,
            self.expandedName, self.expandGen, self.expandSize,
        self.admixedName, self.migrModel, self.migrGen, self.migrRate) \
            = allParam[1:]
        # preparations
        self.createSimulationDir()
        self.saveConfiguration()
        self.seedFile = self.setFile(self.seedName)
        self.expandedFile = self.setFile(self.expandedName)
        self.admixedFile = self.setFile(self.admixedName)
        #
        self.ctrlLoci = self.forCtrlLoci + self.backCtrlLoci
        # adjust parameters startPos, endPos etc.
        self.prepareMarkerParams()

    def createSimulationDir(self):
        '''Create a directory with simulation name'''
        if not os.path.isdir(self.name):
            print 'Creating directory', self.name
            os.makedirs(self.name)
        if not os.path.isdir(self.name):
            raise SystemError('Can not create directory %s, exiting' % self.name)

    def saveConfiguration(self):
        '''Save configuration to $name.cfg'''
        cfgFile = os.path.join(self.name, self.name + '.cfg')
        print 'Save configuration to', cfgFile
        # save current configuration
        saveConfig(options, cfgFile, allParam)

    def setFile(self, filename):
        '''Return $name/filename unless filename is absolute'''
        if os.path.isabs(self.seedName):
            self.seedFile = self.seedName
        else:
            self.seedFile = os.path.join(self.name, self.seedName)

    def setCtrlLociIndex(self, pop):
        '''Translate ctrlLoci to ctrlLociIdx, etc'''
        try:
            self.forCtrlLociIdx = pop.lociByName(self.forCtrlLoci)
            self.backCtrlLociIdx = pop.lociByName(self.backCtrlLoci)
            self.ctrlLociIdx = pop.lociByName(self.ctrlLoci)
        except:
            raise ValueError('''Can not find one of the controlled loci %s in this population
                Please check markers.lst for a list of used markers and their frequency''' % \
                ', '.join(ctrlLoci))
        # this is used for statistical output
        pop.dvars().ctrlLoci = ctrlLoci

    def expandToList(self, par, size, err=''):
        '''If par is a number, return a list of specified size'''
        if type(par) in [IntType, LongType]:
            return [par]*size
        elif type(par) in [TupleType, ListType] and len(par) == 1:
            return list(par)*size
        elif len(par) != size:
            raise ValueError(err)
        else:
            return par

    def prepareMarkerParams(self):
        '''validate marker parameters'''
        if not os.path.isdir(self.HapMap_dir):
            print 'HapMap directory %s does not exist, creating one.' % HapMap_dir
            os.makedirs(self.HapMap_dir)
            if not os.path.isdir(HapMap_dir):
                raise ValueError('Can not create directory %s to store hapmap data, exiting' % HapMap_dir)
        if len(self.chrom) == 0:
            raise ValueError('Please specify one or more chromosomes')
        # in case that chrom is a tuple
        self.chrom = list(self.chrom)
        numChrom = len(self.chrom)
        self.numMarkers = self.expandToList(self.numMarkers, numChrom,
            'Please specify number of marker for each chromosome')
        self.startPos = self.expandToList(self.startPos, numChrom,
            'Wrong starting positions')
        self.endingPos = self.expandToList(self.endingPos, numChrom,
            'Wrong endinging positions')
        self.numMarkers = self.expandToList(self.numMarkers, numChrom,
            'Wrong endinging positions')
        # now, which subpopulations are needed?
        self.popsIdx = []
        for idx,sp in enumerate(HapMap_pops):
            if sp in self.pops:
                print "Using hapmap population %s" % sp
                self.popsIdx.append(idx)
        print "Loading populations ", self.popsIdx
     
#####################################################################
# You have realized how many lines of code is used for parameter
# handling and comments. Now, the real part...
#####################################################################

def expDemoFunc(N0, N1, gen):
    '''
    Return an exponential population expansion demographic function.
    simuUtil::ExponentialExpansion can not actually be used before it expects
        initial and ending total population size.
    
    N0: a list of initial subpopulation sizes.
    N1: ending subpopulation sizes.
    gen: generations to evolve.
    '''
    if type(N1) in [IntType, LongType]:
        NN = [int(N1/len(N0))]*len(N0)
    else:
        NN = N1
    if len(N0) != len(NN):
        raise exceptions.ValueError("Number of subpopulations should be the same")
    rate = [math.log(NN[x]) - math.log(N0[x])/gen for x in range(len(N0))]
    def func(gen, oldSize=[]):
        return [int(N0[x]*math.exp(gen*rate[x])) for x in len(oldSize)]
    return func

    
def printInfo(pop):
    'Print out some information about the population'
    Stat(pop, alleleFreq=range(pop.totNumLoci()))
    # 1. number of population and population size
    print 'Number of populations: %d' % pop.numSubPop()
    print 'Size of populations: %s' % (', '.join([str(x) for x in pop.subPopSizes()]))
    # 2. markers
    print 'Number of chromosomes: %s' % pop.numChrom()
    print 'Number of markers: %s' % (', '.join([str(x) for x in pop.numLoci()]))
    for ch in range(pop.numChrom()):
        print 'Chromosome %d' % ch
        pos = pop.lociPos()[pop.chromBegin(ch):pop.chromEnd(ch)]
        print '    range %.3f - %.3f' % (pos[0], pos[-1])
        dist = [pos[i] - pos[i-1] for i in range(1, len(pos))]
        print '    dist: min %.3f, max %.3f, mean %.3f' % (min(dist), max(dist), sum(dist)/len(dist))
        freq = [pop.dvars().alleleFreq[x][0] for x in range(pop.chromBegin(ch), pop.chromEnd(ch))]
        maf = [min(x, 1-x) for x in freq]
        print '    allele freq: min %.3f, max %.3f, mean %.3f' % (min(maf), max(maf), sum(maf)/len(maf))
        if pop.numSubPop() == 2:
            freq0 = [pop.dvars(0).alleleFreq[x][0] for x in range(pop.chromBegin(ch), pop.chromEnd(ch))]
            freq1 = [pop.dvars(1).alleleFreq[x][0] for x in range(pop.chromBegin(ch), pop.chromEnd(ch))]
            maf0 = [min(x, 1-x) for x in freq0]
            maf1 = [min(x, 1-x) for x in freq1]
            diff = [abs(maf0[x]-maf1[x]) for x in range(len(maf0))]
            print '    diff in allele freq: min %.3f, max %.3f, mean %.3f' % (min(diff), max(diff), sum(diff)/len(diff))


def createInitialPopulation(par):
    '''Create an initial population, with parameters (from the par structure)
    HapMap_dir:     directory that stores hapmap data.
    chrom:          chromosomes to use
    markerList:     list of markers to use
    numMarkers:     number of markers per chromosome
    startPos:       starting position on each chromosome
    endPos:         ending position on each chromosome
    minAF:          minimal allele frequency
    minDiffAF:      minimal allele frequency differences among HapMap populations
    minDist:        minimal distance between adjecent markers
    pops:           hapmap populations to use    
    ''' 
    # load markers!
    for ch in par.chrom:
        if not os.path.isfile(os.path.join(par.HapMap_dir, 'hapmap_%d.bin' % ch)):
            try:
                import loadHapMap
                loadHapMap.loadHapMap([ch], par.HapMap_dir)
            except Exception, e:
                print e
            if not os.path.isfile(os.path.join(par.HapMap_dir, 'hapmap_%d.bin' % ch)):
                raise ValueError('''Failed to load or download hapmap data for chromosome %d
                    Please copy script loadHapMap.py to the current directory, or add
                    path to this script to environmental variable$PYTHONPATH,
                    or run this script manually to download, import, and save HapMap
                    data in simuPOP format''' % ch)
    useHapMapMarker = par.markerList == ''
    ch_pops = []
    if useHapMapMarker:
        for ch, sp, ep, nm in zip(par.chrom, par.startPos, par.endingPos, par.numMarkers):
            ch_pops.append(getMarkersFromRange(par.HapMap_dir, par.popsIdx,
                ch, sp, ep, nm, par.minAF, par.minDiffAF, par.minDist))
        # merge all populations (different chromosomes)
        if len(ch_pops) > 1:
            pop = MergePopulationsByLoci(ch_pops)
        else:
            pop = ch_pops[0].clone()
    else:
        # read the list
        print 'Reading marker list %s' % par.markerList
        mlist = open(par.markerList)
        names = {}
        lastpos = [0]*len(par.chrom)
        for line in mlist.readlines():
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.split()
            ch = int(float(fields[0]))
            pos = float(fields[1])/1000000.
            name = fields[2]
            if ch not in par.chrom:
                continue
            chIdx = par.chrom.index(ch)
            if pos < par.startPos[chIdx]:
                continue
            if par.endingPos[chIdx] != 0 and pos > par.endingPos:
                continue
            if par.minDist > 0 and pos - par.lastpos[chIdx] < par.minDist:
                continue
            if not names.has_key(ch):
                names[ch] = []
            names[ch].append(name)
            lastpos[chIdx] = pos
        pop = getMarkersFromName(par.HapMap_dir, par.names, 
            chroms=par.chrom, hapmap_pops=par.popsIdx, 
            minDiffAF=par.minDiffAF, numMarkers=par.numMarkers)
    # if this population fine?
    if pop.numChrom() != len(par.chrom):
        raise ValueError('Something wrong. The population does not have enough chromosomes')
    return pop


def generateSeedPopulation(par):
    '''Generate seed population, using HapMap dataset and specified marker list'''
    pop = createInitialPopulation(par)
    #
    initPop = pop.clone()
    # print population summary
    print "Initial population has "
    print "    %d subpopulations " % pop.numSubPop()
    print "    %d total number of markers " % pop.totNumLoci()
    for ch in range(pop.numChrom()):
        print "Chromosome %d has %d markers in the range between %.3f and %.3f" \
            % (ch, pop.numLoci(ch), pop.locusPos(pop.chromBegin(ch)),  \
                pop.locusPos(pop.chromEnd(ch)-1))
    # evolve the initial population
    print "Evolving the initial population"
    simu = simulator(pop, randomMating(
        newSubPopSizeFunc = expDemoFunc(pop.subPopSizes(),
            par.seedSize, par.initGen)), rep=1)
    simu.evolve(
        ops = [
            # mutation will be disallowed in the last generation (see later)
            kamMutator(rate=par.mutaRate, loci=range(pop.totNumLoci())),
            recombinator(intensity=par.recIntensity)
        ] + getStatOps(par),
        gen = initGen)
    pop = simu.getPopulation(0, True)
    printInfo(pop)
    #
    # save seed population
    print 'Save seed population to', par.seedFile
    pop.savePopulation(par.seedFile)
    #
    # print marker list fine
    Stat(initPop, alleleFreq=range(0, initPop.totNumLoci()))
    Stat(pop, alleleFreq=range(0, pop.totNumLoci()))
    # write marker information
    print 'Writing a marker list file'
    markers = open(os.path.join(name, 'markers.lst'), 'w')
    print >> markers, 'Name\tchrom\tlocation\t%s\t%s\tseed_freq' % \
        ('\t'.join([x + '_freq' for x in pops]), '\t'.join([x + '_seed_freq' for x in pops]))
    for ch in range(pop.numChrom()):
        for loc in range(pop.chromBegin(ch), pop.chromEnd(ch)):
            print >> markers, '%s\t%d\t%.5f\t%s\t%s\t%.3f' % (pop.locusName(loc), 
                chrom[ch], pop.locusPos(loc), 
                '\t'.join(['%.3f' % initPop.dvars(x).alleleFreq[loc][1] for x in range(initPop.numSubPop())]),
                '\t'.join(['%.3f' % pop.dvars(x).alleleFreq[loc][1] for x in range(pop.numSubPop())]),
                pop.dvars().alleleFreq[loc][1])
    markers.close()
    #
    # write a map file, used by haploview
    print 'Writing a map file ld.map to be used by haploview'
    file = open(os.path.join(name, 'ld.map'), 'w')
    for loc in range(pop.totNumLoci()):
        print >> file, pop.locusName(loc), int(pop.locusPos(loc)*1000000)
    file.close()
    

def migrFunc(gen, curSize):
    ''' return migration rate at each generation'''
    # this is a sample function that migrate to 
    # a third population, with increasing intensity
    return [[0, 0, 0.05*gen], [0, 0, 0.05*gen]]


def getSelector(model, loci, fitness):
    if model in ['additive', 'multiplicative']:
        mlSelModel = {'additive':SEL_Additive, 
            'multiplicative':SEL_Multiplicative}[model]
        return mlSelector(
            # with five multiple-allele selector as parameter
            [ maSelector(locus=loci[x], wildtype=[0], 
                fitness=[fitness[3*x],fitness[3*x+1],fitness[3*x+2]]) \
                    for x in range(len(loci)) ],
            mode=mlSelModel)
    elif model == 'interaction':
        # multi-allele selector can handle multiple DSL case
        return maSelector(loci=loci, fitness=fitness, wildtype=[0])
    else:
        return noneOp()


def getStatOps(par):
    '''Return statistis calculation and progress report operators '''
    # statistics calculation and display
    # 
    ctrlLoci = par.controlledLoci + par.backControlledLoci
    if len(ctrlLoci) > 0 and par.showAlleleFreq:
        statOps = [
            # note: DSL is set when the population is created
            stat(popSize=True, alleleFreq = ctrlLoci, 
                step = step),
            pyEval(r'"gen=%d, size=%s, alleleFreq=%s\n" % (gen, subPopSize,' + \
                r'", ".join([".3f" % alleleFreq[x][1] for x in DSL]))', step=step)
        ],
    else:
        statOps = [
            stat(popSize=True, step = step),
            pyEval(r'"gen=%d, size=%s\n" % (gen, subPopSize)', step=step)
        ]


def getVisualizationOps():
    # ld plot?
    if showLDPlot and figureStep > 0 and ldSampleSize > 0 and len(ldRegionsTmp) > 0:
        ldRegions = []
        if type(ldRegionsTmp[0]) in [IntType, LongType] and len(ldRegionsTmp) == 2:
            ldRegions.append(ldRegionsTmp)
        else:
            ldRegions = ldRegionsTmp
        #
        statOps.append(pyOperator(func=drawLDPlot, 
            param = (haploview, ldRegions, ldSampleSize, pops, name),
            step=figureStep))
        statOps.append(pyOperator(func=drawLDPlot, 
            param = (haploview, ldRegions, ldSampleSize, pops, name),
            at=[-1]))


def drawLDPlot(pop, param):
    '''Draw and display ld plot'''
    (haploview, ldRegions, sampleSize, pops, dir) = param
    sample = RandomSample(pop, [sampleSize]*pop.numSubPop())[0]
    for idx,subPop in enumerate(pops):
        removeSP = range(pop.numSubPop())
        removeSP.remove(idx)
        spSample = sample.clone()
        spSample.removeSubPops(removeSP)
        for reg in ldRegions:
            name = 'LD_%s_%d-%d' % (subPop, reg[0], reg[1])
            filename = os.path.join(dir, name)
            regSample = spSample.clone()
            regSample.removeLoci(range(reg[0], reg[1]))
            SaveMerlinPedFile(pop, output=filename,
                outputAffection=True, affectionCode=['1', '2'])
            os.system('%s -pedfile %s.ped -map %s/ld.map -compressedpng -q -n' % \
                (haploview, filename, dir) )
    return True


def expandSeedPopulation(seedPop, expandGen, expandSize, 
        mutaRate, recIntensity, controlledLoci, controlledFreq, backControlledLoci,
        backControlledFreq, fitness, mlSelModel, step, showAlleleFreq,
        expandedFile):
    '''Expand seed population'''
    # load seed population
    if type(seedPop) == type(''):
        print 'Loading seed population %s...' % seedPop
        pop = LoadPopulation(seedPop)
    else:
        pop = seedPop
    #
    popSizeFunc = expDemoFunc(pop.subPopSizes(), par.expandSize,
        par.expandGen)
    if len(controlledLoci) == 0:
        # evolve it
        #    # evolve the initial population
        #
        print "Evolving the seed population"
        simu = simulator(pop, randomMating(newSubPopSizeFunc = 
            popSizeFunc),rep=1)
        simu.evolve(
            ops = [
                # mutation will be disallowed in the last generation (see later)
                kamMutator(rate = mutaRate, loci=range(pop.totNumLoci())),
                recombinator(intensity = recIntensity),
                getSelector(mlSelModel, controlledLoci + backControlledLoci, fitness),
            ] + getStatOps(step, controlledLoci + backControlledLoci, showAlleleFreq),
            gen = expandGen)
        pop = simu.getPopulation(0, True)
    else:
        # simulate a frequency trajectory
        print "Using controlled random mating on markers %s" % (', '.join(controlledLoci))
        ctrlLoci = pop.lociByNames(controlledLoci)
        Stat(pop, alleleFreq=ctrlLoci)
        currentFreq = []
        # in the order: LOC0: sp0, sp1, sp2, LOC1: sp1, sp2, sp3, ...
        for loc in ctrlLoci:
            print "Current overall frequency %s: %.3f" % (pop.locusName(loc),
                pop.dvars().alleleFreq[loc][1])
            for sp in range(pop.numSubPop()):
                currentFreq.append(pop.dvars(sp).alleleFreq[loc][1])
        print 'Simulating frequency trajectory ...'
        if len(controlledLoci) > 0:
            traj = ForwardFreqTrajectory(
                curGen = 0,
                endGen = expandGen,
                curFreq = currentFreq,
                freq = controlledFreq,
                fitness = fitness,
                NtFunc = popSizeFunc
                )
            intoOps = []
        else:
            # clear these loci
            print 'Clearing mutants at backward-controlled loci'
            backCtrlLoci = pop.lociByNames(backwardControlledLoci)
            for loc in backCtrlLoci:
                for ind in pop.individuals():
                    ind.setAllele(0, loc, 0)
                    ind.setAllele(0, loc, 1)
            # simulate trajectory
            (traj, introGens, trajFunc) = FreqTrajectoryMultiStochWithSubPop(
                curGen=expandGen,
                numLoci=len(backCtrlLoci),
                freq=backControlledFreq,
                NtFunc=popSizeFunc,
                fitness=fitness,
                minMutAge=1, 
                maxMutAge=expandGen, 
                mode='even',
                restartIfFail=True)
            intoOps = [
                pointMutator(atLoci=[backCtrlLoci[i]], toAllele=1, inds=[i],
                    at = [introGens[i]], stage=PreMating)
                for i in range(len(backCtrlLoci))]
        if len(traj) == 0:
            raise ValueError('''Failed to simulate trajectory
                Initial allele frequency:
                Ending allele frequency: %s''' % (currentFreq, controlledFreq))
        # define a trajectory function
        def trajFunc(gen):
            return [x[gen] for x in traj]
        #
        simu = simulator(pop, 
            controlledRandomMating(
                loci=ctrlLoci,
                alleles=[1]*len(ctrlLoci),
                freqFunc=trajFunc,
                newSubPopSizeFunc=popSizeFunc)
        )
        simu.evolve(
            ops =  [
                # mutation will be disallowed in the last generation (see later)
                kamMutator(rate=mutaRate, loci=range(pop.totNumLoci())),
                recombinator(intensity=recIntensity),
            ] + introOps + statOps,
            gen = expandGen
        )
        pop = simu.getPopulation(0, True)
        Stat(pop, alleleFreq=ctrlLoci)
        for i,loc in enumerate(ctrlLoci):
            print "Locus %s: designed freq: (%.3f, %.3f), freq: %.3f" % \
                (pop.alleleName(loc), controlledFreq[i][0],
                controlledFreq[i][1], pop.dvars().alleleFreq[loc][1])
    print 'Saving expanded population to ', expandedFile
    pop.savePopulation(expandedFile)
    return pop

 
def mixExpandedPopulation(pop, migrModel, migrGen, migrRate, mutaRate,
    recIntensity, selLoci, selModel, fitness, admixedFile):
    ''' Evolve the seed population
    '''
    # migration part.
    mergeAt = 1000000  # default not merge
    if migrModel == 'Hybrid Isolation':
        print 'Using %s model' % migrModel
        migr = noneOp()
        mergeAt = 0
    elif migrModel == 'None':
        print 'Do not migrate'
        migr = noneOp()
    elif migrModel == 'Continuous Gene Flow':
        print 'Using %s with migration rate %s' % (migrModel, migrRate)
        migr = migrator(rate=migrRate, mode=MigrByProbability)
    elif migrModel == 'Customized':
        print 'Using customized migration model'
        migr = pyMigrator(rateFunc=migrFunc, mode=MigrByProbability)
    if migrGen > 0:
        simu = simulator(pop, randomMating())
        simu.evolve(
            ops =  [
                # mutation will be disallowed in the last generation (see later)
                kamMutator(rate = mutaRate, loci=range(pop.totNumLoci())),
                recombinator(intensity = recIntensity),
                stat(popSize=True, step=10, begin=9),
                migr,
                getSelector(selLoci, selModel, fitness),
                pyEval(r'"gen=%d, size=%s\n" % (gen, subPopSize)')
            ],
            gen = migrGen
        )
        pop = simu.getPopulation(0, True)
    # save this population
    print "Calculating allele frequency..."
    pop.vars().clear()
    Stat(pop, alleleFreq=range(pop.totNumLoci()))
    print 'Saving admixed population to ', admixedFile
    pop.savePopulation(admixedFile)
    return pop


short_desc = '''This program simulates an admixed population based on 
two or more HapMap populations. Please follow the intructions
of the help message to prepare HapMap population.'''

# determine which script to run.
if __name__ == '__main__':
    # 
    # PARAMETER HANDLING
    # 
    # get all parameters
    allParam = getParam(options, short_desc, __doc__, nCol=2)
    par = admixtureParams(options, allParam)
    #
    # SEED POPULATION GENERATION
    #
    # get seed population
    if par.useSavedExpanded or (par.useSavedSeed and os.path.isfile(par.seedFile)):
        if par.useSavedSeed:
            print "Using existing seed file ", par.seedFile
        seedPop = None
    else:
        seedPop = generateSeedPopulation(par)
    sys.exit(0)
    #
    # step 2:
    # 
    # if both files exists, skip this stage
    if useSavedExpanded and os.path.isfile(expandedFile):
        expandedPop = None
    else:
        if seedPop is None:
            print 'Loading seed population', seedFile
            seedPop = LoadPopulation(seedFile)
        # 
        if len(controlledLoci) != 0 and len(backwardControlledLoci) != 0:
            raise ValueError('This script currently only allows one kind of controlled loci' + \
                'Please specify only one of --controlledLoci and --backwardControlledLoci')
        #
        # fitness
        numDSL = len(controlledLoci) + len(backControlledLoci)
        if mlSelModel == 'none':
            fitness = []
        elif mlSelModel == 'interaction':
            if numDSL == 1:
                raise ValueError("Interaction model can only be used with more than one DSL");
            if len(fitnessTmp) != 3**numDSL:
                raise ValueError("Please specify 3^n fitness values for n DSL");
            fitness = fitnessTmp
        else:
            if fitnessTmp == []:    # neutral process
                fitness = [1,1,1]*numDSL
            else:
                # for a single DSL
                if len(fitnessTmp) == 3:
                    fitness = fitnessTmp*numDSL
                elif len(fitnessTmp) != numDSL*3:
                    raise ValueError("Please specify fitness for each DSL")
                else:
                    fitness = fitnessTmp
        #
        controlledFreq = []
        if len(controlledFreqTmp) > 0:
            if type(controlledFreqTmp[0]) in [TupleType, ListType]:
                if len(controlledFreqTmp) != len(controlledLoci):
                    raise ValueError('Please specify frequency range for each controlled locus')
                for rng in controlledFreqTmp:
                    if len(rng) != 2:
                        print "Wrong allele frequency range", rng
                    controlledFreq.append(rng)
            # give only one
            else:
                if len(controlledFreqTmp) != 2:
                    print "Wrong allele frequency range", controlledFreqTmp
                for i in range(len(controlledLoci)):
                    controlledFreq.append(controlledFreqTmp)
        # backward controlled freq
        backControlledFreq = []
        if len(backControlledFreqTmp) == 1:
            for lpc in backControlledLoci:
                backControlledFreq.append(backControlledFreqTmp)
        elif len(backControlledFreqTmp) != len(backControlledLoci):
            raise ValueError('Number of backward controlled freq does not match the number of such loci')
        else:
            backControlledFreq = backControlledFreqTmp
        #
        expandedPop = expandSeedPopulation(par)
    # admixture
    if expandedPop is None:
        print 'Loading expanded population from file ', expandedFile
        expandedPop = LoadPopulation(expandedFile)
    admixedPop = mixExpandedPopulation(par)
  
