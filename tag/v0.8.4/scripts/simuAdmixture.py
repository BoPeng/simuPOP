#!/usr/bin/env python
'''
This script simulates an admixed population based on the HapMap dataset. Please
read this help message carefully, making sure your know how this script works,
then run a few test commands, before you explore the capacity of this script.


Step 0: Prepare HapMap dataset. (call scripts/loadHapMap.py if available)
=========================================================================

This script makes use of the HapMap dataset. The dataset is downloaded, imported
and saved in simuPOP format automatically, using script scripts/loadHapMap.py.
If loadHapMap.py can not be imported (not in the working directory or in $PYTHONPATH),
please try to run loadHapMap.py manually and provide a path with files hapmap_XX.bin
to parameter --HapMap_dir.


Step 1: Generate a seed population (scripts/simuAdmixture.py)
==============================================================

Determine which markers to use.
-----------------------------------

Two methods are provided.

a) Use a marker list in the format of 'chromosome position name'. In
addition, you can specify which chromosome(s) to use, starting and ending
position(s), number of markers on each chromosome, minimal distance
between adjacent markers, and minimal allele frequency difference between
two HapMap populations. 'position' in the marker list is assumed to be
in base pair (as in Affymetrix or Illumina annotation data). Marker positions
are converted to centiMorgan by dividing the positions by 10^6 (i.e. 
1 centiMorgan = 1 M basepairs). Starting and ending position should be
inputted in centiMorgan.

b) Use selected HapMap markers. In addition to which chromosome(s) to use,
starting and ending position(s), number of markers on each chromosome, 
minimal distance between adjacent markers, you can specify minimal minor
allele frequency, and minimal allele frequency difference between HapMap
populations.

The first method is suitable when you would like to simulate a population
with markers in a real dataset, such as the Affymetrix 100k chipset.
The 'allele frequency difference' criteria can be used to find markers
of high ancestray information content.


Generate a seed population. 
----------------------------

Two or three HapMap populations can be used. The population is expanded
instantly by copying individual 10 (default value of parameter initCopy)
time to avoid quick loss of heterogenity due to small population sizes. 
and evolve for a relatively long period of time (200 generations by default)
without migration. This stage is designed to make HapMap populations 
more distinct, because the HapMap populations are already admixed.

During the evolution, recombination at a rate of 0.01 per cM, and mutation
at a rate of 1e-7 per nucleotide per generation is applied. Recombination
and mutation will also be applied at the subsequent stages.

This population will be saved as 'seed.bin' (configurable through --seed
parameter) that will be used for the subsequent replicate simulations. We
consider the subpopulations of the seed population as the populations around
2000 years ago, that have been evolved separately from the Out-Of-Africa 
population for abut 2000 generations

Each simulation will be given a name and all files will be saved to a directory
named after it. This step will by default produce file $name/seed.bin.


Step 2. Evolve the seed population.
=====================================

After a seed population is generated (or loaded if already exists and parameter
--useSavedSeed is given), this script will evolve it as follows:

1. Evolve the seed population subject to rapid population expansion for
100 (default value of expandGen) generations. This is to mimic the rapid
population expansion in the past 2000 years. Exponential population 
growth model is used. No migration is allowed.

2. Mix the populations using given migration model. Two basic models,
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
  10% of individuals will migrate from population 2 to 1 at each 
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


The result of this stage will be saved to two files
  $name/expanded.bin (after population expandion)
  $name/admixed.bin (after population admixture)
The script will skip population expansion and/or admixture if the corresponding
file is found and parameter useSavedExpanded and/or useSavedAdmixed is
specified. In the useSavedExpanded case, the simulation starts mixing population
$name/expanded.bin directly. In the useSavedAdmixed case, the simulation draw
sample from the admixed population directly.


Step 3: Sample from admixed population.
========================================

This stage does not really belong to this script because the goal of the
script is to simulate admixed populations. Arbitrary quantitative trait,
penetrance models can be applied to the simulated population, and arbitrary
ascertain scheme can then be applied to the population. Ideally, one should
load $name/admixed.bin and draw sample from it. However, for the sake of
convenience and completeness, this script allows three kinds of quantitative
trait and sample models.


A multi-locus quantitative trait model
---------------------------------------

Using this method, a list of chromosomes with disease pre-disposing locus
is given (if there are more than one disease locus, repeat the chromosome
number), along with desired allele frequency (and allowed window around this
frequency). From given positions (default to 0.5, from the middle of 
each chromosome), this script will locate loci that satify allele frequency
requirement.

Then, a list of 'percentage of variation that is explained by each disease
locus' is required. An additive quantitative trait model is used which 
yields trait value -a, 0 and a for genotype AA, Aa and aa respectively.
Assume allele frequency p, a locus would contribute 2p(1-p)a^2 variance to
the total trait variance. With given variance v, a is calculated using
sqrt(v/(2p(1-p))).

The rest of the variance is simulated using a normal distribution.


A trait-derived affection status model
---------------------------------------

Using a above quantitative trait model, a cutoff value can be given which
assign affection status to individuals according to their trait value.


With a given penetrance model
------------------------------

It is also possible to use an explicit penetrance model to assign affection
status. Using this approach, a list of markers (by names) is expected and 
a penetrance matrix is given. For a two-locus model, the penetrance matrix
should be given as p_AABB, p_AABb, p_AAbb, p_AaBB, p_AaBb, p_Aabb, p_aaBB,
p_aaBb, p_aabb, which p_XXXX is the probability of affected for genotype
XXXX.


Draw a Case-control sample
------------------------------

If affection status is assigned to each individual, case-control samples
can be drawn from the resulting population. One can simply specify number
of cases and controls, or a list of cases and controls from each population.
For example, [0,0,500], [0,0,500] would ignore the first two populations
and draw 500 cases and 500 controls.


Draw a random sample
---------------------

If a quantitative trait is specified, a random sample can be drawn from
the admixed population. Given number of random individuals can be drawn from
from all populations (a single number is given), or from each population
(a list is given).

The samples are saved in Merlin QTDT format, which can be easily converted
to other formats.

If parameter --resample is given, quantitative trait or penetrance models are 
applied to $name/admixed.bin.

The following are a few examples of using this script. If you use a parameter
dialog, you can get the command line from the last line of the configuration 
file.


Test scripts
==============

The following test scripts demonstrate the use of this script using a small
number of loci. Note that most parameters have default parameters so the
command can be shortened considerably in practise. I list (almost) the full
command for clarity purposes.

simuAdmixture.py --noDialog  --reseed --HapMap_dir='../../HapMap' \
    --chrom="range(1,3)"  --markerList='../../Affy/mapAffySNPs.txt' \
    --startPos="[0]" --endingPos='[0]' --numMarkers="[100,100]" --minAF='0' --minDist='0'  \
    --pops="['CEU', 'YRI', 'JPT+CHB']" --initCopy='10' --gen='20' --size='3600'  \
    --seed='test_seed.bin' --saveConfig='test_seed.cfg'

simuAdmixture.py --noDialog  --seed=test_seed.bin --migrGen='5' \
    --migrRate="([1, 0, 0], [1, 0, 0], [1, 0, 0])" --remix
    
simuAdmixture.py --noDialog  --HapMap_dir='../../HapMap' --chromWithDSL="[1,2]" \
    --sampleSize="(5, 5)" --freqDSL=0.1 --cutoff=1 \
    --dslVar="[0.05, 0.1]"   --seed='test_seed.bin'


simulation for XJ Gu et al (2008)
====================================

# 1: Hybrid isolation.
simuAdmixture.py --noDialog  --name='IH' --initPop='../../Affy/affyAll_CEU.bin' \
    --HapMap_dir='../../HapMap' --pops="['CEU']" --chrom="range(1, 23)" \
    --markerList='../../Affy/mapAffySNPs.txt'  --numMarkers="[0]*22" \
    --startPos='[0]' --endingPos='[0]' --minAF='0' --minDiffAF='0' --minDist='0' --initCopy='10' \
    --gen='200' --size='4800' --expandGen='100' --expandSize='24000' \
    --useSavedAdmixed=True --migrModel='None' --migrGen='1' \
    --migrRate="()" --chromWithDSL="(1, 2, 3, 4)" \
    --freqDSL='0.2' --freqDev='0.02' --dslVar="(0.005, 0.01, 0.03, 0.05)" \
    --cutoff='-0.5' --DSLpene='[]' --peneFunc='None' --parameter='[0.5]' --ccSampleSize="(600, 600)" \
    --ccSampleName='case-control' --randomSampleSize='800' --randomSampleName='random'

simuAdmixture.py --noDialog  --name='admix' --useSavedSeed=True --initPop='' --HapMap_dir='../../HapMap' \
    --pops="['CEU', 'YRI', 'JPT+CHB']" --chrom="(2, 3)" --markerList='' --numMarkers="(1000, 1000)" \
    --startPos='[0]' --endingPos='[0]' --minAF='0.05' --minDiffAF='0' --minDist='0' --initCopy='10' \
    --gen='200' --size='4800' --useSavedExpanded=True --expandGen='100' --expandSize='24000' \
    --useSavedAdmixed=True --migrModel='Continuous Gene Flow' --migrGen='5' \
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

simuAdmixture.py --noDialog  --HapMap_dir='../HapMap' --chrom='[2]' \
   --numMarkers=500 --startPos=100  --minAF=0.1 --minDiffAF=0.2 \
   --minDist=0.05 --pops="['CEU', 'YRI']" --reseed

250 sample from CEU, 250 sample from YRI, 500 sample from admixed population
when 10% of the CEU popopulation migrate to YRI for 5 generations.

Round 1: expand and get sample from CEU

simuAdmixture.py  --noDialog --migrModel='None' --migrGen='0' --sampleType='random'\
  --sampleSize="(250, 0)" --sampleName='CEU' --name='simu2'

Round 2: load expanded population and get sample from YRI

simuAdmixture.py  --noDialog --remix=True --migrModel='None' --migrGen='0' --sampleType='random'\
  --sampleSize="(0, 250)" --sampleName='YRI' --name='simu2'

Bo Peng (bpeng@mdanderson.org)

$LastChangedDate: 2007-10-12 14:48:57 -0500 (Fri, 12 Oct 2007) $
$Rev: 1090 $

'''

from simuOpt import *
setOptions(alleleType='binary')
from simuPOP import *
from hapMapUtil import getMarkersFromName, getMarkersFromRange, evolveHapMap

import os, sys, types, exceptions, getopt, math, random
from simuUtil import SaveQTDT

HapMap_pops = ['CEU', 'YRI', 'JPT+CHB']

options = [
    {'arg': 'h',
     'longarg': 'help',
     'default': False, 
     'description': 'Print this usage message.',
     'allowedTypes': [types.NoneType, type(True)],
     'jump': -1                    # if -h is specified, ignore any other parameters.
    },
    {'longarg': 'name=',
     'default': 'simu',
     'useDefault': True,
     'allowedTypes': [types.StringType],
     'label': 'Name of the simulation',
     'description': '''Name of this simulation. A directory with this name
                will be created. Configuration file (.cfg) and samples will be
                saved to this directory''',
    },
    {'separator': 'Generate seed population'},
    {'arg': 's:',
     'longarg': 'seed=',
     'default': 'seed.bin',
     'useDefault': True,
     'description': '''Name of the seed population''',
     'allowedTypes': [types.StringType],
    },
    {'longarg': 'useSavedSeed=',
     'default': False,
     'label': 'Use saved seed population',
     'useDefault': True,
     'allowedTypes': [types.StringType, types.BooleanType],
     'description': '''Use specified or a default seed population, if available.
        The default seed population is seed.bin under the simulation directory.
        '''
    },
    {'longarg': 'initPop=',
     'default': '',
     'useDefault': True,
     'label': 'Initial population',
     'allowedTypes': [types.NoneType, types.StringType],
     'description': '''If you already have an initial population, load this one instead
                of loading from HapMap. Setting of this parameter will skip parameters
                chrom, markerList, numMarkers etc.'''
    },       
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
     'allowedTypes': [types.StringType],
     #'validate': valueValidDir(),
    },
    {'longarg': 'pops=',
     'default' : ['CEU', 'YRI', 'JPT+CHB'],
     'useDefault': True,
     'label' : 'HapMap populations',
     'description': '''Which HapMap populations to use?''',
     'allowedTypes': [types.TupleType, types.ListType],
     'chooseFrom': HapMap_pops,
     'validate': valueListOf(valueOneOf(HapMap_pops)),
    },
    {'longarg': 'chrom=',
     'default': [2, 3],
     'label': 'Chromosomes to use',
     'description': '''A list of chromosomes to use from the HapMap data. When multiple
                chromosomes are involves, numMarkers, if used, should be a list that specicy
                number of markers on each chromosome. The same rule applies to startPos
                and endingPos as well.''',
     'allowedTypes': [types.TupleType, types.ListType],
     'validate': valueListOf(valueBetween(1, 22)),
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
                compare marker location.''',
     'allowedTypes': [types.StringType],
     'validate': valueOr(valueEqual(''), valueValidFile()),
    },
    {'longarg': 'numMarkers=',
     'default': [1000, 1000],
     'label': 'Number of markers to use',
     'description': '''Number of markers to use from the marker list file. If 0 is used,
                all markers that satisfy conditions startPos, endingPos, minDist will
                be used.''',
     'allowedTypes': [types.TupleType, types.ListType],
     'validate': valueOr(valueGT(0), valueListOf(valueGE(0)))
    },
    {'longarg': 'startPos=',
     'default': 0,
     'useDefault': True,
     'label': 'staring position',
     'description': '''Starting position of the markers. If multiple
                chromosomes are used, the positions for each 
                chromosome can be specified as a list.''',
     'allowedTypes': [types.TupleType, types.ListType],
     'validate': valueOr(valueGE(0), valueListOf(valueGE(0)))
    },
    {'longarg': 'endingPos=',
     'default': 0,
     'useDefault': True,
     'label': 'Ending position',
     'description': '''Ending position of the markers. Ignored if its value 
                is 0.  If multiple chromosomes are used, the positions for each 
                chromosome can be specified as a list. ''',
     'allowedTypes': [types.TupleType, types.ListType],
     'validate': valueOr(valueGE(0), valueListOf(valueGE(0)))
    },
    {'longarg': 'minAF=',
     'default': 0.05,
     'useDefault': True,
     'label': 'Minimal allele frequency',
     'description': '''Minimal allele frequency, only used for picking markers
                from the HapMap dataset''',
     'allowedTypes': [types.IntType, types.LongType, types.FloatType],
     'validate': valueGE(0)
    },
    {'longarg': 'minDiffAF=',
     'default': 0,
     'useDefault': True,
     'label': 'Minimal allele frequency difference',
     'description': '''Minimal allele frequency difference between two HapMap population,
                , can only be used when two HapMap populations are used. This options can be
                used to choose markers with high ancestry information content.''',
     'allowedTypes': [types.IntType, types.LongType, types.FloatType],
     'validate': valueGE(0)
    },
    {'longarg': 'minDist=',
     'default': 0,
     'useDefault': True,
     'label': 'Minimal distance between markers (cM)',
     'allowedTypes': [types.IntType, types.LongType, types.FloatType],
     'description': '''Minimal distance between markers (in the unit of cM).
                Can be used for both methods.''',
    },
    {'longarg': 'initCopy=',
     'default': 10,
     'useDefault': True,
     'label': 'Initial propagation',
     'description': '''How to expand the initial small HapMap sample to
                 avoid quick loss of heterogenity''',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'gen=',
     'default': 200,
     'useDefault': True,
     'label': 'Generations to evolve',
     'description': '''Number of generations to evolve to get the seed
                population.''',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'size=',
     'default': 4800,
     'useDefault': True,
     'label': 'Ending population size',
     'description': 'Size of the seed population',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': valueGE(100)
    },
    #
    {'separator': 'Population expansion'},
    {'longarg': 'useSavedExpanded=',
     'default': False,
     'useDefault': True,
     'allowedTypes': [types.BooleanType, types.StringType],
     'label': 'Use saved expanded population',
     'description': '''If set to true, load specified or saved $name/expanded.bin and 
                skip population expansion'''
    },
    {'longarg': 'expandGen=',
     'default': 100,
     'useDefault': True,
     'label': 'Expansion generations',
     'description': '''Number of generations to evolve during the population
                expansion stage''',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': valueGT(0)
    },
    {'longarg': 'expandSize=',
     'default': 24000,
     'useDefault': True,
     'label': 'Expanded population size',
     'description': 'Size of the expanded population',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': valueGE(100)
    },
    {'separator': 'Population admixture'},
    {'longarg': 'useSavedAdmixed=',
     'default': False,
     'useDefault': True,
     'allowedTypes': [types.StringType, types.BooleanType],
     'label': 'Use saved admixed population',
     'description': '''If set to true, the program will use saved $name/admixed.bin''',
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
     'allowedTypes': [types.StringType],
     'validate': valueOneOf(['Hybrid Isolation', 'Continuous Gene Flow', 'Customized', 'None'])
    },
    {'longarg': 'migrGen=', 
     'default': 5,
     'useDefault': True,
     'label': 'Migration generations',
     'description': '''Length of migration stage. If set to zero, the migration stage
                is ignored''',
     'allowedTypes': [types.IntType, types.LongType],
     'validate': valueGE(0),
    },
    {'longarg': 'migrRate=',
     'default': [[0.9, 0.1], [0., 1.]],
     'label': 'Migration rate matrix',
     'useDefault': True,
     'description': '''Migration rate matrix. Use only for the continuous gene flow model.
                A_ij of this matrix represents the probability of moving from population i 
                to j, and A_ii is the probability of staying in the same population, and 
                is calculated as 1-sum_(j \ne i) A_ij. It is possible to create another 
                subpopulation in this way, like sending some individuals from both parental
                populations to a new subpopulation. ''',
     'allowedTypes': [types.TupleType, types.ListType],
     'validate': valueListOf(valueListOf(valueBetween(0,1))),    
    },
    {'separator': 'Quantitative trait model'},
    {'longarg': 'DSLtrait=',
     'default': [],
     'label': 'Disease susceptibility loci',
     'useDefault': True,
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Names of disease susceptibility loci for the quantitative trait model.
                If given, parameters chromWithDSL, freqDSL and freqDev will be ignored. 
                Otherwise, these DSL will be determined by allele frequency automatically.'''
    },
    {'longarg': 'chromWithDSL=',
     'default': [1,2],
     'label': 'Chromosomes with DSL',
     'useDefault': True,
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''Chromosomes with DSL. The chromosomes are counted from the loaded
                population and indexed from 1. For example, if chromosomes 5, 6, 8 are
                used in the simulation, [1, 3] here would refer to chromosome 5 and 8.'''
    },
    {'longarg': 'freqDSL=',
     'default': 0.1,
     'label': 'MAF of DSL',
     'useDefault': True,
     'allowedTypes': [types.FloatType],
    },
    {'longarg': 'freqDev=',
     'default': 0.01,
     'useDefault': True,
     'label': 'Allowed deviation from dslFreq',
     'allowedTypes': [types.FloatType],
    },
    {'longarg': 'dslVar=',
     'default': [0.05, 0.1],
     'label': 'proportion of Variance to explain',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Proption of variance explained by each DSL,
                The a in the additive formula is determined by
                this variable and allale frequency'''
    },
    {'longarg': 'cutoff=',
     'default': 0,
     'useDefault': True,
     'label': 'Cutoff value to determine affection status',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': 'Cutoff value used to determine affection status'
    },      
    {'separator': 'Penetrance model'},
    {'longarg': 'DSLpene=',
     'label': 'Disease susceptibility loci',
     'default': [],
     'useDefault': True,
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''A list of markers (by name) that will be used to 
                determine affection status'''
    },
    {'longarg': 'peneFunc=',
     'default': 'None',
     'label': 'Penetrance function',
     'allowedTypes': [types.StringType],
     'description': ''' Penetrance functions to be applied to the final
        population. Two penetrance fucntions are provided, namely recessive
        or additive single-locus model with heterogeneity multi-locus model. 
        You can define another customized penetrance functions by modifying
        this script. ''',
     'validate':    valueOneOf(['recessive', 'additive', 'custom', 'None']),
     'chooseOneOf': [ 'recessive', 'additive', 'custom', 'None']
    },
    {'longarg': 'peneParam=',
     'default': [0.5],
     'label': 'Penetrance parameters',
     'description': '''Penetrance parameter for all DSL. An array of parameter 
        can be given to each DSL. The meaning of this parameter differ by 
        penetrance model. For a recessive model, the penetrance is 0,0,p for 
        genotype AA,Aa,aa (a is disease allele) respectively. For an additive 
        model, the penetrance is 0,p/2,p respectively.''',
     'allowedTypes': [types.ListType, types.TupleType],
     'validate':   valueOneOf([ 
             valueBetween(0,1), valueListOf(valueBetween(0,1))] )
    },
    {'separator': 'Draw samples'},
    {'longarg': 'ccSampleSize=',
     'default': [100,100],
     'label': 'Case-control sample size',
     'description': '''Sample size. If sampleType == 'random' and a single number s
                is given, a random sample of size s will be drawn from the whole
                population. If a list [s1, s2, ...] is given, s1 random individuals will
                be drawn from the first population, s2 random individuals will be
                drawn from the second population, and so on.
                If sampleType == 'case-control' and a two-number list 
                [s1, s2] is given, s1 case and s2 controls will be drawn from the
                whole population. If [list1, list2] is given, list1 is the number
                of cases from each subpopulation, and list2 is the number of controls
                from each subpopulation. The second form obviously can not be
                used for the Hybrid Isolation model because there is no subpopulation
                in the final population. Note that elements in list1 and list2 can
                be zero so you can sample from only one of the subpopulations.''',
      'allowedTypes': [types.IntType, types.LongType, types.TupleType, types.ListType],
      'validate': valueOr(valueGE(0), valueListOf(valueOr(valueGE(0), valueListOf(valueGE(0)))))
    },
    {'longarg': 'ccSampleName=',
     'default': 'case-control',
     'useDefault': True,
     'allowedTypes': [types.StringType],
     'label': 'Name of the case-control sample',
     'description': 'Name of the case-control sample',
    },
    {'longarg': 'randomSampleSize=',
     'default': 100,
     'label': 'Random sample size',
     'description': '''Size of a random sample''',
     'allowedTypes': [types.IntType, types.LongType, types.TupleType, types.ListType],
     'validate': valueOr(valueGE(0), valueListOf(valueOr(valueGE(0), valueListOf(valueGE(0)))))
    },
    {'longarg': 'randomSampleName=',
     'default': 'random',
     'useDefault': True,
     'allowedTypes': [types.StringType],
     'label': 'Name of the random sample',
     'description': 'Name of the random sample',
    },
]



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


def createInitialPopulation(HapMap_dir, chrom, markerList, numMarkers, startPos,
            endingPos, minAF, minDiffAF, minDist, pops):
    '''Create an initial population''' 
    # process parameters
    if not os.path.isdir(HapMap_dir):
        print 'HapMap directory %s does not exist, creating one.' % HapMap_dir
        os.makedirs(HapMap_dir)
    if not os.path.isdir(HapMap_dir):
        print 'Can not create directory %s to store hapmap data, exiting' % HapMap_dir
        sys.exit(1)
    if len(chrom) == 0:
        print 'Please specify one or more chromosomes'
        sys.exit(1)
    for ch in chrom:
        if not os.path.isfile(os.path.join(HapMap_dir, 'hapmap_%d.bin' % ch)):
            try:
                import loadHapMap
                loadHapMap.loadHapMap([ch], HapMap_dir)
            except Exception, e:
                print e
            if not os.path.isfile(os.path.join(HapMap_dir, 'hapmap_%d.bin' % ch)):
                print 'Failed to load or download hapmap data for chromosome %d' % ch
                print 'Please copy script loadHapMap.py to the current directory, or add'
                print 'path to this script to environmental variable$PYTHONPATH,'
                print 'or run this script manually to download, import, and save HapMap'
                print 'data in simuPOP format'
                sys.exit(1)
    # in case that chrom is a tuple
    chrom = list(chrom)
    useHapMapMarker = markerList == ''
    if type(numMarkers) not in [types.TupleType, types.ListType] or \
        len(numMarkers) != len(chrom):
        print "Please specify number of marker for each chromosome: ", numMarkers
        sys.exit(1)
    if type(startPos) not in [types.TupleType, types.ListType] or \
        len(startPos) not in [1, len(chrom)]:
        print "Wrong starting positions"
        sys.exit(1)
    if len(startPos) == 1:
        startPos = startPos * len(chrom)
    if type(endingPos) not in [types.TupleType, types.ListType] or \
        len(endingPos) not in [1, len(chrom)]:
        print "Wrong endinging positions"
        sys.exit(1)
    if len(endingPos) == 1:
        endingPos = endingPos * len(chrom)
    if type(numMarkers) not in [types.TupleType, types.ListType] or \
        len(numMarkers) not in [1, len(chrom)]:
        print "Wrong endinging positions"
        sys.exit(1)
    if len(numMarkers) == 1:
        numMarkers = numMarkers * len(chrom)
    # now, which subpopulations are needed?
    load_sp = []
    for idx,sp in enumerate(HapMap_pops):
        if sp in pops:
            load_sp.append(idx)
    # load markers!
    ch_pops = []
    if useHapMapMarker:
        for ch, sp, ep, nm in zip(chrom, startPos, endingPos, numMarkers):
            ch_pops.append(getMarkersFromRange(HapMap_dir, load_sp, ch, sp, 
                ep, nm, minAF, minDiffAF, minDist))
        # merge all populations (different chromosomes)
        if len(ch_pops) > 1:
            pop = MergePopulationsByLoci(ch_pops)
        else:
            pop = ch_pops[0].clone()
    else:
        # read the list
        print 'Reading marker list %s' % markerList
        mlist = open(markerList)
        names = {}
        lastpos = [0]*len(chrom)
        for line in mlist.readlines():
            (ch, pos, name) = line.split()
            ch = int(float(ch))
            pos = float(pos)/1000000.
            if ch not in chrom:
                continue
            chIdx = chrom.index(ch)
            if pos < startPos[chIdx]:
                continue
            if endingPos[chIdx] != 0 and pos > endingPos:
                continue
            if minDist > 0 and pos - lastpos[chIdx] < minDist:
                continue
            if not names.has_key(ch):
                names[ch] = []
            names[ch].append(name)
            lastpos[chIdx] = pos
        pop = getMarkersFromName(HapMap_dir, names, chroms=chrom,
            hapmap_pops=load_sp, minDiffAF=minDiffAF, numMarkers=numMarkers)
    # if this population fine?
    if pop.numChrom() != len(chrom):
        print "Something wrong. The population does not have enough chromosomes"
        sys.exit(1)
    return pop


def generateSeedPopulation(HapMap_dir, initPop, chrom, markerList, numMarkers, startPos,
            endingPos, minAF, minDiffAF, minDist, pops, initCopy, gen, size, popName):
    'Generate seed population, using HapMap dataset and specified marker list'
    if os.path.isfile(initPop):
        pop = LoadPopulation(initPop)
    else:
        pop = createInitialPopulation(HapMap_dir, chrom, markerList, numMarkers, startPos,
            endingPos, minAF, minDiffAF, minDist, pops)
        if initPop is not None and initPop != '':
            pop.savePopulation(initPop)
    # print population summary
    print "Initial population has "
    print "    %d subpopulations " % pop.numSubPop()
    print "    %d total number of markers " % pop.totNumLoci()
    for ch in range(len(chrom)):
        print "Chromosome %d has %d markers in the range between %.3f and %.3f" \
            % (ch, pop.numLoci(ch), pop.locusPos(pop.chromBegin(ch)),  \
                pop.locusPos(pop.chromEnd(ch)-1))
    #
    # evolve the HapMap population
    pop = evolveHapMap(pop, 
        initMultiple=initCopy,
        gen=gen,
        endingSize=size,
        step=10
    )
    printInfo(pop)
    # save seed population
    print 'Save seed population to', popName
    pop.savePopulation(popName)


def migrFunc(gen, curSize):
    ''' return migration rate at each generation'''
    # this is a sample function that migrate to 
    # a third population, with increasing intensity
    return [[0, 0, 0.05*gen], [0, 0, 0.05*gen]]


def expandSeedPopulation(seedPop, expandGen, expandSize,
        expandedFile):
    '''Expand seed population'''
    # load seed population
    if type(seedPop) == type(''):
        print 'Loading seed population %s...' % seedPop
        pop = LoadPopulation(seedPop)
    else:
        pop = seedPop
    # evolve it
    pop = evolveHapMap(pop,
        gen=expandGen, 
        endingSize=expandSize,
        step=10,
    )
    print 'Saving expanded population to ', expandedFile
    pop.savePopulation(expandedFile)
    return pop

 
def mixExpandedPopulation(pop, migrModel, migrGen, migrRate, admixedFile):
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
        pop = evolveHapMap(pop,
            mergeAt=mergeAt,
            gen=migrGen,
            # constant population size
            endingSize=pop.popSize(),
            step=1,
            migr=migr)
    # save this population
    print "Calculating allele frequency..."
    pop.vars().clear()
    Stat(pop, alleleFreq=range(pop.totNumLoci()))
    print 'Saving admixed population to ', admixedFile
    pop.savePopulation(admixedFile)
    return pop


def setQuanTrait(pop, DSLtrait, chromWithDSL, p, sd, vars, cutoff, name):
    '''Set quantitative trait and affection status if a cutoff value
    is given
    
    DSLtrait: names of each disease susceptibility loci. If given,
        ignore parameters chromWithDSL, p, sd.
    
    chromWithDSL: chromosomes with DSL, chromosomes should be
        indexed from 0

    p: target allele frequency of the DSL

    sd: allowed deviation from p

    vars: variance for each marker

    cutoff: set an individual as 'affected' if its trait value is greater
        than this cutoff value.
    '''
    #
    if len(DSLtrait) == 0:
        numDSL = len(chromWithDSL)
        DSL = []
        sign = []
        for ch1 in chromWithDSL:
            ch = ch1 - 1
            DSL.append(pop.chromBegin(ch) + pop.numLoci(ch)/2)
            for i in range(DSL[-1], pop.chromEnd(ch)):
                af = 1 - pop.dvars().alleleFreq[i][0]
                # allele 1 is the minor allele
                if af > p - sd and af < p + sd:
                    DSL[-1] = i
                    sign.append(1)
                    break
                # allele 0 is the minor allele
                elif 1 - af > p - sd and 1 - af < p + sd:
                    DSL[-1] = i
                    sign.append(-1)
                    break
    else:
        try:
            DSL = pop.lociByNames(DSLtrait)
            numDSL = len(DSL)
        except:
            print 'Can not find one of the DSL %s in this population' % \
                ', '.join(DSLtrait)
            sys.exit(1)
        sign = []
        for i in DSL:
            af = 1 - pop.dvars().alleleFreq[i][0]
            # allele 1 is the minor allele
            if af < 0.5:
                sign.append(1)
            # allele 0 is the minor allele
            else:
                sign.append(-1)
    maf = [min(pop.dvars().alleleFreq[x][0], 1-pop.dvars().alleleFreq[x][0]) for x in DSL]
    print 'Using DSL %s with minor allele frequency %s' % (DSL, maf)
    #
    dsl = open(os.path.join(name, 'DSL.lst'), 'w')
    print >> dsl, "name\tchrom\tindex\tlocation\tminor allele frequency\tpercentage of variation"
    for i,d in enumerate(DSL):
        print >> dsl, "%s\t%d\t%d\t%.5f\t%.3f\t%.3f" % (pop.locusName(d),
            pop.chromLocusPair(d)[0] + 1, d - pop.chromBegin(pop.chromLocusPair(d)[0]), pop.locusPos(d),
            pop.dvars().alleleFreq[d][0], vars[i])
    dsl.close()
    # applying quantitative trait
    table = [[-math.sqrt(vars[i]/(2.*maf[i]*(1-maf[i])))*sign[i], 
        0, sign[i]*math.sqrt(vars[i]/(2.*maf[i]*(1-maf[i])))] \
            for i in range(numDSL)]
    def traitFunc(geno):
        return sum([table[x][geno[2*x]+geno[2*x+1]] for x in range(numDSL)]) + \
            random.normalvariate(0, math.sqrt(1-sum(vars)))
    pop.addInfoField('qtrait')
    print "Assigning quantitative trait..."
    PyQuanTrait(pop, loci=DSL, func=traitFunc)
    print "Mean quantitative trait is ", sum(pop.indInfo('qtrait', False))/pop.popSize()
    print "Minimal quantitative trait is ", min(pop.indInfo('qtrait', False))
    print "Maximal quantitative trait is ", max(pop.indInfo('qtrait', False))
    # set affection status
    qidx = pop.infoIdx('qtrait')
    for ind in pop.individuals():
        if ind.info(qidx) > cutoff:
            ind.setAffected(True)
        else:
            ind.setAffected(False)
    print 'Using cutoff value %.2f' % cutoff
    Stat(pop, numOfAffected=True)
    print "There are %d (%.2f percent) affected individuals." % (pop.dvars().numOfAffected, pop.dvars().numOfAffected*100.0/pop.popSize())


# penetrance generator functions. They will return a penetrance function
# with given penetrance parameter
def recessive(para):
    ''' recessive single-locus, heterogeneity multi-locus '''
    def func(geno):
        val = 1
        for i in range(len(geno)/2):
            if geno[i*2] + geno[i*2+1] == 2:
                val *= 1 - para[i]
        return 1-val
    return func
    

def additive(para):
    ''' additive single-locus, heterogeneity multi-locus '''
    def func(geno):
        val = 1
        for i in range(len(geno)/2):
            val *= 1 - (geno[i*2]+geno[i*2+1])*para[i]/2.
        return 1-val
    return func


# if you need some specialized penetrance function, modify this
# function here.
# NOTE:
# 
# 1. geno is the genptype at DSL. For example, if your DSL is [5,10]
#     geno will be something like [0,1,1,1] where 0,1 is the genotype at 
#     locus 5 and 1,1 is the genotype at locus 10.
# 2. in simuComplexDisease.py, 0 is wild type, 1 is disease allele.
def custom(para):
    ''' quantitative trait '''
    def func(geno):
        return 1
    return func


short_desc = '''This program simulates an admixed population based on 
two or more HapMap populations. Please follow the intructions
of the help message to prepare HapMap population.'''

# determine which script to run.
if __name__ == '__main__':
    # 
    # seed population does not exist, generate it
    allParam = getParam(options, short_desc, __doc__, nCol=2)
    # when user click cancel ...
    if len(allParam) == 0:
        sys.exit(1)
    # -h or --help
    if allParam[0]:    
        print usage(seed_options, __doc__)
        sys.exit(0)
    # 
    (name, 
      seed, useSavedSeed, initPop, HapMap_dir, pops, chrom,
        markerList, numMarkers, startPos, endingPos, minAF,
        minDiffAF, minDist, initCopy, gen, size,
      useSavedExpanded, expandGen, expandSize,
      useSavedAdmixed, migrModel, migrGen, migrRate,
      DSLtrait, chromWithDSL, freqDSL, freqDev, dslVar, cutoff,
      DSLpene, peneFunc, peneParam, 
      ccSampleSize, ccSampleName,
      randomSampleSize, randomSampleName) = allParam[1:]
    # simulation name?
    if not os.path.isdir(name):
        print 'Creating directory', name
        os.makedirs(name)
    if not os.path.isdir(name):
        print 'Can not create directory %s, exiting' % name
        sys.exit(1)
    cfgFile = os.path.join(name, name + '.cfg')
    print 'Save configuration to', cfgFile
    # save current configuration
    saveConfig(options, cfgFile, allParam)
    if os.path.isabs(seed):
        seedFile = seed
    else:
        seedFile = os.path.join(name, seed)
    expandedFile = os.path.join(name, 'expanded.bin')
    admixedFile = os.path.join(name, 'admixed.bin')
    # specified seed file?
    if type(useSavedSeed) == type(''):
        seedFile = useSavedSeed
        useSavedSeed = True
    if type(useSavedExpanded) == type(''):
        expandedFile = useSavedExpanded
        useSavedExpanded = True
    if type(useSavedAdmixed) == type(''):
        expandedFile = useSavedAdmixed
        useSavedAdmixed = True
    #
    # get seed population
    if useSavedSeed and os.path.isfile(seedFile):
        seedPop = None
    else:
        seedPop = generateSeedPopulation(HapMap_dir, initPop, 
            chrom, markerList, numMarkers, startPos,
            endingPos, minAF, minDiffAF, minDist, pops, initCopy,
            gen, size, seedFile)
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
        expandedPop = expandSeedPopulation(seedPop, expandGen,
            expandSize, expandedFile)
    # admixture
    if useSavedAdmixed and os.path.isfile(admixedFile):
        admixedPop = None
    else:
        if expandedPop is None:
            print 'Loading expanded population from file ', expandedFile
            expandedPop = LoadPopulation(expandedFile)
        admixedPop = mixExpandedPopulation(expandedPop, migrModel, migrGen,
            migrRate, admixedFile)
    #
    # step 3:
    # 
    if admixedPop is None:
        print 'Loading admixed population file ', admixedFile
        admixedPop = LoadPopulation(admixedFile)
    #
    print 'ch with dsl', chromWithDSL
    if len(chromWithDSL) > 0:
        # assign case/control status and quantitative trait
        setQuanTrait(admixedPop, DSLtrait, chromWithDSL, freqDSL, freqDev, 
            dslVar, cutoff, name)
    #
    if len(DSLpene) > 0:
        if 'recessive' == peneFunc:
            pene_func = recessive(para)
        elif 'additive' == peneFunc:
            pene_func = additive(para)
        elif 'custom' == peneFunc:
            pene_func = custom(para)
        else:
            pene_func = None
        #
        if pene_func is not None:
            PyPenetrance(admixedPop, admixedPop.lociByNames(DSLpene), func=pene_func)
    #
    # draw sample
    def comb(geno):
        return geno[0]+geno[1]
    if ccSampleSize != [0, 0]:
        print 'Drawing a case control sample with size ', ccSampleSize
        (samples,) = CaseControlSample(admixedPop, cases=ccSampleSize[0],
            controls=ccSampleSize[1])
        print "Saving case control sample to %s " % os.path.join(name, 'case_control')
        SaveQTDT(samples, output=os.path.join(name, 'case_control'), affectionCode=['1', '2'], 
                fields=['affection', 'qtrait'], combine=comb, header=True)
    if randomSampleSize != 0:
        print 'Drawing a random sample with size ', randomSampleSize
        (ran,) = RandomSample(admixedPop, size=randomSampleSize)
        # random sample
        print "Saving random sample to %s ..." % os.path.join(name, 'random')
        SaveQTDT(ran, output=os.path.join(name, 'random'), affectionCode=['1', '2'], 
            fields=['qtrait'], combine=comb, header=True) 


