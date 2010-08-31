#!/usr/bin/env python
#
# Author: Bo Peng
#
# Purpose: Simulation of samples with rare variants.
#
'''
Simulating a population of sequences forward in time, subject to mutation,
natural selection and population expansion. Because most mutants are introduced
during the rapid population expansion, most of the alleles will be rare at the
end of the simulation. Samples simulated using this script can be used to study
genetic diseases caused by a large number of rare variants.

'''

#
# This script simulates mutants in mutational space. That is to say, mutants
# are stored as locations on chromosomes, instead of allele number. Because
# indiviudals have different number of mutants, zeros are used to fill the
# rest of the markers.
#

import simuOpt
# 1.0.3 is required for ALL_AVAIL in PySelector to work
simuOpt.setOptions(alleleType='long', optimized=True, quiet=True, version='1.0.3')

import simuPOP as sim
from simuPOP.utils import ProgressBar
from simuPOP.sandbox import InfSitesSelector, InfSitesMutator

import types, os, sys, exceptions, logging, math, time

options = [
    {'separator': 'Genetic structure'},
    {'longarg': 'initPop=',
     'default': '',
     'label': 'Initial population',
     'description': '''Name of an initial population. If this file exists, it
        will be loaded and the evolution will start from this population, instead
        of a blank population. If this file does not exist, the initial
        population after the burnin stage will be saved to this file. This 
        option is often used to reduce long burn-in generations for replicate
        simulations.
     '''
    },
    {'longarg': 'regions=',
     'default': 'chr1:1..50000',
     'label': 'Regions',
     'description': '''A region (in basepair) means a piece of chromosome in
        which mutations can happen. A region should be expressed as chrXX:YYYY..ZZZZ
        where XX is chromosome number, YYYY is the starting position in basepair and
        ZZZZ is the ending position in basepair. The starting position should be at
        least one. If multiple regions are specified as a list of regions, they
        are assumed to be unlinked and will segregate independently even if they are
        on the same chromosome. Recombination within each region is not considered in
        this script because we assume that each region represents a narrow region of a
        chromosome that contains a gene.''',
     'allowedTypes': (types.ListType, types.TupleType),
    },
    {'longarg': 'addMutantsFrom=',
     'default': '',
     'label': 'Add mutants from population',
     'description': '''If a population is given, mutants from this population
        will be added to the population just before population expansion. Only
        loci that are within the specified regions will be selected. This
        population will be resized to bottleneck population size N1 before it
        is merged to the simulated population. This population is usually
        prepared using selectMarkers.py, using HapMap populations loaded using
        scripts loadHapMap2.py and loadHapMap3.py. These scripts are available
        from the simuPOP cookbook.''',
     'allowedTypes': types.StringType,
    },
    {'longarg': 'addMutantsGen=',
     'default': 0,
     'label': 'Generation at which to add mutants',
     'description': '''Generation number at which mutants from an external
        population will be inserted to the evolving population.''',
     'allowedTypes': [types.IntType, types.LongType],
    },
    {'separator': 'Initial population and demographic model'},
    {'longarg': 'N=',
     'default': [8100, 7900, 900000],
     'label': 'Population Sizes',
     'description': '''The size of population at the end of each stage. We assume
        that the first stage is a burn-in stage with constant population size. For
        other stages, a bottleneck model (N1->N2 in one generation followed by a
        constant population size of N2) will be used for population shrinkage and
        an exponential population exponential model will be used for population
        expansion.''',
     'allowedTypes': [types.ListType, types.TupleType],
    },
    {'longarg': 'numGen=',
     'default': [5000, 10, 370],
     'label': 'Numbers of generations at each stage',
     'description': '''Numbers of generations of different stages.''',
     'allowedTypes': [types.ListType, types.TupleType],
    },
    {'longarg': 'numStep=',
     'default': [],
     'label': 'Frequency of outputing simulation results',
     'description': '''Calculate and output statistics at intervals of specified
        number of generations. A single number or a list of numbers for each stage
        can be specified. If left unspecified, statistics at the beginning of
        each stage will be printed.''',
     'allowedTypes': [types.ListType, types.TupleType],
    },
    {'separator': 'Genetic forces'},
    {'longarg': 'mutationModel=',
     'default': 'infinite_sites',
     'label': 'Mutation model',
     'allowedTypes': types.StringType,
     'description': '''Mutation model. The default mutation model is an infinite-site
        model in which mutations that happen at previously mutated loci, even if those
        mutants have got lost, will be ignored. This model makes sure that all mutants
        in the simulated population can be traced to a single mutation event. In a
        finite sites model, mutations at all loci are allowed. If a mutant is mutated,
        it will be mutated from allele 1 to allele 0. Because both recurrent and
        reverse mutations are allowed, the histories of mutants in the simulated
        populations are more difficult to trace. Note that if the region is reasonably
        wide and the mutation rate is reasonable small, very few mutants will hit the
        same locus so these two models will yield very similar results.''',
     'chooseOneOf': ['infinite_sites', 'finite_sites'],
     'validate': simuOpt.valueOneOf(['infinite_sites', 'finite_sites']),
    },
    {'longarg': 'mu=',
     'default': 1.8e-8,
     'label': 'Mutation rate',
     'description': '''Mutation rate''',
     'allowedTypes': (types.IntType, types.LongType, types.FloatType),
     'validate': simuOpt.valueBetween(0., 1.),
    },
    {'longarg': 'selModel=',
     'default': 'multiplicative',
     'label': 'Multi-locus selection model',
     'chooseOneOf': ('multiplicative', 'additive', 'exponential'),
     'validate': simuOpt.valueOneOf(('multiplicative', 'additive', 'exponential')),
     'description': '''Multi-locus selection model, namely how to obtain an
        overall individual fitness after obtaining selection coefficients
        for all mutants (s_i, s_i < 0 when deleterious and s_i > 0 when adaptive).
        This script supports three models:
        |  multiplicative: prod (1 + s_i) Product of one minus selection coefficient.
        |  additive: max(0, 1 + sum(s_i)) One minus the sum of all selection 
            coefficients.
        |  exponential: exp(sum(s_i)) Exponential of minus sum of all selection
            coefficients.
        |Note that these models are equivalent if s_i are small and number of
        mutants are small.'''
    },
    {'longarg': 'selDist=',
     'default': 'constant',
     'label': 'Distribution of selection coefficient',
     'allowedTypes': types.StringType,
     'chooseOneOf': ['constant'] + ['gamma%d' % x for x in range(1,5)] + ['mixed_gamma'],
     'validate': simuOpt.valueOneOf(['constant'] + ['gamma%d' % x for x in range(1,5)] + ['mixed_gamma']),
     'description': '''Distribution of selection coefficient for new mutants.
        This script currently support the following distributions:
        |* constant: A single selection coefficient that gives each mutant a
            constant value s. The default parameter for this model is 
            -0.005. You can set selCoef to 0 to simulate neutral cases.
        |* gamma1: A basic gamma distribution (for -s) assuming a constant
            population size model (Eyre-Walker et al, 2006). The default
            parameters for this model is Pr(s=-x)=Gamma(0.23, 0.185).
        |* gamma2: A gamma distribution (for -s) assuming a two-epoch population
            size change model (Eyre-Walker et al, 2006). The default
            parameters for this model is Pr(s=-x)=Gamma(0.28, 0.05).
        |* gamma3: A gamma distribution (for -s) assuming a two-epoch population
            size change model for African population (Boyko et al, 2008). The
            default parameters for this model is Pr(s=-x)=Gamma(0.184, 0.160).
        |* gamma4: A gamma distribution (for -s) assuming a complex bottleneck
            model for European population (Boyko et al, 2008). The default
            parameters for this model is Pr(s=-x)=Gamma(0.206, 0.146).
        |* mixed_gamma: Parameter of this model should be a list of
            (a, p, k, theta) where a is the probability of having s=p (neutral
            or adptive sites), k, theta are the parameter of a gamma distribution
            Recomended parameter is (0.0186, 0.0001, 0.184, 0.160) for
            P(s=0.0001)=0.0186 and P(s=-x)=(1-0.0186)*Gamma(0.184,0.160).
        |If you would like to define your own selection model, please have a
            look at class InfiniteSitesSelector. You can either modify function
            getGamma or define your own gammaFitness, getGamma and set
            self.func properly.
        |NOTE:
        |The mixed_gamma model is implemented as a user-defined Python function
        that return a random selection coefficient. You can define your own function
        to define arbitrary distribution of selection coefficients and pass this function
        to parameter selDist of function simuRareVariants() if you import this script
        and call function simuRareVariants directly. Furthermore, if this
        function has a parameter loc, the positions of the mutant will be passed
        to this function so that you can return site-specific coefficients.'''
    },
    {'longarg': 'selCoef=',
     'default': None,
     'label': 'Customized selection coefficient',
     'allowedTypes': (types.NoneType, types.IntType,
        types.LongType, types.FloatType, types.ListType, types.TupleType),
     'validate': simuOpt.valueOr(simuOpt.valueEqual(None), simuOpt.valueOr(
        simuOpt.valueIsNum(), simuOpt.valueListOf(simuOpt.valueIsNum()))),
     'description': '''Selection coefficient with its meaning determined by
        parameter selDist. If None is given, the default parameter for the
        selected distribution will be used.''',
    },
    {'separator': 'Output'},
    {'longarg': 'popFile=',
     'default': 'output.pop',
     'label': 'Save simulated population as',
     'allowedTypes': types.StringType,
     'description': '''Filename to which the simulated population will be saved
        in simuPOP format. The file can be loaded later for post processing
        purposes.'''
    },
    {'longarg': 'markerFile=',
     'default': 'output.map',
     'label': 'Save marker information to',
     'allowedTypes': types.StringType,
     'description': '''Filename to which the marker information, including
        marker name (reg+index), chromosome, location, allele frequency and
        selection coefficient are saved. Monomorphic markers are ignored.'''
    },
    {'longarg': 'mutantFile=',
     'default': 'output.mut',
     'label': 'Save mutants to',
     'allowedTypes': types.StringType,
     'description': '''Filename to which the mutants are outputed. The file
        will be saved in the format of 
        |  ind_id mut1 mut2 ...
        |where ind_id is the index of individual (1, 2, ...), mut1 and mut2 are
        locations of mutants. Haplotypes for different regions and homologous
        chromosomes are saved in different lines in the order of 
        | reg1_ploidy1
        | reg1_ploidy2
        | reg2_ploidy1
        | reg2_ploidy2
        | ...
        ''',
    },
    {'longarg': 'genotypeFile=',
     'default': '',
     'label': 'Save genotype (in .ped format) to',
     'allowedTypes': types.StringType,
     'description': '''Filename to which the genotypes of all individuals are
        saved. The file will be saved in the format of
        |    famid id fa mo sex aff loc1_a1 loc1_a2 loc2_a1 loc2_a2 ...
        |where famid is 1, 2, 3, ... id is always 1, fa, mo is always 0.
        Wildtype and mutant alleles are denoted by 0 and 1 respectively. This
        option is turned off by default because this format is not efficient in
        storing a small number of mutants.'''
    },
    {'longarg': 'verbose=',
     'default': 1,
     'allowedTypes': types.IntType,
     'validate': simuOpt.valueBetween(0, 2),
     'description': '''0 for quiet, 1 for regular output, 2 for debug output.
        In the debug output, a file 'mutations.lst' will be saved with all 
        mutation events. This option is not visible from gui.''',
    },
]


class NumSegregationSites(sim.PyOperator):
    '''A Python operator to count the number of segregation sites (number of
    distinct mutants), average number of segreagation sites of individuals,
    and average allele frequency of these mutants. The results are saved
    in variables ``numSites``, ``avgSites`` and ``avgFreq``.
    '''
    def __init__(self, *args, **kwargs):
        sim.PyOperator.__init__(self, func=self.countSites, *args, **kwargs)

    def countSites(self, pop):
        '''Count the number of segregation sites, average sites per individual,
        average allele frequency.'''
        nFixed = self.countFixedMutants(pop)
        geno = pop.genotype()
        numMutants = float(len(geno) - geno.count(0) - nFixed * pop.popSize() * 2) 
        numSites = len(set(geno)) - nFixed - 1
        if numMutants == 0:
            avgFreq = 0
        else:
            avgFreq = numMutants / numSites / (2*pop.popSize())
        pop.dvars().numSites = numSites
        pop.dvars().avgSites = float(numMutants) / pop.popSize()
        pop.dvars().avgFreq = avgFreq
        return True

    def countFixedMutants(self, pop):
        '''Find mutants that exist in all individuals.
        '''
        # initial set...
        commonMutants = set(pop.individual(0).genotype(0))
        commonMutants.discard(0)
        for ind in pop.individuals():
            commonMutants &= set(ind.genotype(0))
            commonMutants &= set(ind.genotype(1))
            if len(commonMutants) == 0:
                return 0
        #print '%d fixed mutants are found' % len(commonMutants)
        return len(commonMutants)


def mutantsToAlleles(pop, logger):
    '''Convert a population from mutational space to allele space. Monomorphic
    markers are ignored.
    '''
    # figure out chromosomes and markers
    markers = {}
    for ch,region in enumerate(pop.chromNames()):
        chNumber = region.split(':')[0][3:]
        loci = set()
        for ind in pop.individuals():
            loci |= set(ind.genotype(0, ch))
            loci |= set(ind.genotype(1, ch))
        if markers.has_key(chNumber):
            markers[chNumber] |= loci
        else:
            markers[chNumber] = loci
    # create a population for each chromosome
    pops = []
    chroms = markers.keys()
    chroms.sort()
    for ch in chroms:
        markers[ch] = list(markers[ch])
        markers[ch].remove(0)
        markers[ch].sort()
        if logger:
            logger.info('Chromosome %s has %d markers' % (ch, len(markers[ch])))
        apop = sim.Population(pop.popSize(), loci=len(markers[ch]),
            lociPos=markers[ch])
        # get a dictionary of loci position
        lociPos = {}
        for idx,loc in enumerate(apop.lociPos()):
            lociPos[loc] = idx
        for aind,mind in zip(apop.individuals(), pop.individuals()):
            for p in range(2):
                for mutant in mind.genotype(p):
                    if mutant != 0:
                        aind.setAllele(1, lociPos[mutant], p)
        pops.append(apop)
    for pop in pops[1:]:
        pops[0].addChromFrom(pop)
    return pops[0]


def allelesToMutants(pop, regions, logger=None):
    '''Convert a population from allele space to mutational space, using
    specified regions.
    '''
    pops = []
    for region in regions:
        loci = []
        ch_name = region.split(':')[0][3:]
        start, end = [int(x) for x in region.split(':')[1].split('..')]
        try:
            ch = pop.chromByName(ch_name)
        except:
            raise exceptions.ValueError('Chromosome %s is not available in passed population.' % ch_name)
        for loc in range(pop.chromBegin(ch), pop.chromEnd(ch)):
            pos = pop.locusPos(loc)
            if pos >= start and pos <= end:
                loci.append(loc)
        # get the mutants for each individual
        allAlleles = []
        for ind in pop.individuals():
            alleles0 = []
            alleles1 = []
            for loc in loci:
                if ind.allele(loc, 0) != 0:
                    alleles0.append(int(pop.locusPos(loc)))
                if ind.allele(loc, 1) != 0:
                    alleles1.append(int(pop.locusPos(loc)))
            allAlleles.extend([alleles0, alleles1])
        # maximum number of mutants
        maxMutants = max([len(x) for x in allAlleles])
        if logger is not None:
            logger.info('%d loci are identified with at most %d mutants in region %s.' % (len(loci), maxMutants, region))
        # create a population
        mpop = sim.Population(pop.popSize(), loci=maxMutants, chromNames=region)
        # put in mutants
        for idx,ind in enumerate(mpop.individuals()):
            geno = ind.genotype(0)
            for loc,mut in enumerate(allAlleles[idx*2]):
                geno[loc] = mut
            geno = ind.genotype(1)
            for loc,mut in enumerate(allAlleles[idx*2+1]):
                geno[loc] = mut
        pops.append(mpop)
    # merge all populations into one
    for pop in pops[1:]:
        pops[0].addChromFrom(pop)
    return pops[0]

def addMutants(pop, param):
    # Adding mutants
    addMutantsFrom, regions, logger = param
    mPop = sim.loadPopulation(addMutantsFrom)
    # convert allele-based population to mutation based population.
    mPop = allelesToMutants(mPop, regions, logger)
    #
    mPop.resize(pop.popSize())
    # Add loci to pop
    for ch in range(mPop.numChrom()):
        pop.addLoci([ch]*mPop.numLoci(ch), range(pop.numLoci(ch) + 1, pop.numLoci(ch) + mPop.numLoci(ch) + 1))
    if logger:
        # if an initial population is given
        logger.info('Adding mutants to population after bottleneck')
    # Add mutants to pop
    for ind, mInd in zip(pop.individuals(), mPop.individuals()):
        for p in range(2):
            for ch in range(pop.numChrom()):
                geno = ind.genotype(p, ch)
                mGeno = mInd.genotype(p, ch)
                idx = geno.index(0)
                for i,m in enumerate(mGeno):
                    if m == 0:
                        break
                    geno[idx + i] = m
    return True


###
### The container.Counter class only exist in Python 2.7 so I put it here
###
from operator import itemgetter
from heapq import nlargest
from itertools import repeat, ifilter

class Counter(dict):
    '''Dict subclass for counting hashable objects.  Sometimes called a bag
    or multiset.  Elements are stored as dictionary keys and their counts
    are stored as dictionary values.

    >>> Counter('zyzygy')
    Counter({'y': 3, 'z': 2, 'g': 1})

    '''

    def __init__(self, iterable=None, **kwds):
        '''Create a new, empty Counter object.  And if given, count elements
        from an input iterable.  Or, initialize the count from another mapping
        of elements to their counts.

        >>> c = Counter()                           # a new, empty counter
        >>> c = Counter('gallahad')                 # a new counter from an iterable
        >>> c = Counter({'a': 4, 'b': 2})           # a new counter from a mapping
        >>> c = Counter(a=4, b=2)                   # a new counter from keyword args

        '''        
        self.update(iterable, **kwds)

    def __missing__(self, key):
        return 0


    def elements(self):
        '''Iterator over elements repeating each as many times as its count.

        >>> c = Counter('ABCABC')
        >>> sorted(c.elements())
        ['A', 'A', 'B', 'B', 'C', 'C']

        If an element's count has been set to zero or is a negative number,
        elements() will ignore it.

        '''
        for elem, count in self.iteritems():
            for _ in repeat(None, count):
                yield elem

    # Override dict methods where the meaning changes for Counter objects.

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError(
            'Counter.fromkeys() is undefined.  Use Counter(iterable) instead.')

    def update(self, iterable=None, **kwds):
        '''Like dict.update() but add counts instead of replacing them.

        Source can be an iterable, a dictionary, or another Counter instance.

        >>> c = Counter('which')
        >>> c.update('witch')           # add elements from another iterable
        >>> d = Counter('watch')
        >>> c.update(d)                 # add elements from another counter
        >>> c['h']                      # four 'h' in which, witch, and watch
        4

        '''        
        if iterable is not None:
            if hasattr(iterable, 'iteritems'):
                if self:
                    self_get = self.get
                    for elem, count in iterable.iteritems():
                        self[elem] = self_get(elem, 0) + count
                else:
                    dict.update(self, iterable) # fast path when counter is empty
            else:
                self_get = self.get
                for elem in iterable:
                    self[elem] = self_get(elem, 0) + 1
        if kwds:
            self.update(kwds)

    def __delitem__(self, elem):
        'Like dict.__delitem__() but does not raise KeyError for missing values.'
        if elem in self:
            dict.__delitem__(self, elem)

    def __repr__(self):
        if not self:
            return '%s()' % self.__class__.__name__
        items = ', '.join(map('%r: %r'.__mod__, self.most_common()))
        return '%s({%s})' % (self.__class__.__name__, items)
 
#
# End of copied code
#


def saveMarkerInfoToFile(pop, filename, logger=None):
    '''Save a map file with an additional column of allele frequency. The
    population has to be in mutational space. This function assumes that
    there is a variable selCoef in this population which contains selection
    coefficients for all mutants.
    '''
    allCounts = [Counter() for x in range(pop.numChrom())]
    prog = ProgressBar('Counting number of mutants', pop.popSize())
    for ind in pop.individuals():
        # there can be memory problem....
        for ch in range(pop.numChrom()):
            allCounts[ch].update(ind.genotype(0, ch))
            allCounts[ch].update(ind.genotype(1, ch))
        prog.update()
    allMutants = []
    selCoef = pop.dvars().selCoef
    if filename:
        map = open(filename, 'w')
        print >> map, 'name\tchrom\tposition\tfrequency\tselCoef'
    for ch,region in enumerate(pop.chromNames()):
        # real chromosome number
        chName = region.split(':')[0][3:]
        counts = allCounts[ch]
        # get markers
        mutants = counts.keys()
        mutants.sort()
        # allele 0 is fake
        if mutants[0] == 0:
            mutants = mutants[1:]
        allMutants.append(mutants)
        if filename:
            # find all markers
            sz = pop.popSize() * 2.
            for idx,marker in enumerate(mutants):
                print >> map, 'loc%d_%d\t%s\t%d\t%.8f\t%.8f' % (ch + 1, idx + 1, chName, marker,
                    counts[marker] / sz, selCoef[marker] if type(selCoef) == type({}) else selCoef)
    if filename:
        map.close()
    return allMutants
        

def saveMutantsToFile(pop, filename, infoFields=[], logger=None):
    '''Save haplotypes as a list of mutant locations to file, in the format of
       ind_idx reg_id FIELDS mut1 mut2 ...
    where FIELDS are information fields.
    '''
    mut = open(filename, 'w')
    prog = ProgressBar('Writing mutants of %d individuals to %s' % (pop.popSize(), filename), pop.popSize())
    for idx,ind in enumerate(pop.allIndividuals()):
        fields = ' '.join([str(ind.info(x)) for x in infoFields])
        for ch in range(pop.numChrom()):
            geno = list(ind.genotype(0, ch))
            geno.sort()
            print >> mut, idx+1, fields, ' '.join([str(x) for x in geno if x != 0])
            geno = list(ind.genotype(1, ch))
            geno.sort()
            if geno[0] == 0:
                geno = geno[1:]
            print >> mut, idx+1, fields, ' '.join([str(x) for x in geno if x != 0])
        prog.update()
    mut.close()

def saveGenotypeToFile(pop, filename, allMutants, logger=None):
    '''Save genotype in .ped file format. Because there is no family structure, we have
        famid = 1, 2, 3, ...
        id = 1
        fa = 0
        ma = 0
        sex = 1 for male and 2 for female
        aff = 1 for unaffected and 2 for affected
        genotype 

    allMutants:
        lists of mutants returned by function markerFile
    '''
    if logger:
        logger.info('Saving genotype to %s in standard .ped format.' % filename)
    ped = open(filename, 'w')
    # marker index...
    markerPos = []
    for mutants in allMutants:
        pos = {}
        for idx,m in enumerate(mutants):
            pos[m] = idx
        markerPos.append(pos)
    prog = ProgressBar('Writing genotype of %d individuals to %s' % (pop.popSize(), filename), pop.popSize())
    sexCode = {sim.MALE: 1, sim.FEMALE: 2}
    affCode = {False: 1, True: 2}
    for cnt, ind in enumerate(pop.individuals()):
        print >> ped, '%s 0 0 %d %d' % (cnt + 1, sexCode[ind.sex()], affCode[ind.affected()]),
        for ch in range(pop.numChrom()):
            # a blank genotype
            geno = [0]*(len(markerPos[ch])*2)
            # add 1 according to mutant location (first ploidy)
            for m in ind.genotype(0, ch):
                if m == 0:
                    break
                geno[2*markerPos[ch][m]] = 1
            # add 1 according to mutant location (second ploidy)
            for m in ind.genotype(1, ch):
                if m == 0:
                    break
                geno[2*markerPos[ch][m]+1] = 1
            print >> ped, ' '.join([str(x) for x in geno]),
        print >> ped
        prog.update()
    ped.close()

class fitnessCollector:
    '''This is a simple connection class that gets output from 
    a InfSiteSelector and collect mutant fitness'''
    def __init__(self):
        self.selCoef = {}

    def getCoef(self, lines):
        for line in lines.strip().split('\n'):
            mut, sel = line.split()
            self.selCoef[int(mut)] = float(sel)


def mixedGamma():
    '''This function returns a random fitness value for a new mutant
    according to a mixed_gamma distribution. If a parameter loc is defined,
    locus index will be passed so that you can return different selection
    coefficient for different locations.
    '''
    if sim.getRNG().randUniform() < selCoef[0]:
        return selCoef[1]
    while True:
        s = sim.getRNG().randGamma(selCoef[2], selCoef[3])
        if s > 0.0001 and s < 0.01:
            return -s


def simuRareVariants(regions, N, numGen, mu, selModel, selDist, selCoef,
        numStep=[], mutationModel='finite_sites', initPop='', addMutantsFrom='', addMutantsGen=0,
        popFile='', markerFile='', mutantFile='', genotypeFile='', verbose=1, logger=None):
    '''
    Please refer to simuRareVariants.py -h for a detailed description of all parameters.
    Note that a user-defined function can be passed to parameter selDist to specify
    arbitrary distribution of fitness.
    '''
    #
    # convert regions to start/end positions
    ranges = []
    for region in regions:
        start, end = [int(x) for x in region.split(':')[1].split('..')]
        ranges.append((start, end+1))
    if logger:
        logger.info('%s regions with a total length of %d basepair.' % (len(ranges), sum([x[1]-x[0] for x in ranges])))
    #
    # set default parameter
    if selCoef is None:
        # set default parameters
        if selDist == 'mixed_gamma':
            selCoef = (0.0186, 0.0001, 0.184, 0.160)
        elif selDist == 'gamma1':
            selCoef = (0.23, 0.185)
        elif selDist == 'gamma2':
            selCoef = (0.28, 0.005)
        elif selDist == 'gamma3':
            selCoef = (0.184, 0.160)
        elif selDist == 'gamma4':
            selCoef = (0.206, 0.146)
        elif selDist == 'constant':
            selCoef = -0.005
        elif not callable(selDist):
            raise exceptions.ValueError("Unsupported random distribution")
    # check if lists obj 'N', 'numGen' and 'numStep' have the same length
    if len(numStep) == 0:
        numStep = numGen
    elif len(numStep) == 1:
        numStep = [numStep] * len(numGen)
    if len(N) != len(numGen) or len(numGen) !=len(numStep):
        raise ValueError("N:sizes of population, numGen:numbers of generations and numStep:numbers of steps "
            "of saving simulation results, should be lists that have the equal length")
    # use a right selection operator.
    collector = fitnessCollector()
    mode = {'multiplicative': sim.MULTIPLICATIVE,
        'additive': sim.ADDITIVE,
        'exponential': sim.EXPONENTIAL}[selModel]
    #
    if callable(selDist):
        mySelector = InfSitesSelector(selDist=selDist, mode=mode, output=collector.getCoef)
    elif selDist.startswith('gamma'):
        mySelector = InfSitesSelector(selDist=(sim.GAMMA_DISTRIBUTION, selCoef[0], selCoef[1]),
            mode=mode, output=collector.getCoef)
    elif selDist == 'constant':
        if selCoef == 0:
            mySelector = sim.NoneOp()
        else:
            mySelector = InfSitesSelector(selDist=(sim.CONSTANT, selCoef),
                mode=mode, output=collector.getCoef)
    elif selDist == 'mixed_gamma':
        mySelector = InfSitesSelector(selDist=mixedGamma, mode=mode, output=collector.getCoef)
    #
    # 0, G[0], G[0]+G[1], ..., sum(G)
    sumNumGen = [sum(numGen[:i]) for i in range(len(numGen)+1)]
    # the demographic model: N[0] = the population size of the burnin generation
    def demoFunc(gen):
        for i in range(len(numGen)):
            if sumNumGen[i] <= gen < sumNumGen[i+1]:
                # at the burnin stage
                if i == 0:
                    return N[i]
                # at any bottleneck stage
                elif N[i] <= N[i-1]:
                    return N[i]
                # at any expansion stage
                else:
                    # to make sure that the last generation of this expansion stage 
                    # has the exact required number of individuals.
                    if gen == sumNumGen[i+1] - 1:
                        return N[i]
                    else:
                        r = math.log(N[i] * 1.0 / N[i-1]) / numGen[i]
                        return int(N[i-1] * math.exp(r*(gen - sumNumGen[i])))
        # after the last generation, this should not happen, but...
        return N[-1]
    # Evolve
    if os.path.isfile(initPop):
        if logger:
            logger.info('Loading initial population %s...' % initPop)
        pop = sim.loadPopulation(initPop)
        if pop.numChrom() != len(regions):
            raise exceptions.ValueError('Initial population %s does not have specified regions.' % initPop)
        for ch,reg in enumerate(regions):
            if pop.chromName(ch) != reg:
                raise exceptions.ValueError('Initial population %s does not have region %s' % (initPop, reg))
    else:
        pop = sim.Population(size=N[0], loci=[10]*len(regions), chromNames=regions,
            infoFields='fitness')
    if logger:
        startTime = time.clock()

    atGen = []
    for i in range(len(sumNumGen)-1):
        atGen += range(sumNumGen[i], sumNumGen[i+1], numStep[i])
    pop.evolve(
        initOps=sim.InitSex(),
        preOps=[
            sim.PyOutput('Statistics outputted are Generation number, population size, '
                         'number of segregation sites, average number of sites per individual, '
                         'average allele frequency, average fitness value, and minimal fitness value '
                         'of the parental population\n', at = 0)] + \
            [sim.PyOutput('Entering stage %d\n' % i, at = sumNumGen[i]) for i in range(1, len(sumNumGen))] + \
            # mutate in a region at rate mu, if verbose > 2, save mutation events to a file
            [InfSitesMutator(mu, ranges, {'infinite_sites':1, 'finite_sites':2}[mutationModel], output='' if verbose < 2 else '>>mutations.lst'),
            # selection on all loci
            mySelector,
            # output statistics in verbose mode
            sim.IfElse(verbose > 0, ifOps=[
                sim.Stat(popSize=True, meanOfInfo='fitness', minOfInfo='fitness'),
                NumSegregationSites(),
                sim.PyEval(r'"%5d %8d %5d %.6f %.6f %.6f %.6f\n" '
                    '% (gen, popSize, numSites, avgSites, avgFreq*100, meanOfInfo["fitness"], minOfInfo["fitness"])'),
                ], at = atGen
            ),
        ],
        matingScheme=sim.RandomMating(subPopSize=demoFunc),
        postOps = [
            # save a population if an initial population is given
            sim.IfElse(isinstance(initPop, str) and initPop != '',
                ifOps = [
                    sim.PyOutput('Saving initial population (after burn-in) tp %s' % initPop),
                    sim.SavePopulation(output=initPop if initPop else 'dummy'),
                ], at = N[0],
            ),
            # add alleles from an existing population before population expansion
            sim.IfElse(addMutantsFrom != '' and addMutantsGen != 0,
                ifOps = [
                    sim.PyOutput('Loading and converting population %s' % addMutantsFrom),
                    sim.PyOperator(func=addMutants, param=(addMutantsFrom, regions, logger)),
                ], at = addMutantsGen),
        ],
        finalOps=[
            # statistics after evolution
            sim.Stat(popSize=True),
            NumSegregationSites(),
            sim.PyEval(r'"Simulated population has %d individuals, %d segregation sites.'
                       r'There are on average %.1f sites per individual. Mean allele frequency is %.4f%%.\n"'
                       r'% (popSize, numSites, avgSites, avgFreq*100)'),
        ],
        gen = sumNumGen[-1]
    )
    # record selection coefficients to population
    if len(collector.selCoef) == 0:
        # this must be the neutral case where a NonOp has been used.
        pop.dvars().selCoef = 0
    else:
        pop.dvars().selCoef = collector.selCoef
    #
    if logger:
        logger.info('Population simulation takes %.2f seconds' % (time.clock() - startTime))
    if popFile:
        if logger:
            logger.info('Saving population in %s (in simuPOP format)' % popFile)
        pop.save(popFile)
    if markerFile or genotypeFile:
        if logger:
            logger.info('Saving marker information to file %s' % markerFile)
        mutants = saveMarkerInfoToFile(pop, markerFile, logger)
        if genotypeFile:
            if logger:
                logger.info('Saving genotype in .ped format to file %s' % genotypeFile)
            saveGenotypeToFile(pop, genotypeFile, mutants, logger)
    if mutantFile:
        if logger:
            logger.info('Saving mutants to file %s' % mutantFile)
        saveMutantsToFile(pop, mutantFile, logger)
    return pop


if __name__ == '__main__':
    pars = simuOpt.Params(options, 'A simulator of rare variants.', __doc__)
    if not pars.getParam():
        sys.exit(1)
    # if only one region is given.
    if isinstance(pars.regions[0], (int, long)):
        pars.regions = (pars.regions, )
    # if you change level to logging.DEBUG, a lot of debug information, including
    # location and generation of mutants will be outputed. You can also use parameter
    # filename of function basicConfig to output log to a file.
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('simuRareVariants')
    simuRareVariants(logger=logger, **pars.asDict())


