#!/usr/bin/env python

#
# $File: __init__.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit http://simupop.sourceforge.net
# for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
"""
simuPOP core module.
"""

__all__ = [
    '__version__',
    # Constants
    #
    'MALE',
    'FEMALE',
    #
    'CUSTOMIZED',
    'AUTOSOME',
    'CHROMOSOME_X',
    'CHROMOSOME_Y',
    'MITOCHONDRIAL',
    #
    'CONSTANT',
    'BINOMIAL_DISTRIBUTION',
    'EXPONENTIAL_DISTRIBUTION',
    'GEOMETRIC_DISTRIBUTION',
    'POISSON_DISTRIBUTION',
    'UNIFORM_DISTRIBUTION',
    'NORMAL_DISTRIBUTION',
    'GAMMA_DISTRIBUTION',
    #
    'NO_SEX',
    'RANDOM_SEX',
    'PROB_OF_MALES',
    'NUM_OF_MALES',
    'NUM_OF_FEMALES',
    'SEQUENCE_OF_SEX',
    'GLOBAL_SEQUENCE_OF_SEX',
    #
    'NO_CONVERSION',
    'NUM_MARKERS',
    'TRACT_LENGTH',
    #
    'SPOUSE',
    'OUTBRED_SPOUSE',
    'FULLSIBLING',
    'SIBLING',
    'OFFSPRING',
    'COMMON_OFFSPRING',
    #
    'ANY_SEX',
    'MALE_ONLY',
    'FEMALE_ONLY',
    'SAME_SEX',
    'OPPOSITE_SEX',
    'PAIR_ONLY',
    #
    'MULTIPLICATIVE',
    'ADDITIVE',
    'HETEROGENEITY',
    'EXPONENTIAL',
    #
    'BY_IND_INFO',
    'BY_PROBABILITY',
    'BY_PROPORTION',
    'BY_COUNTS',
    #
    'PATERNAL',
    'MATERNAL',
    #
    'MEAN',
    'MAXIMUM',
    'MINIMUM',
    'SUMMATION',
    'MULTIPLICATION',
    #
    'PER_LOCI',
    'PER_CHROMOSOME',
    'PER_PLOIDY',
    'PER_INDIVIDUAL',
    'FROM_INFO',
    'FROM_INFO_SIGNED',
    #
    'HAPLODIPLOID',
    'ALL_AVAIL',
    'UNSPECIFIED',
    #
    # type
    'defdict',
    #
    # Major simuPOP classes
    'Population',
    # This is just to make help(Individual) available to users.
    'Individual',
    'Simulator',
    'Pedigree',
    # splitters
    'SexSplitter',
    'AffectionSplitter',
    'CombinedSplitter',
    'ProductSplitter',
    'ProportionSplitter',
    'InfoSplitter',
    'GenotypeSplitter',
    'RangeSplitter',
    # mating schemes
    'MatingScheme',
    'HomoMating',
    'HeteroMating',
    'ConditionalMating',
    'PedigreeMating',
    'OffspringGenerator',
    'ControlledOffspringGenerator',
    'RandomParentChooser',
    'PyParentsChooser',
    'CombinedParentsChooser',
    'SequentialParentChooser',
    'RandomParentsChooser',
    'SequentialParentsChooser',
    'ParentChooser',
    'PolyParentsChooser',
    #
    'RandomMating',
    'RandomSelection',
    'MonogamousMating',
    'SelfMating',
    'HermaphroditicMating',
    'CloneMating',
    'HaplodiploidMating',
    #'consanguineousMating',
    'ControlledRandomMating',
    'PolygamousMating',
    #
    #
    # Operators
    #
    'InitSex',
    'InitInfo',
    'InitGenotype',
    'InitLineage',
    #
    'PyOutput',
    'PyEval',
    'PyExec',
    'InfoEval',
    'InfoExec',
    #
    'Migrator',
    'BackwardMigrator',
    'MergeSubPops',
    'SplitSubPops',
    'ResizeSubPops',
    #
    'GenoTransmitter',
    'HaplodiploidGenoTransmitter',
    'MendelianGenoTransmitter',
    'MitochondrialGenoTransmitter',
    'SelfingGenoTransmitter',
    'CloneGenoTransmitter',
    'Recombinator',
    #
    'PointMutator',
    'MatrixMutator',
    'MixedMutator',
    'PyMutator',
    'StepwiseMutator',
    'SNPMutator',
    'AcgtMutator',
    'AminoAcidMutator',
    'ContextMutator',
    'KAlleleMutator',
    'RevertFixedSites',
    'FiniteSitesMutator',
    #
    'MapSelector',
    'MaSelector',
    'MlSelector',
    'PySelector',
    'PyMlSelector',
    #
    'MaPenetrance',
    'MapPenetrance',
    'MlPenetrance',
    'PyPenetrance',
    'PyMlPenetrance',
    #
    'PyQuanTrait',
    #
    'Stat',
    #
    'IdTagger',
    'InheritTagger',
    'ParentsTagger',
    'OffspringTagger',
    'PedigreeTagger',
    'PyTagger',
    'SummaryTagger',
    #
    'TerminateIf',
    'RevertIf',
    'DiscardIf',
    'discardIf',
    #
    'PyOperator',
    #
    'NoneOp',
    'Dumper',
    'SavePopulation',
    'IfElse',
    'Pause',
    'TicToc',
    #
    # Function form of operators
    'dump',
    #
    'initSex',
    'initInfo',
    'initGenotype',
    'initLineage',
    #
    'pyEval',
    'pyExec',
    'infoEval',
    'infoExec',
    #
    'tagID',
    #
    'migrate',
    'backwardMigrate',
    'resizeSubPops',
    'splitSubPops',
    'mergeSubPops',
    #
    'matrixMutate',
    'contextMutate',
    'kAlleleMutate',
    'mixedMutate',
    'pointMutate',
    'pyMutate',
    'stepwiseMutate',
    'snpMutate',
    'acgtMutate',
    'revertFixedSites',
    #
    'stat',
    #
    # Global functions
    'WithArgs',
    'WithMode',
    # RNG Related
    'RNG',
    'getRNG',
    'setRNG',
    #
    'closeOutput',
    'describeEvolProcess',
    'loadPopulation',
    'loadPedigree',
    'moduleInfo',
    'turnOffDebug',
    'turnOnDebug',
    'setOptions',
    #
    'maPenetrance',
    'mapPenetrance',
    'mlPenetrance',
    'pyPenetrance',
    #
    'mapSelect',
    'maSelect',
    'mlSelect',
    'pySelect',
    'pyMlSelect',
    #
    'pyQuanTrait',
    #
    # For testing only
    'applyDuringMatingOperator',
    'Bernullitrials',
    'Bernullitrials_T',
    'WeightedSampler',
    #
    # modules are not loaded by default because they require rpy or matplotlib
    #
    #'utils',
    #'sampling',
]

# get options
from simuOpt import simuOptions
import os, sys, re

if simuOptions['Optimized']:
    if simuOptions['AlleleType'] == 'short':
        from simuPOP.simuPOP_op import *
    elif simuOptions['AlleleType'] == 'long':
        from simuPOP.simuPOP_laop import *
    elif simuOptions['AlleleType'] == 'binary':
        from simuPOP.simuPOP_baop import *
    elif simuOptions['AlleleType'] == 'mutant':
        from simuPOP.simuPOP_muop import *
    elif simuOptions['AlleleType'] == 'lineage':
        from simuPOP.simuPOP_linop import *
    else:
        from simuPOP.simuPOP_op import *
else:
    if simuOptions['AlleleType'] == 'short':
        from simuPOP.simuPOP_std import *
    elif simuOptions['AlleleType'] == 'long':
        from simuPOP.simuPOP_la import *
    elif simuOptions['AlleleType'] == 'binary':
        from simuPOP.simuPOP_ba import *
    elif simuOptions['AlleleType'] == 'mutant':
        from simuPOP.simuPOP_mu import *
    elif simuOptions['AlleleType'] == 'lineage':
        from simuPOP.simuPOP_lin import *
    else:
        from simuPOP.simuPOP_std import *

__version__ = moduleInfo()['version']

    
if simuOptions['Version'] is not None:
    expMajor, expMinor, expRelease = [
        int(x) for x in re.match(r'^(\d+)\.(\d+)\.(\d+)',
                                 simuOptions['Version']).groups()
    ]
    myMajor, myMinor, myRelease = [
        int(x) for x in re.match(r'^(\d+)\.(\d+)\.(\d+)', __version__).groups()
    ]
    if (expMajor > myMajor) or (expMajor == myMajor and expMinor > myMinor) or \
        (expMajor == myMajor and expMinor == myMinor and expRelease > myRelease):
        raise ImportError('simuPOP version %s is installed but version >= %s is required. ' % \
            (__version__, simuOptions['Version']) +
            'Please upgrade your simuPOP installation.')

if simuOptions['Revision'] is not None:
    rev = moduleInfo()['revision']
    if rev < simuOptions['Revision']:
        raise ImportError('simuPOP version %s (revision %d) is installed ' %
                          (__version__, rev) +
                          'but simuPOP revision >= %d is required. ' %
                          simuOptions['Revision'] +
                          'Please upgrade your simuPOP installation.')
    if 'DBG_COMPATIBILITY' in simuOptions['Debug']:
        print(
            'WARNING: parameter revision in simuOpt.setOptions is obsolete and might be removed from a future version of simuPOP.', file=sys.stderr
        )

# set number of threads in openMP
if simuOptions['NumThreads'] is not None:
    setOptions(numThreads=simuOptions['NumThreads'])

if not simuOptions['Quiet']:
    info = moduleInfo()
    print("simuPOP Version %s : Copyright (c) 2004-2016 Bo Peng" %
          (__version__), file=sys.stderr)
    # compile date, compiler etc are macros that are replaced during compile time.
    print("Revision %d (%s) for Python %s (%dbit, %d%s)" % \
            (info['revision'], info['date'], info['python'], info['wordsize'], info['threads'],
                    'thread' if info['threads'] <= 1 else 'threads'), file=sys.stderr)
    print("Random Number Generator is set to %s with random seed 0x%08x." %
          (getRNG().name(), getRNG().seed()), file=sys.stderr)
    # MaxAllele + 1 since 0 is one of the allelic states
    print("This is the %s %s allele version with %d maximum allelic states." % \
            ('optimized' if info['optimized'] else 'standard', info['alleleType'], info['maxAllele']+1), file=sys.stderr)
    print("For more information, please visit http://simupop.sourceforge.net,", file=sys.stderr)
    print(
        "or email simupop-list@lists.sourceforge.net (subscription required).", file=sys.stderr)
    # Turn on general debug information when not in 'quiet' mode
    # This will print out error messages when needed.
    turnOnDebug('DBG_GENERAL')
    # only display banner once
    simuOptions['Quiet'] = True

if simuOptions['Debug'] != []:
    for g in simuOptions['Debug']:
        if g not in ['', None]:
            print("Turn on debug '%s'" % g, file=sys.stderr)
            turnOnDebug(g)

ALL_AVAIL = True
UNSPECIFIED = False


def evolve_pop(self,
               initOps=[],
               preOps=[],
               matingScheme=MatingScheme(),
               postOps=[],
               finalOps=[],
               gen=-1,
               dryrun=False):
    '''Evolve the current population *gen* generations using mating scheme
    *matingScheme* and operators *initOps* (applied before evolution), *preOps*
    (applied to the parental population at the beginning of each life cycle),
    *postOps* (applied to the offspring population at the end of each life
    cycle) and *finalOps* (applied at the end of evolution). More specifically,
    this function creates a *Simulator* using the current population, call its
    *evolve* function using passed parameters and then replace the current
    population with the evolved population. Please refer to function
    ``Simulator.evolve`` for more details about each parameter.'''
    if dryrun:
        print(
            describeEvolProcess(initOps, preOps, matingScheme, postOps,
                                finalOps, gen, 1), file=sys.stderr)
        return (0,)
    if isinstance(self, Pedigree):
        raise ValueError("Evolving a pedigree object directly is not allowed.")
    # create a simulator with self
    simu = Simulator(self)
    # evolve
    gen = simu.evolve(initOps, preOps, matingScheme, postOps, finalOps, gen)
    # get the evolved population
    self.swap(simu.population(0))
    return gen[0]


Population.evolve = evolve_pop


def all_individuals(self, subPops=ALL_AVAIL, ancGens=ALL_AVAIL):
    '''Return an iterator that iterat through all (virtual) subpopulations in
    all ancestral generations. A list of (virtual) subpopulations (*subPops*)
    and a list of ancestral generations (*ancGens*, can be a single number)
    could be specified to iterate through only selected subpopulation and
    generations. Value ``ALL_AVAIL`` is acceptable in the specification of
    ``sp`` and/or ``vsp`` in specifying a virtual subpopulation ``(sp, vsp)``
    for the iteration through all or specific virtual subpopulation in all or
    specific subpopulations.
    '''
    if ancGens is ALL_AVAIL:
        gens = range(self.ancestralGens() + 1)
    elif hasattr(ancGens, '__iter__'):
        gens = ancGens
    else:
        gens = [ancGens]
    #
    curGen = self.curAncestralGen()
    for gen in gens:
        self.useAncestralGen(gen)
        if subPops is ALL_AVAIL:
            for ind in self.individuals():
                yield ind
        else:
            for subPop in subPops:
                if hasattr(subPop, '__iter__'):
                    if len(subPop) != 2:
                        raise ValueError('Invalid subpopulation ID %s' % subPop)
                    if subPop[0] is ALL_AVAIL:
                        # (ALL_AVAIL, ALL_AVAIL)
                        if subPop[1] is ALL_AVAIL and self.numVirtualSubPop(
                        ) > 0:
                            for sp in range(self.numSubPop()):
                                for vsp in range(self.numVirtualSubPop()):
                                    for ind in self.individuals([sp, vsp]):
                                        yield ind
                        # (ALL_AVAIL, vsp)
                        else:
                            for sp in range(self.numSubPop()):
                                for ind in self.individuals([sp, subPop[1]]):
                                    yield ind
                    else:
                        # (sp, ALL_AVAIL)
                        if subPop[1] is ALL_AVAIL and self.numVirtualSubPop(
                        ) > 0:
                            for vsp in range(self.numVirtualSubPop()):
                                for ind in self.individuals([subPop[0], vsp]):
                                    yield ind
                        # (sp, vsp)
                        else:
                            for ind in self.individuals(subPop):
                                yield ind
                else:
                    for ind in self.individuals(subPop):
                        yield ind
    self.useAncestralGen(curGen)


Population.allIndividuals = all_individuals


def as_pedigree(self,
                idField='ind_id',
                fatherField='father_id',
                motherField='mother_id'):
    '''Convert the existing population object to a pedigree. After this function
    pedigree function should magically be usable for this function.
    '''
    if isinstance(self, Pedigree):
        return
    ped = Pedigree(
        self,
        loci=ALL_AVAIL,
        infoFields=ALL_AVAIL,
        ancGens=ALL_AVAIL,
        idField=idField,
        fatherField=fatherField,
        motherField=motherField,
        stealPop=True)
    # swap ped and this object. (I do not know if this is the right thing to do)
    self.__class__, ped.__class__ = ped.__class__, self.__class__
    self.this, ped.this = ped.this, self.this


Population.asPedigree = as_pedigree


def as_population(self):
    '''Convert the existing pedigree object to a population. This function will
    behave like a regular population after this function call.'''
    if isinstance(self, Population):
        return
    pop = population(0)
    # the pedigree data has been swapped to pop
    pop.swap(self)
    # swap ped and this object. (I do not know if this is the right thing to do)
    self.__class__, pop.__class__ = pop.__class__, self.__class__
    self.this, pop.this = pop.this, self.this


Pedigree.asPopulation = as_population


# Other definitions that does not really belong to simuUtil.py
class _dw(object):

    def __init__(self, var):
        try:
            self.__dict__ = var
        except TypeError:
            raise TypeError(
                "The returned value is not a dictionary.\nNote: simu.vars() is a list so simu.dvars() is not allowed. \n    Use simu.dvars(rep) for population namespace."
            )

    def clear(self):
        self.__dict__.clear()

    def __repr__(self):
        return str(self.__dict__)


def dvars(self, *args, **kwargs):
    return _dw(self.vars(*args, **kwargs))


Population.dvars = dvars
Simulator.dvars = dvars


# expose the clone() method to Python copy module.
def _deepcopy(self, memo):
    return self.clone()


Population.__deepcopy__ = _deepcopy
Simulator.__deepcopy__ = _deepcopy
BaseOperator.__deepcopy__ = _deepcopy

Population.__deepcopy__ = _deepcopy
Simulator.__deepcopy__ = _deepcopy
BaseOperator.__deepcopy__ = _deepcopy


def ind_setInfo2(self, field, value):
    self.setInfo(value, field)


def ind_setInfo3(self, field, value):
    if field == 'this':
        self.__dict__['this'] = value
    else:
        self.setInfo(value, field)


def ind_getInfo3(self, field):
    if field == 'this':
        return self.__dict__['this']
    else:
        return self.info(field)


if sys.version_info[0] == 3:
    Individual.__setattr__ = ind_setInfo3
    Individual.__getattr__ = ind_getInfo3
else:
    Individual.__setattr__ = ind_setInfo2
    Individual.__getattr__ = Individual.info


def obj_equal(self, obj):
    return self.__cmp__(obj) == 0


if sys.version_info[0] == 3:
    Population.__eq__ = obj_equal
    Pedigree.__eq__ = obj_equal
    Simulator.__eq__ = obj_equal

##
## This is the result of an attemp to move describeEvolProcess from the C++
## level to the Python level so that it can correctly call the 'describe'
## function of a derived operator or mating scheme if they define a 'describe'
## function.
##
## This change has been temporarily reverted because many simuPOP components
## such as HomoMating, OffspringGenerator accepts such objects but cannot
## call the Python level 'describe' function correctly. Because this function
## can only describe customized pre or post mating operators and homogeneous
## mating schemes, it does not appear to be a good idea to expose (although
## still hidden) functions such as BaseOperator.applicability() to the Python
## interface.
##

## def describeEvolProcess(initOps = [], preOps = [], matingScheme = None,
##     postOps = [], finalOps = [], gen = -1, numRep = 1):
##     '''This function takes the same parameters as ``Simulator.evolve`` and
##     output a description of how an evolutionary process will be executed. It is
##     recommended that you call this function if you have any doubt on how your
##     simulation will proceed.'''
##     allDesc = [''] * numRep
##
##     # handle single inputs
##     if not hasattr(initOps, '__iter__'):
##         initOps = [initOps]
##     if not hasattr(preOps, '__iter__'):
##         preOps = [preOps]
##     if not hasattr(postOps, '__iter__'):
##         postOps = [postOps]
##     if not hasattr(finalOps, '__iter__'):
##         finalOps = [finalOps]
##
##     for curRep in range(numRep):
##         desc = ''
##         if not initOps:
##             desc += 'No operator is used to initialize Population (initOps).\n'
##         else:
##             desc += 'Apply pre-evolution operators to the initial population (initOps).\n<ul>\n'
##             for op in initOps:
##                 desc += '<li>' + op.describe(False) + ' ' + op.applicability(True, False) + '\n'
##             desc += '</ul>\n'
##         if gen < 0:
##             desc += '\nEvolve a population indefinitely until an operator determines it.\n'
##         else:
##             desc += '\nEvolve a population for %s generations\n' % gen
##         desc += '<ul>\n'
##         if not preOps:
##             desc += '<li>No operator is applied to the parental generation (preOps).\n'
##         else:
##             desc += '<li>Apply pre-mating operators to the parental generation (preOps)\n<ul>\n'
##             for op in preOps:
##                 if op.isActive(curRep, 0):
##                     desc += '<li>' + op.describe(False) + ' ' + op.applicability() + '\n'
##             desc += '</ul>\n'
##         #
##         desc += '\n<li>Populate an offspring populaton from the parental population using mating scheme ' \
##             + matingScheme.describe(False) + '\n'
##         #
##         if not postOps:
##             desc += '\n<li>No operator is applied to the offspring population (postOps).\n'
##         else:
##             desc += '\n<li>Apply post-mating operators to the offspring population (postOps).\n<ul>\n'
##             for op in postOps:
##                 if op.isActive(curRep, 0):
##                     desc += '<li>' + op.describe(False) + ' ' + op.applicability() + '\n'
##             desc += '</ul>\n'
##         desc += '</ul>\n\n'
##         #
##         if not finalOps:
##             desc += 'No operator is applied to the final population (finalOps).\n'
##         else:
##             desc += 'Apply post-evolution operators (finalOps)\n<ul>\n'
##             for op in finalOps:
##                 desc += '<li>' + op.describe(False) + ' ' + op.applicability(True, False) + '\n'
##             desc += '</ul>\n'
##         #
##         allDesc[curRep] = desc
##     #
##     reps = []
##     for curRep in range(numRep):
##         if not reps:
##             reps.append(str(curRep))
##         else:
##             if allDesc[curRep] == allDesc[curRep - 1]:
##                 reps.append(str(curRep))
##             else:
##                 desc += 'Replicate' + ' '.join(reps)
##                 desc += ':\n' + allDesc[curRep - 1] + '\n'
##                 reps = [str(curRep)]
##     #
##     desc += 'Replicate' + ' '.join(reps) + ':\n' + allDesc[-1]
##     return formatDescription(desc)


class WithArgs:
    '''This class wraps around a user-provided function and provides an
    attribute ``args`` so that simuPOP knows which parameters to send to the
    function. This is only needed if the function can not be defined with
    allowed parameters.
    '''

    def __init__(self, func, args):
        '''Return a callable object that wraps around function ``func``.
        Parameter ``args`` should be a list of parameter names.
        '''
        self.__args__ = args
        self.func = func

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)


class WithMode:
    '''This class wraps around a user-provided output string, function
    or file handle (acceptable by parameter ``output`` of operators) so
    that simuPOP knows which mode the output should be written to. For
    example, if the output of the operator is a binary compressed stream,
    ``WithMode(output, 'b')`` could be used to tell the operators to
    output bytes instead of string. This is most needed for Python 3
    because files in Python 2 accepts string even if they are opened in
    binary mode.'''

    def __init__(self, output, mode=''):
        '''Return an object that wraps around ``output`` and tells simuPOP
        to output string in ``mode``. This class currently only support
        ``mode=''`` for text mode and ``mode='b'`` for binary output.'''
        self._with_output = output
        self._with_mode = mode


# mating schemes


class SequentialParentsChooser(CombinedParentsChooser):
    '''This parent chooser chooses two parents (a father and a mother)
    sequentially from their respective sex groups. Selection is not considered.
    If all fathers (or mothers) are exhausted, this parent chooser will choose
    fathers (or mothers) from the beginning of the (virtual) subpopulation
    again.'''

    def __init__(self):
        '''Create a parent chooser that chooses two parents sequentially from a
        parental (virtual) subpopulation.'''
        CombinedParentsChooser.__init__(self,
                                        SequentialParentChooser(MALE_ONLY),
                                        SequentialParentChooser(FEMALE_ONLY))


class CloneMating(HomoMating):
    '''A homogeneous mating scheme that uses a sequential parent chooser and
    a clone offspring generator.'''

    def __init__(self,
                 numOffspring=1,
                 sexMode=None,
                 ops=CloneGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField=None):
        '''Create a clonal mating scheme that clones parents to offspring using
        a ``CloneGenoTransmitter``. Please refer to class ``OffspringGenerator``
        for parameters *ops* and *numOffspring*, and to class ``HomoMating`` for
        parameters  *subPopSize*, *subPops* and *weight*. Parameters *sexMode* and
        *selectionField* are ignored because this mating scheme does not support
        natural selection, and ``CloneGenoTransmitter`` copies sex from parents
        to offspring. Note that ``CloneGenoTransmitter`` by default also copies
        all parental information fields to offspring.
        '''
        HomoMating.__init__(
            self,
            chooser=SequentialParentChooser(),
            generator=OffspringGenerator(ops, numOffspring, RANDOM_SEX),
            subPopSize=subPopSize,
            subPops=subPops,
            weight=weight)


class RandomSelection(HomoMating):
    '''A homogeneous mating scheme that uses a random single-parent parent
    chooser with replacement, and a clone offspring generator. This mating
    scheme is usually used to simulate the basic haploid Wright-Fisher model
    but it can also be applied to diploid populations.'''

    def __init__(self,
                 numOffspring=1,
                 sexMode=None,
                 ops=CloneGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField='fitness'):
        '''Create a mating scheme that select a parent randomly and copy him or
        her to the offspring population. Please refer to class
        ``RandomParentChooser`` for parameter *selectionField*, to class
        ``OffspringGenerator`` for parameters *ops* and *numOffspring*, and to
        class ``HomoMating`` for parameters *subPopSize*, *subPops* and *weight*.
        Parameter *sexMode* is ignored because ``cloneOffspringGenerator`` copies
        sex from parents to offspring.
        '''
        HomoMating.__init__(
            self,
            chooser=RandomParentChooser(True, selectionField),
            generator=OffspringGenerator(ops, numOffspring, RANDOM_SEX),
            subPopSize=subPopSize,
            subPops=subPops,
            weight=weight)


class RandomMating(HomoMating):
    '''A homogeneous mating scheme that uses a random parents chooser with
    replacement and a Mendelian offspring generator. This mating scheme is
    widely used to simulate diploid sexual Wright-Fisher random mating.'''

    def __init__(self,
                 numOffspring=1,
                 sexMode=RANDOM_SEX,
                 ops=MendelianGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField='fitness'):
        '''Creates a random mating ssheme that selects two parents randomly and
        transmit genotypes according to Mendelian laws. Please refer to class
        ``RandomParentsChooser`` for parameter *selectionField*, to class
        ``OffspringGenerator`` for parameters *ops*, *sexMode* and
        *numOffspring*, and to class ``HomoMating`` for parameters
        *subPopSize*, *subPops* and *weight*.
        '''
        HomoMating.__init__(
            self,
            chooser=RandomParentsChooser(True, selectionField),
            generator=OffspringGenerator(ops, numOffspring, sexMode),
            subPopSize=subPopSize,
            subPops=subPops,
            weight=weight)


class MonogamousMating(HomoMating):
    '''A homogeneous mating scheme that uses a random parents chooser without
    replacement and a Mendelian offspring generator. It differs from the basic
    random mating scheme in that each parent can mate only once so there is no
    half-sibling in the population.'''

    def __init__(self,
                 numOffspring=1,
                 sexMode=RANDOM_SEX,
                 ops=MendelianGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField=None):
        '''Creates a monogamous mating scheme that selects each parent only
        once. Please refer to class ``OffspringGenerator`` for parameters
        *ops*, *sexMode* and *numOffspring*, and to class ``HomoMating`` for
        parameters *subPopSize*, *subPops* and *weight*. Parameter
        *selectionField* is ignored because this mating scheme does not
        support natural selection.
        '''
        HomoMating.__init__(
            self,
            chooser=RandomParentsChooser(replacement=False),
            generator=OffspringGenerator(ops, numOffspring, sexMode),
            subPopSize=subPopSize,
            subPops=subPops,
            weight=weight)


class PolygamousMating(HomoMating):
    '''A homogeneous mating scheme that uses a multi-spouse parents chooser
    and a Mendelian offspring generator. It differs from the basic random
    mating scheme in that each parent of sex *polySex* will have *polyNum*
    spouses.'''

    def __init__(self,
                 polySex=MALE,
                 polyNum=1,
                 numOffspring=1,
                 sexMode=RANDOM_SEX,
                 ops=MendelianGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField='fitness'):
        '''Creates a polygamous mating scheme that each parent mates with
        multiple spouses. Please refer to class ``PolyParentsChooser`` for
        parameters *polySex*, *polyNum* and *selectionField*, to class
        ``OffspringGenerator`` for parameters *ops*,  *sexMode* and
        *numOffspring*, and to class ``HomoMating`` for parameters
        *subPopSize*, *subPops* and *weight*. '''
        HomoMating.__init__(
            self,
            chooser=PolyParentsChooser(polySex, polyNum),
            generator=OffspringGenerator(ops, numOffspring, sexMode),
            subPopSize=subPopSize,
            subPops=subPops,
            weight=weight)


class HaplodiploidMating(HomoMating):
    '''A homogeneous mating scheme that uses a random parents chooser with
    replacement and a haplodiploid offspring generator. It should be used
    in a haplodiploid population where male individuals only have one set
    of homologous chromosomes.'''

    def __init__(self,
                 numOffspring=1.,
                 sexMode=RANDOM_SEX,
                 ops=HaplodiploidGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField='fitness'):
        '''Creates a mating scheme in haplodiploid populations. Please refer
        to class ``RandomParentsChooser`` for parameter *selectionField*, to
        class ``OffspringGenerator`` for parameters *ops*, *sexMode* and
        *numOffspring*, and to class ``HomoMating`` for parameters
        *subPopSize*, *subPops* and *weight*.
        '''
        HomoMating.__init__(
            self,
            chooser=RandomParentsChooser(True, selectionField),
            generator=OffspringGenerator(ops, numOffspring, sexMode),
            subPopSize=subPopSize,
            subPops=subPops,
            weight=weight)


class SelfMating(HomoMating):
    '''A homogeneous mating scheme that uses a random single-parent parent
    chooser with or without replacement (parameter *replacement*) and a
    selfing offspring generator. It is used to mimic self-fertilization
    in certain plant populations.'''

    def __init__(self,
                 replacement=True,
                 numOffspring=1,
                 sexMode=RANDOM_SEX,
                 ops=SelfingGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField='fitness'):
        '''Creates a selfing mating scheme where two homologous copies of
        parental chromosomes are transmitted to offspring according to
        Mendelian laws. Please refer to class ``RandomParentChooser`` for
        parameter *replacement* and  *selectionField*, to class
        ``OffspringGenerator`` for parameters *ops*, *sexMode* and
        *numOffspring*, and to class ``HomoMating`` for parameters
        *subPopSize*, *subPops* and *weight*. '''
        HomoMating.__init__(
            self,
            chooser=RandomParentChooser(replacement, selectionField),
            generator=OffspringGenerator(ops, numOffspring, sexMode),
            subPopSize=subPopSize,
            subPops=subPops,
            weight=weight)


class HermaphroditicMating(HomoMating):
    '''A hermaphroditic mating scheme that chooses two parents randomly
    from the population regardless of sex. The parents could be chosen
    with or without replacement (parameter *replacement*). Selfing (if
    the same parents are chosen) is allowed unless *allowSelfing* is
    set to *False* '''

    def __init__(self,
                 replacement=True,
                 allowSelfing=True,
                 numOffspring=1,
                 sexMode=RANDOM_SEX,
                 ops=MendelianGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField='fitness'):
        '''Creates a hermaphroditic mating scheme where individuals can
        serve as father or mother, or both (self-fertilization). Please
        refer to class ``CombinedParentsChooser`` for parameter *allowSelfing``,
        to ``RandomParentChooser`` for parameter *replacement* and
        *selectionField*, to class ``OffspringGenerator`` for parameters *ops*,
        *sexMode* and *numOffspring*, and to class ``HomoMating`` for parameters
        *subPopSize*, *subPops* and *weight*. '''
        HomoMating.__init__(
            self,
            chooser=CombinedParentsChooser(
                RandomParentChooser(replacement, selectionField),
                RandomParentChooser(replacement, selectionField),
                allowSelfing=allowSelfing),
            generator=OffspringGenerator(ops, numOffspring, sexMode),
            subPopSize=subPopSize,
            subPops=subPops,
            weight=weight)


##
## def consanguineousMating(infoFields = [], func = None, param = None,
##         replacement = False, numOffspring = 1.,    sexMode = RANDOM_SEX,
##         ops = MendelianGenoTransmitter(), subPopSize = [],
##         subPops = ALL_AVAIL, weight = 0, selectionField = 'fitness'):
##     '''A homogeneous mating scheme that uses an information parents chooser and
##     a Mendelian offspring generator. A function *func* should be defined to
##     locate certain types of relative to each individual and save their indexes
##     to information fields *infoFields*. This mating scheme will then choose a
##     parent randomly and then another parent from his/her relatives using their
##     saved indexes. Please refer to class ``infoParentsChooser`` for parameters
##     *infoFields*, *func*, *param* and  *selectionField*, to class
##     ``OffspringGenerator`` for parameters *ops*, *sexMode* and *numOffspring*,
##     and to class ``HomoMating`` for parameters *subPopSize*, *subPops* and
##     *weight*.
##     '''
##     return HomoMating(
##         chooser = infoParentsChooser(infoFields, func, param, selectionField),
##         generator = OffspringGenerator(ops, numOffspring, sexMode),
##         subPopSize = subPopSize,
##         subPops = subPops,
##         weight = weight)
##


class ControlledRandomMating(HomoMating):
    '''A homogeneous mating scheme that uses a random sexual parents chooser
    with replacement and a controlled offspring generator using Mendelian
    genotype transmitter. It falls back to a regular random mating scheme
    if there is no locus to control or no trajectory is defined.'''

    def __init__(self,
                 loci=[],
                 alleles=[],
                 freqFunc=None,
                 numOffspring=1,
                 sexMode=RANDOM_SEX,
                 ops=MendelianGenoTransmitter(),
                 subPopSize=[],
                 subPops=ALL_AVAIL,
                 weight=0,
                 selectionField='fitness'):
        '''Creates a random mating scheme that controls allele frequency at
        loci *loci*. At each generation, function *freqFunc* will be called to
        called to obtain intended frequencies of alleles *alleles* at loci
        *loci*. The controlled offspring generator will control the acceptance
        of offspring so that the generation reaches desired allele frequencies
        at these loci. If *loci* is empty or *freqFunc* is ``None``, this mating
        scheme works identically to a ``RandomMating scheme``. Rationals and
        applications of this mating scheme is described in details in a paper *Peng
        et al, 2007 (PLoS Genetics)*. Please refer to class ``RandomParentsChooser``
        for parameters *selectionField*, to class ``ControlledOffspringGenerator``
        for parameters *loci*, *alleles*, *freqFunc*, to class
        ``OffspringGenerator`` for parameters *ops*, *sexMode* and *numOffspring*,
        and to class ``HomoMating`` for parameters *subPopSize*, *subPops* and
        *weight*.
        '''
        if (type(loci) in [type([]), type(
            ())] and len(loci) == 0) or (freqFunc is None):
            HomoMating.__init__(
                self,
                chooser=RandomParentsChooser(True, selectionField),
                generator=OffspringGenerator(ops, numOffspring, sexMode),
                subPopSize=subPopSize,
                subPops=subPops,
                weight=weight)
        else:
            HomoMating.__init__(
                self,
                chooser=RandomParentsChooser(True, selectionField),
                generator=ControlledOffspringGenerator(loci, alleles, freqFunc,
                                                       ops, numOffspring,
                                                       sexMode),
                subPopSize=subPopSize,
                subPops=subPops,
                weight=weight)


class SNPMutator(MatrixMutator):
    '''A mutator model that assumes two alleles 0 and 1 and accepts mutation
    rate from 0 to 1, and from 1 to 0 alleles. '''

    def __init__(self,
                 u=0,
                 v=0,
                 loci=ALL_AVAIL,
                 mapIn=[],
                 mapOut=[],
                 output='',
                 begin=0,
                 end=-1,
                 step=1,
                 at=[],
                 reps=ALL_AVAIL,
                 subPops=ALL_AVAIL,
                 infoFields=['ind_id'],
                 lineageMode=FROM_INFO):
        '''Return a ``MatrixMutator`` with proper mutate matrix for a two-allele
        mutation model using mutation rate from allele 0 to 1 (parameter ``u``)
        and from 1 to 0 (parameter ``v``)'''
        MatrixMutator.__init__(self, [[1 - u, u], [v, 1 - v]], loci, mapIn,
                               mapOut, output, begin, end, step, at, reps,
                               subPops, infoFields, lineageMode)


class AcgtMutator(MatrixMutator):
    '''This mutation operator assumes alleles 0, 1, 2, 3 as nucleotides ``A``,
    ``C``, ``G`` and ``T`` and use a 4 by 4 mutation rate matrix to mutate them.
    Although a general model needs 12 parameters, less parameters are needed
    for specific nucleotide mutation models (parameter ``model``). The length
    and meaning of parameter ``rate`` is model dependent.'''

    def __init__(self,
                 rate=[],
                 model='general',
                 loci=ALL_AVAIL,
                 mapIn=[],
                 mapOut=[],
                 output='',
                 begin=0,
                 end=-1,
                 step=1,
                 at=[],
                 reps=ALL_AVAIL,
                 subPops=ALL_AVAIL,
                 infoFields=['ind_id'],
                 lineageMode=FROM_INFO):
        '''Create a mutation model that mutates between nucleotides ``A``,
        ``C``, ``G``, and ``T`` (alleles are coded in that order as 0, 1, 2
        and 3). Currently supported models are Jukes and Cantor 1969 model
        (``JC69``), Kimura's 2-parameter model (``K80``), Felsenstein 1981
        model (``F81``), Hasgawa, Kishino and Yano 1985 model (``HKY85``),
        Tamura 1992 model (``T92``), Tamura and Nei 1993 model (``TN93``),
        Generalized time reversible model (``GTR``), and a general model
        (``general``) with 12 parameters. Please refer to the simuPOP user's
        guide for detailed information about each model.
        '''
        if model == 'JC69':
            if type(rate) in [type(()), type([])]:
                if len(rate) != 1:
                    raise ValueError(
                        'A Jukes and Cantor 1969 model needs one parameter mu.')
                mu = rate[0]
            else:
                mu = rate
            m = [[0, mu / 4., mu / 4., mu / 4.], [mu / 4., 0, mu / 4., mu / 4.],
                 [mu / 4., mu / 4., 0, mu / 4.], [mu / 4., mu / 4., mu / 4., 0]]
        elif model == 'K80':
            if len(rate) != 2:
                raise ValueError(
                    'A Kimura 2-parameter model requires two parameters mu and k'
                )
            mu, k = rate
            m = [[0, mu / 4., mu * k / 4., mu / 4.],
                 [mu / 4., 0, mu / 4., mu * k / 4.],
                 [mu * k / 4., mu / 4., 0, mu / 4.],
                 [mu / 4., mu * k / 4., mu / 4., 0]]
        elif model == 'F81':
            if len(rate) != 4:
                raise ValueError(
                    'A Felsenstein 1981 model requires four parameters mu, pi_A, pi_C and pi_G'
                )
            mu, piA, piC, piG = rate
            piT = 1 - piA - piC - piG
            if piA < 0 or piA > 1 or piC < 0 or piC > 1 or \
                piG < 0 or piG > 1 or piT < 0 or piT > 1:
                raise ValueError('Basic frequencies should be between 0 and 1')
            m = [[0, mu * piC, mu * piG, mu * piT],
                 [mu * piA, 0, mu * piG, mu * piT],
                 [mu * piA, mu * piC, 0, mu * piT],
                 [mu * piA, mu * piC, mu * piG, 0]]
        elif model == 'HKY85':
            if len(rate) != 5:
                raise ValueError(
                    'A Hasegawa, Kishino and Yano 1985 model requires five parameters mu, k, pi_A, pi_C and pi_G'
                )
            mu, k, piA, piC, piG = rate
            piT = 1 - piA - piC - piG
            if piA < 0 or piA > 1 or piC < 0 or piC > 1 or \
                piG < 0 or piG > 1 or piT < 0 or piT > 1:
                raise ValueError('Basic frequencies should be between 0 and 1')
            m = [[0, mu * piC, mu * k * piG, mu * piT],
                 [mu * piA, 0, mu * piG, mu * k * piT],
                 [mu * k * piA, mu * piC, 0, mu * piT],
                 [mu * piA, mu * k * piC, mu * piG, 0]]
        elif model == 'T92':
            if len(rate) != 2:
                raise ValueError(
                    'A Tamura 1992 model requires two parameters mu and pi_GC')
            mu, piGC = rate
            piG = piC = piGC / 2.
            piA = piT = (1 - piGC) / 2.
            if piA < 0 or piA > 1 or piC < 0 or piC > 1 or \
                piG < 0 or piG > 1 or piT < 0 or piT > 1:
                raise ValueError('Basic frequencies should be between 0 and 1')
            m = [[0, mu * piC, mu * piG, mu * piT],
                 [mu * piA, 0, mu * piG, mu * piT],
                 [mu * piA, mu * piC, 0, mu * piT],
                 [mu * piA, mu * piC, mu * piG, 0]]
        elif model == 'TN93':
            if len(rate) != 6:
                raise ValueError(
                    'A Tamura and Nei 1993 model requires six parameters mu, k1, k2, pi_A, pi_C and pi_G'
                )
            mu, k1, k2, piA, piC, piG = rate
            piT = 1 - piA - piC - piG
            if piA < 0 or piA > 1 or piC < 0 or piC > 1 or \
                piG < 0 or piG > 1 or piT < 0 or piT > 1:
                raise ValueError('Basic frequencies should be between 0 and 1')
            m = [[0, mu * piC, mu * k1 * piG, mu * piT],
                 [mu * piA, 0, mu * piG, mu * k2 * piT],
                 [mu * k1 * piA, mu * piC, 0, mu * piT],
                 [mu * piA, mu * k2 * piC, mu * piG, 0]]
        elif model == 'GTR':
            if len(rate) != 9:
                raise ValueError(
                    'A generalized time reversible model requires nine parameters x1, ..., x6, pi_A, pi_C and pi_G'
                )
            x1, x2, x3, x4, x5, x6, piA, piC, piG = rate
            piT = 1 - piA - piC - piG
            if piA < 0 or piA > 1 or piC < 0 or piC > 1 or \
                piG < 0 or piG > 1 or piT < 0 or piT > 1:
                raise ValueError('Basic frequencies should be between 0 and 1')
            m = [[0, piA * x1 / piC, piA * x2 / piG, piA * x3 / piT],
                 [x1, 0, piC * x4 / piG, piC * x5 / piT],
                 [x2, x4, 0, piG * x6 / piT], [x3, x5, x6, 0]]
        elif model == 'general':
            if len(rate) != 12:
                raise ValueError(
                    'Please specify 12 parameters for this general nucleotide mutation model'
                )
            m = [[0, rate[0], rate[1], rate[2]], [rate[3], 0, rate[4], rate[5]],
                 [rate[6], rate[7], 0, rate[8]],
                 [rate[9], rate[10], rate[11], 0]]
        else:
            raise ValueError('Unrecognized nucleotide mutation model %s' %
                             model)
        MatrixMutator.__init__(self, m, loci, mapIn, mapOut, output, begin, end,
                               step, at, reps, subPops, infoFields, lineageMode)


class AminoAcidMutator(MatrixMutator):
    '''
    This operator has not been implemented.
    '''

    def __init__(self,
                 rate=[],
                 model='general',
                 loci=ALL_AVAIL,
                 mapIn=[],
                 mapOut=[],
                 output='',
                 begin=0,
                 end=-1,
                 step=1,
                 at=[],
                 reps=ALL_AVAIL,
                 subPops=ALL_AVAIL,
                 infoFields=['ind_id'],
                 lineageMode=FROM_INFO):
        MatrixMutator.__init__(self, rate, loci, mapIn, mapOut, output, begin,
                               end, step, at, reps, subPops, infoFields,
                               lineageMode)


#
# functions to corresponding operators
def dump(pop, *args, **kwargs):
    'Apply operator ``Dumper`` to population *pop*.'
    Dumper(*args, **kwargs).apply(pop)


def initSex(pop, *args, **kwargs):
    'Apply operator ``InitSex`` to population *pop*.'
    InitSex(*args, **kwargs).apply(pop)


def initInfo(pop, *args, **kwargs):
    'Apply operator ``InitInfo`` to population *pop*.'
    InitInfo(*args, **kwargs).apply(pop)


def initGenotype(pop, *args, **kwargs):
    'Apply operator ``InitGenotype`` to population *pop*.'
    InitGenotype(*args, **kwargs).apply(pop)


def initLineage(pop, *args, **kwargs):
    'Apply operator ``InitLineage`` to population *pop*.'
    InitLineage(*args, **kwargs).apply(pop)


def pyEval(pop, *args, **kwargs):
    '''Evaluate statements *stmts* (optional) and expression *expr* in
    population *pop*\ 's local namespace and return the result of *expr*.
    If *exposePop* is given, population *pop* will be exposed in its local
    namespace as a variable with a name specified by *exposePop*. Unlike its
    operator counterpart, this function returns the result of *expr* rather
    than writting it to an output.
    '''
    return PyEval(*args, **kwargs).evaluate(pop)


def pyExec(pop, *args, **kwargs):
    '''Execute *stmts* in population *pop*\ 's local namespace.'''
    PyExec(*args, **kwargs).apply(pop)


def infoEval(pop, *args, **kwargs):
    '''Evaluate *expr* for each individual, using information fields as variables.
    Please refer to operator ``InfoEval`` for details.
    '''
    InfoEval(*args, **kwargs).apply(pop)


def infoExec(pop, *args, **kwargs):
    '''Execute *stmts* for each individual, using information fields as variables.
    Please refer to operator ``InfoExec`` for details.
    '''
    InfoExec(*args, **kwargs).apply(pop)


def migrate(pop, *args, **kwargs):
    'Function form of operator ``Migrator``.'
    Migrator(*args, **kwargs).apply(pop)


def backwardMigrate(pop, *args, **kwargs):
    'Function form of operator ``BackwardMigrator``.'
    BackwardMigrator(*args, **kwargs).apply(pop)


def splitSubPops(pop, *args, **kwargs):
    '''Split subpopulations (*subPops*) of population *pop* according to either
    *sizes* or *proportions* of the resulting subpopulations, or an information
    field. Please refer to the operator form of this function (``splitSubPop``)
    for details.'''
    SplitSubPops(*args, **kwargs).apply(pop)


def mergeSubPops(pop, *args, **kwargs):
    '''Merge subpopulations *subPops* of population *pop* into a single
    subpopulation. Please refer to the operator form of this funciton
    (``MergeSubPops``) for details'''
    MergeSubPops(*args, **kwargs).apply(pop)


def resizeSubPops(pop, *args, **kwargs):
    '''Resize subpopulations *subPops* of population *pop* into new sizes
    *size*. Individuals will be added or removed accordingly. Please refer to
    the operator form of this funciton (``ResizeSubPops``) for details'''
    ResizeSubPops(*args, **kwargs).apply(pop)


def matrixMutate(pop, *args, **kwargs):
    'Function form of operator ``MatrixMutator``'
    MatrixMutator(*args, **kwargs).apply(pop)


def snpMutate(pop, *args, **kwargs):
    'Function form of operator ``SNPMutator``'
    SNPMutator(*args, **kwargs).apply(pop)


def acgtMutate(pop, *args, **kwargs):
    'Function form of operator ``AcgtMutator``'
    AcgtMutator(*args, **kwargs).apply(pop)


def kAlleleMutate(pop, *args, **kwargs):
    'Function form of operator ``KAlleleMutator``'
    KAlleleMutator(*args, **kwargs).apply(pop)


def stepwiseMutate(pop, *args, **kwargs):
    'Function form of operator ``StepwiseMutator``'
    StepwiseMutator(*args, **kwargs).apply(pop)


def pyMutate(pop, *args, **kwargs):
    'Function form of operator ``PyMutator``'
    PyMutator(*args, **kwargs).apply(pop)


def mixedMutate(pop, *args, **kwargs):
    'Function form of operator ``MixedMutator``'
    MixedMutator(*args, **kwargs).apply(pop)


def contextMutate(pop, *args, **kwargs):
    'Function form of operator ``ContextMutator``'
    ContextMutator(*args, **kwargs).apply(pop)


def pointMutate(pop, *args, **kwargs):
    'Function form of operator ``PointMutator``'
    PointMutator(*args, **kwargs).apply(pop)


def revertFixedSites(pop, *args, **kwargs):
    'Function form of operator ``RevertFixedSites``'
    RevertFixedSites(*args, **kwargs).apply(pop)


def stat(pop, *args, **kwargs):
    '''Apply operator ``Stat`` with specified parameters to population *pop*.
    Resulting statistics could be accessed from the local namespace of ``pop``
    using functions ``pop.vars()`` or ``pop.dvars()``'''
    Stat(*args, **kwargs).apply(pop)


def tagID(pop, reset=False, *args, **kwargs):
    '''Apply operator ``IdTagger`` to population *pop* to assign a unique ID
    to all individuals in the population. Individuals ID will starts from a
    system wide index. You can reset this start ID using parameter ``reset``
    which can be ``True`` (reset to 1) or a non-negative number (start from
    this number).'''
    if reset != False:
        # True is 1.
        IdTagger().reset(reset)
    IdTagger(*args, **kwargs).apply(pop)


def mapPenetrance(pop, loci, penetrance, ancGens=ALL_AVAIL, *args, **kwargs):
    '''Apply opertor ``MapPenetrance`` to population *pop*. Unlike the
    operator form of this operator that only handles the current generation,
    this function by default assign affection status to all generations.'''
    MapPenetrance(loci, penetrance, ancGens, *args, **kwargs).apply(pop)


def maPenetrance(pop,
                 loci,
                 penetrance,
                 wildtype=0,
                 ancGens=ALL_AVAIL,
                 *args,
                 **kwargs):
    '''Apply opertor ``MaPenetrance`` to population *pop*. Unlike the
    operator form of this operator that only handles the current generation,
    this function by default assign affection status to all generations.'''
    MaPenetrance(loci, penetrance, wildtype, ancGens, *args,
                 **kwargs).apply(pop)


def mlPenetrance(pop, ops, mode, ancGens=ALL_AVAIL, *args, **kwargs):
    '''Apply opertor ``MapPenetrance`` to population *pop*. Unlike the
    operator form of this operator that only handles the current generation,
    this function by default assign affection status to all generations.'''
    MlPenetrance(ops, mode, ancGens, *args, **kwargs).apply(pop)


def pyPenetrance(pop, func, loci=[], ancGens=ALL_AVAIL, *args, **kwargs):
    '''Apply opertor ``PyPenetrance`` to population *pop*. Unlike the
    operator form of this operator that only handles the current generation,
    this function by default assign affection status to all generations.'''
    PyPenetrance(func, loci, ancGens, *args, **kwargs).apply(pop)


def pyMlPenetrance(pop,
                   func,
                   mode,
                   loci=[],
                   ancGens=ALL_AVAIL,
                   *args,
                   **kwargs):
    '''Apply opertor ``PyMlPenetrance`` to population *pop*. Unlike the
    operator form of this operator that only handles the current generation,
    this function by default assign affection status to all generations.'''
    PyMlPenetrance(func, loci, mode, ancGens, *args, **kwargs).apply(pop)


def mapSelect(pop, loci, fitness, *args, **kwargs):
    '''Apply opertor ``MapSelector`` to population *pop*.'''
    MapSelector(loci, fitness, *args, **kwargs).apply(pop)


def maSelect(pop, loci, fitness, wildtype=0, *args, **kwargs):
    '''Apply opertor ``MaSelector`` to population *pop*. '''
    MaSelector(loci, fitness, wildtype, *args, **kwargs).apply(pop)


def mlSelect(pop, ops, mode=MULTIPLICATIVE, *args, **kwargs):
    '''Apply opertor ``MlSelector`` to population *pop*.'''
    MlSelector(ops, mode, *args, **kwargs).apply(pop)


def pySelect(pop, func, loci=[], *args, **kwargs):
    '''Apply opertor ``PySelector`` to population *pop*.'''
    PySelector(func, loci, *args, **kwargs).apply(pop)


def pyMlSelect(pop, func, mode=EXPONENTIAL, loci=[], *args, **kwargs):
    '''Apply opertor ``PyMlSelector`` to population *pop*.'''
    PyMlSelector(func, mode, loci, *args, **kwargs).apply(pop)


def pyQuanTrait(pop, func, loci=[], ancGens=ALL_AVAIL, *args, **kwargs):
    '''Apply opertor ``PyQuanTrait`` to population *pop*. Unlike the
    operator form of this operator that only handles the current generation,
    this function by default assign affection status to all generations.'''
    PyQuanTrait(func, loci, ancGens, *args, **kwargs).apply(pop)


def discardIf(pop, *args, **kwargs):
    '''Apply operator ``DiscardIf`` to population *pop* to remove individuals according
    to an expression or a Python function.'''
    DiscardIf(*args, **kwargs).apply(pop)


def setRNG(name='', seed=0):
    '''Set random number generator. This function is obsolete but is provided
    for compatibility purposes. Please use setOptions instead'''
    setOptions(name=name, seed=seed)
