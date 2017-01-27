#!/usr/bin/env python

#
# $File: sampling.py $
# $LastChangedDate: 2009-10-11 17:41:30 -0500 (Sun, 11 Oct 2009) $
# $Rev: 3033 $
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit http://simupop.sourceforge.net
# for details.
#
# Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
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
This module provides classes and functions that could be used to draw samples
from a simuPOP population. These functions accept a list of parameters such
as ``subPops`` ((virtual) subpopulations from which samples will be drawn) and
``numOfSamples`` (number of samples to draw) and return a list of populations. Both
independent individuals and dependent individuals (Pedigrees) are supported.

Independent individuals could be drawn from any Population. pedigree
information is not necessary and is usually ignored. Unique IDs are not needed
either although such IDs could help you identify samples in the parent
Population.

Pedigrees could be drawn from multi-generational populations or age-structured
populations. All individuals are required to have a unique ID (usually tracked
by operator ``IdTagger`` and are stored in information field ``ind_id``).
Parents of individuals are usually tracked by operator ``PedigreeTagger`` and
are stored in information fields ``father_id`` and ``mother_id``. If parental
information is tracked using operator ``ParentsTagger`` and information fields
``father_idx`` and ``mother_idx``, a function ``sampling.indexToID`` can be
used to convert index based pedigree to ID based Pedigree. Note that
``ParentsTagger`` can not be used to track Pedigrees in age-structured
populations because they require parents of each individual resides in a
parental generation.

All sampling functions support virtual subpopulations through parameter
``subPops``, although sample size specification might vary. This feature
allows you to draw samples with specified properties. For example, you
could select only female individuals for cases of a female-only disease,
or select individuals within certain age-range. If you specify a list
of (virtual) subpopulations, you are usually allowed to draw certain
number of individuals from each subpopulation.
"""

__all__ = [
    #
    'indexToID',
    # Classes that can be derived to implement more complicated
    # sampling scheme
    'BaseSampler',
    'RandomSampler',
    'CaseControlSampler',
    'PedigreeSampler',
    'AffectedSibpairSampler',
    'NuclearFamilySampler',
    'ThreeGenFamilySampler',
    'CombinedSampler',
    # Functions to draw samples
    'drawRandomSample',
    'drawRandomSamples',
    'drawCaseControlSample',
    'drawCaseControlSamples',
    'drawAffectedSibpairSample',
    'drawAffectedSibpairSamples',
    'drawNuclearFamilySample',
    'drawNuclearFamilySamples',
    'drawThreeGenFamilySample',
    'drawThreeGenFamilySamples',
    'drawCombinedSample',
    'drawCombinedSamples',
    #
]

from simuPOP import ALL_AVAIL, Pedigree, OUTBRED_SPOUSE, COMMON_OFFSPRING, FEMALE_ONLY, \
    MALE, AFFECTED, tagID, getRNG

import random
def random_shuffle(x):
    random.shuffle(x, getRNG().randUniform)

def isSequence(obj):
    return hasattr(obj, '__iter__')

def isNumber(obj):
    return isinstance(obj, (int, float))

def indexToID(pop, idField='ind_id', fatherField='father_id', motherField='mother_id',
              fatherIndex='father_idx', motherIndex='mother_idx', reset=False):
    '''This function adds information field idField (default to ``'ind_id'``)
    to population ``pop`` and assigns an unique ID to every individual in this
    Population. It then adds information fields fatherField (default to
    ``'fatherField'``) and motherField (default to ``'motherField'``) and set
    their values with IDs according to the established index based
    parents-children relationship. Existing information fields will be used if
    idField, fatherField or motherField already exist. This function uses a
    system-wide ID generator for unique IDs, which does not have to start from
    zero. A parameter ``reset`` could be used to reset starting ID to zero
    (if ``reset=True``) or a specified number (``reset=number``).
    '''
    pop.addInfoFields([idField, fatherField, motherField], -1)
    # set each individual's unique ID to idField
    tagID(pop, reset=reset, infoFields=idField)
    # save each individual's parents' IDs to fatherField and motherField
    for gen in range(pop.ancestralGens()-1, -1, -1):
        pop.useAncestralGen(gen)
        for ind in pop.individuals():
            if ind.info(fatherIndex) != -1:
                father = pop.ancestor(ind.info(fatherIndex), gen+1)
                ind.setInfo(father.info(idField), fatherField)
            if ind.info(motherIndex) != -1:
                mother = pop.ancestor(ind.info(motherIndex), gen+1)
                ind.setInfo(mother.info(idField), motherField)


# Sampling classes and functions

class BaseSampler:
    '''
    A sampler extracts individuals from a simuPOP population and return them
    as separate populations. This base class defines the common interface of
    all sampling classes, including how samples prepared and returned.
    '''
    def __init__(self, subPops=ALL_AVAIL):
        '''Create a sampler with parameter ``subPops``, which will be used
        to prepare population for sampling. ``subPops`` should be a list of
        (virtual) subpopulations from which samples are drawn. The default
        value is ALL_AVAIL, which means all available subpopulations of a
        Population.
        '''
        self.subPops=subPops
        self.pop = None

    def prepareSample(self, pop, rearrange):
        '''Prepare passed population object for sampling according to parameter
        ``subPops``. If samples are drawn from the whole population, a
        Population will be trimmed if only selected (virtual) subpopulations
        are used. If samples are drawn separately from specified subpopulations,
        Population ``pop`` will be rearranged (if ``rearrange==True``) so that
        each subpoulation corresponds to one element in parameter ``subPops``.
        '''
        if self.subPops == ALL_AVAIL:
            self.pop = pop
        else:
            self.pop = pop.extractSubPops(self.subPops, rearrange)
        return True

    def drawSample(self, pop):
        '''
        Draw and return a sample.
        '''
        raise SystemError('Please re-implement this drawSample function in the derived class.')

    def drawSamples(self, pop, numOfSamples):
        '''
        Draw multiple samples and return a list of populations.
        '''
        if numOfSamples < 0:
            raise ValueError("Negative number of samples are unacceptable")
        # 
        return [self.drawSample(pop) for x in range(numOfSamples)]


class RandomSampler(BaseSampler):
    '''A sampler that draws individuals randomly.
    '''
    def __init__(self, sizes, subPops=ALL_AVAIL):
        '''Creates a random sampler with specified number of individuals.
        '''
        BaseSampler.__init__(self, subPops)
        self.sizes = sizes

    def drawSample(self, input_pop):
        '''Draw a random sample from passed population.
        '''
        if self.pop is None:
            # this will produce self.pop.
            self.prepareSample(input_pop, isSequence(self.sizes))
        #
        if not isSequence(self.sizes):
            size = self.sizes
            if size > self.pop.popSize():
                print('Warning: sample size %d is greater than population size %d.' % (size, self.pop.popSize()))
                size = self.pop.popSize()
            # randomly choose size individuals
            values = list(range(self.pop.popSize()))
            random_shuffle(values)
            indexes = values[:size]
        else:
            indexes = []
            for sp in range(self.pop.numSubPop()):
                size = self.sizes[sp]
                if size > self.pop.subPopSize(sp):
                    print('Warning: sample size (%d) at subpopulation %d is greater than subpopulation size %d ' \
                        % (size, sp, self.pop.subPopSize(sp)))
                values = list(range(self.pop.subPopBegin(sp), self.pop.subPopEnd(sp)))
                random_shuffle(values)
                indexes.extend(values[:size])
        return self.pop.extractIndividuals(indexes = indexes)


def drawRandomSample(pop, sizes, subPops=ALL_AVAIL):
    '''Draw ``sizes`` random individuals from a population. If a single ``sizes``
    is given, individuals are drawn randomly from the whole population or
    from specified (virtual) subpopulations (parameter ``subPops``). Otherwise,
    a list of numbers should be used to specify number of samples from each
    subpopulation, which can be all subpopulations if ``subPops=ALL_AVAIL``
    (default), or from each of the specified (virtual) subpopulations. This
    function returns a population with all extracted individuals.
    '''
    return RandomSampler(sizes=sizes, subPops=subPops).drawSample(pop)


def drawRandomSamples(pop, sizes, numOfSamples=1, subPops=ALL_AVAIL):
    '''Draw ``numOfSamples`` random samples from a population and return a list of
    populations. Please refer to function ``drawRandomSample`` for more details
    about parameters ``sizes`` and ``subPops``.'''
    return RandomSampler(sizes=sizes, subPops=subPops).drawSamples(pop, numOfSamples=numOfSamples)


class CaseControlSampler(BaseSampler):
    '''A sampler that draws affected and unaffected individuals randomly.
    '''
    def __init__(self, cases, controls, subPops=ALL_AVAIL):
        '''Ceates a case-control sampler with specified number of cases
        and controls.
        '''
        BaseSampler.__init__(self, subPops)
        self.cases = cases
        self.controls = controls
        if type(self.cases) != type(controls):
            raise ValueError("Parameter cases and controls should have the same type.")
        if isSequence(self.cases) and isSequence(self.controls) and \
            len(self.cases) != len(self.controls):
            raise ValueError("Cases and controls should have the same type")

    def prepareSample(self, input_pop):
        '''Find out indexes all affected and unaffected individuales.
        '''
        BaseSampler.prepareSample(self, input_pop, isSequence(self.cases))
        if self.pop is None:
            # this will produce self.pop and self.cases and self.controls
            self.prepareSample(input_pop)
        #
        if not isSequence(self.cases):
            # find affected individuals
            self.affected = []
            self.unaffected = []
            for idx,ind in enumerate(self.pop.individuals()):
                if ind.affected():
                    self.affected.append(idx)
                else:
                    self.unaffected.append(idx)
            #
            if self.cases > len(self.affected):
                print('Warning: number of cases %d is greater than number of affectedected individuals %d.' \
                    % (self.cases, len(self.affected)))
            #
            if self.controls > len(self.unaffected):
                print('Warning: number of controls %d is greater than number of affectedected individuals %d.' \
                    % (self.controls, len(self.unaffected)))
        else:
            if len(self.cases) != self.pop.numSubPop():
                raise ValueError('If an list of cases is given, it should be specified for all subpopulations')
            self.affected = []
            self.unaffected = []
            for sp in range(self.pop.numSubPop()):
                # find self.affectedected individuals
                aff = []
                unaff = []
                for idx in range(self.pop.subPopBegin(sp), self.pop.subPopEnd(sp)):
                    if self.pop.individual(idx).affected():
                        aff.append(idx)
                    else:
                        unaff.append(idx)
                #
                if self.cases[sp] > len(aff):
                    print('Warning: number of cases %d is greater than number of affectedected individuals %d in subpopulation %d.' \
                        % (self.cases[sp], len(aff), sp))
                #
                if self.controls[sp] > len(unaff):
                    print('Warning: number of controls %d is greater than number of affectedected individuals %d in subpopulation %d.' \
                        % (self.controls[sp], len(unaff), sp))
                self.affected.append(aff)
                self.unaffected.append(unaff)

    def drawSample(self, input_pop):
        '''Draw a case control sample
        '''
        if self.pop is None:
            # this will produce self.pop and self.affected and self.unaffected
            self.prepareSample(input_pop)
        #
        if not isSequence(self.cases):
            random_shuffle(self.affected)
            random_shuffle(self.unaffected)
            indexes = self.affected[:self.cases] + self.unaffected[:self.controls]
        else:
            indexes = []
            for sp in range(self.pop.numSubPop()):
                random_shuffle(self.affected[sp])
                random_shuffle(self.unaffected[sp])
                indexes.extend(self.affected[sp][:self.cases[sp]])
                indexes.extend(self.unaffected[sp][:self.controls[sp]])
        return self.pop.extractIndividuals(indexes = indexes)


def drawCaseControlSample(pop, cases, controls, subPops=ALL_AVAIL):
    '''Draw a case-control samples from a population with ``cases``
    affected and ``controls`` unaffected individuals. If single ``cases`` and
    ``controls`` are given, individuals are drawn randomly from the whole
    Population or from specified (virtual) subpopulations (parameter
    ``subPops``). Otherwise, a list of numbers should be used to specify
    number of cases and controls from each subpopulation, which can be all
    subpopulations if ``subPops=ALL_AVAIL`` (default), or from each of the
    specified (virtual) subpopulations. This function returns a population with
    all extracted individuals.
    '''
    return CaseControlSampler(cases, controls, subPops).drawSample(pop) 


def drawCaseControlSamples(pop, cases, controls, numOfSamples=1, subPops=ALL_AVAIL):
    '''Draw ``numOfSamples`` case-control samples from a population with ``cases``
    affected and ``controls`` unaffected individuals and return a list of
    populations. Please refer to function ``drawCaseControlSample`` for a
    detailed descriptions of parameters.
    '''
    return CaseControlSampler(cases, controls, subPops).drawSamples(pop, numOfSamples) 


class PedigreeSampler(BaseSampler):
    '''The base class of all pedigree based sampler.
    '''
    def __init__(self, families, subPops=ALL_AVAIL, idField='ind_id',
        fatherField='father_id', motherField='mother_id'):
        '''Creates a pedigree sampler with parameters

        families
            number of families. This can be a number or a list of numbers. In
            the latter case, specified families are drawn from each
            subpopulation.

        subPops
            A list of (virtual) subpopulations from which samples are drawn.
            The default value is ALL_AVAIL, which means all available
            subpopulations of a population.
        '''
        BaseSampler.__init__(self, subPops)
        self.families = families
        self.idField = idField
        self.fatherField = fatherField
        self.motherField = motherField
        self.pedigree = None

    def prepareSample(self, pop, loci=[], infoFields=[], ancGens=ALL_AVAIL):
        '''
        Prepare self.pedigree, some pedigree sampler might need additional loci and
        information fields for this sampler.
        '''
        # create self.pop
        BaseSampler.prepareSample(self, pop, isSequence(self.families))
        # get self.pedigree
        self.pedigree = Pedigree(self.pop, loci, infoFields,
            ancGens, self.idField, self.fatherField, self.motherField)
        # Your prepareSample should define selected_IDs, which should be the
        # anchor individuals of each Pedigree, typically father or grandfather.
        # Other family members will be looked up through function self.family.
        self.selected_IDs = []

    def family(self, id):
        'Get the family of individual with id.'
        return [id]

    def drawSample(self, input_pop):
        'Randomly select Pedigrees'
        if self.pedigree is None:
            # this will give us self.pop, self.pedigree, and self.selectedIDs
            self.prepareSample(input_pop)
        #
        if not isSequence(self.families):
            if self.families > len(self.selectedIDs):
                print('Warning: number of requested Pedigrees %d is greater than what exists (%d).' \
                    % (self.families, len(self.selectedIDs)))
            # a tuple might be returned
            self.selectedIDs = list(self.selectedIDs)
            random_shuffle(self.selectedIDs)
            # we select families one by one to exclude overlapping families.
            IDs = set()
            cnt = 0
            for id in self.selectedIDs:
                fam = set(self.family(id))
                if len(IDs & fam) == 0:
                    IDs.update(fam)
                    cnt += 1
                if cnt == self.families:
                    break
            if cnt != self.families:
                print('Warning: not enough non-overlapping Pedigrees are found (requested %d, found %d).' \
                    % (self.families, cnt))
        else:
            IDs = set()
            for sp in range(self.pop.numSubPop()):
                if self.families[sp] > len(self.selectedIDs[sp]):
                    print('Warning: number of requested Pedigrees %d is greater than what exists (%d) in subpopulation %d.' \
                        % (self.families[sp], len(self.selectedIDs[sp]), sp))
                #
                # a tuple might be returned
                self.selectedIDs[sp] = list(self.selectedIDs[sp])
                random_shuffle(self.selectedIDs[sp])
                # we select families one by one to exclude overlapping families.
                cnt = 0
                for id in self.selectedIDs[sp]:
                    fam = set(self.family(id))
                    if len(IDs & fam) == 0:
                        IDs.update(fam)
                        cnt += 1
                    if cnt == self.families[sp]:
                        break
                if cnt != self.families[sp]:
                    print('Warning: not enough non-overlapping Pedigrees are found (requested %d, found %d).' \
                        % (self.families, cnt))
        # get family members
        return self.pop.extractIndividuals(IDs = list(IDs), idField = self.idField)


class AffectedSibpairSampler(PedigreeSampler):
    '''A sampler that draws a nuclear family with two affected offspring.
    '''
    def __init__(self, families, subPops=ALL_AVAIL, idField='ind_id',
        fatherField='father_id', motherField='mother_id'):
        '''Initialize an affected sibpair sampler.'''
        PedigreeSampler.__init__(self, families, subPops, idField, fatherField, motherField)

    def family(self, id):
        'Return id, its spouse and their children'
        ind = self.pedigree.indByID(id)
        return id, ind.spouse, ind.off1, ind.off2

    def prepareSample(self, input_pop):
        'Find the father or all affected sibpair families'
        # this will give us self.pop and self.pedigree
        PedigreeSampler.prepareSample(self, input_pop, isSequence(self.families))
        if isSequence(self.families) and len(self.families) != self.pop.numSubPop():
            raise ValueError('If an list of family counts is given, it should be specified for all subpopulations')
        #
        # locate all affected siblings
        self.pedigree.addInfoFields(['spouse', 'off1', 'off2'])
        # only look for wife so families will not overlap
        self.pedigree.locateRelatives(OUTBRED_SPOUSE, ['spouse'], FEMALE_ONLY)
        # look for affected offspring
        self.pedigree.locateRelatives(COMMON_OFFSPRING, ['spouse', 'off1', 'off2'], affectionStatus=AFFECTED)
        # find qualified families
        if not isSequence(self.families):
            self.selectedIDs = self.pedigree.individualsWithRelatives(['spouse', 'off1', 'off2'])
        else:
            self.selectedIDs = []
            for sp in range(self.pedigree.numSubPop()):
                self.selectedIDs.append(self.pedigree.individualsWithRelatives(['spouse', 'off1', 'off2'],
                    subPops=sp))


def drawAffectedSibpairSample(pop, families, subPops=ALL_AVAIL, 
    idField='ind_id', fatherField='father_id', motherField='mother_id'):
    '''Draw affected sibpair samples from a population. If a single
    ``families`` is given, affected sibpairs and their parents are drawn
    randomly from the whole population or from specified (virtual)
    subpopulations (parameter ``subPops``). Otherwise, a list of numbers should
    be used to specify number of families from each subpopulation, which can be
    all subpopulations if ``subPops=ALL_AVAIL`` (default), or from each of the
    specified (virtual) subpopulations. This function returns a population that
    contains extracted individuals.
    '''
    return AffectedSibpairSampler(families, subPops, idField, fatherField,
        motherField).drawSample(pop)
 

def drawAffectedSibpairSamples(pop, families, numOfSamples=1, subPops=ALL_AVAIL, 
    idField='ind_id', fatherField='father_id', motherField='mother_id'):
    '''Draw ``numOfSamples`` affected sibpair samplesa from population ``pop`` and
    return a list of populations. Please refer to function
    ``drawAffectedSibpairSample`` for a description of other parameters.
    '''
    return affectedSibpairSample(families, subPops, idField, fatherField,
        motherField).drawSamples(pop, numOfSamples)


class NuclearFamilySampler(PedigreeSampler):
    '''A sampler that draws nuclear families with specified number of affected
    parents and offspring.
    '''
    def __init__(self, families, numOffspring, affectedParents=0, affectedOffspring=0,
        subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id'):
        '''Creates a nuclear family sampler with parameters

        families
            number of families. This can be a number or a list of numbers. In the latter
            case, specified families are drawn from each subpopulation.

        numOffspring
            number of offspring. This can be a fixed number or a range [min, max].

        affectedParents
            number of affected parents. This can be a fixed number or a range [min, max].

        affectedOffspring
            number of affected offspring. This can be a fixed number of a range [min, max].

        subPops
            A list of (virtual) subpopulations from which samples are drawn.
            The default value is ALL_AVAIL, which means all available
            subpopulations of a population.
        '''
        if isNumber(numOffspring):
            if numOffspring < 1:
                raise ValueError('Number of offsprings must be equal to or larger than 1.')
            self.numOffspring = numOffspring, numOffspring
        elif isSequence(numOffspring):
            if len(numOffspring) != 2:
                raise ValueError('Two boundary numbers are needed for the allowed range of number of offsprings')
            if numOffspring[0] < 1 or numOffspring[0] > numOffspring[1]:
                raise ValueError('Minimum number of offsprings must not be smaller than 1 or larger than maximum number of offsprings.')
            self.numOffspring = numOffspring
        else:
            raise ValueError('Number of offsprings should be an integer number or a range of allowed values.')
        #
        if isNumber(affectedParents):
            if affectedParents not in [0, 1, 2]:
                raise ValueError('Number of affected individuals in parents can only take 0 or 1 or 2.')
            self.affectedParents = affectedParents, affectedParents
        elif isSequence(affectedParents):
            if len(affectedParents) != 2:
                raise ValueError('Two boundary numbers are needed for the range of number of affected parents.')
            if affectedParents[0] not in [0,1,2] or affectedParents[1] not in [0,1,2] or affectedParents[0] > affectedParents[1]:
                raise ValueError('Range of number of affected parents must be within (0, 2).')
            self.affectedParents = affectedParents
        else:
            raise ValueError('Number of affected parents should be an integer number (<= 2) or a range of allowed values.')
        #
        if isNumber(affectedOffspring):
            if affectedOffspring > self.numOffspring[1]:
                raise ValueError('Number of affected offsprings cannot be larger than number of offsprings.')
            self.affectedOffspring = affectedOffspring, affectedOffspring
        elif isSequence(affectedOffspring):
            if len(affectedOffspring) != 2:
                raise ValueError('Two boundary numbers are needed for the range of number of affected offsprings.')
            if affectedOffspring[0] > self.numOffspring[1]:
                raise ValueError('Minimum number of affected offsprings cannot be larger than number of offsprings.')
            self.affectedOffspring = affectedOffspring
        else:
            raise ValueError('Number of affected offsprings should be a proper integer nubmer or a range of allowed values.')
        #
        PedigreeSampler.__init__(self, families, subPops, idField, fatherField, motherField)

    def family(self, id):
        'Return id, its spouse and their children'
        ind = self.pedigree.indByID(id)
        offIDs = [ind.info('off%d' % x) for x in range(self.numOffspring[1])]
        return [id, ind.spouse] + [x for x in offIDs if x >= 0]

    def prepareSample(self, input_pop):
        # this will give us self.pop and self.pedigree
        PedigreeSampler.prepareSample(self, input_pop, isSequence(self.families))
        if isSequence(self.families) and len(self.families) != self.pop.numSubPop():
            raise ValueError('If an list of family counts is given, it should be specified for all subpopulations')
        #
        # locate all affected siblings
        minOffFields = ['off%d' % x for x in range(self.numOffspring[0])]
        offFields = ['off%d' % x for x in range(self.numOffspring[1])]
        self.pedigree.addInfoFields(['spouse'] + offFields)
        # only look for wife so families will not overlap
        self.pedigree.locateRelatives(OUTBRED_SPOUSE, ['spouse'], FEMALE_ONLY)
        # look for offspring
        self.pedigree.locateRelatives(COMMON_OFFSPRING, ['spouse'] + offFields)
        # check number of affected individuals and filter them out.
        def qualify(id):
            father = self.pedigree.indByID(id)
            mother = self.pedigree.indByID(father.info('spouse'))
            parAff = father.affected() + mother.affected()
            if parAff < self.affectedParents[0] or parAff > self.affectedParents[1]:
                return False
            offID = [father.info('off%d' % x) for x in range(self.numOffspring[1])]
            offAff = sum([self.pedigree.indByID(id).affected() for id in offID if id >= 0])
            if offAff < self.affectedOffspring[0] or offAff > self.affectedOffspring[1]:
                return False
            return True
        # find all families with at least minOffFields...
        if not isSequence(self.families):
            self.selectedIDs = list(filter(qualify, self.pedigree.individualsWithRelatives(['spouse'] + minOffFields)))
        else:
            self.selectedIDs = []
            for sp in range(self.pedigree.numSubPop()):
                self.selectedIDs.append(list(filter(qualify, self.pedigree.individualsWithRelatives(['spouse'] + minOffFields, subPops=sp))))


def drawNuclearFamilySample(pop, families, numOffspring, affectedParents=0,
    affectedOffspring=0, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id',
    motherField='mother_id'):
    '''Draw nuclear families from a population. Number of offspring, number of
    affected parents and number of affected offspring should be specified using
    parameters ``numOffspring``, ``affectedParents`` and ``affectedOffspring``,
    which can all be a single number, or a range ``[a, b]`` (``b`` is incldued).
    If a single ``families`` is given, Pedigrees are drawn randomly from the
    whole population or from specified (virtual) subpopulations (parameter
    ``subPops``). Otherwise, a list of numbers should be used to specify
    numbers of families from each subpopulation, which can be all
    subpopulations if ``subPops=ALL_AVAIL`` (default), or from each of the
    specified (virtual) subpopulations. This function returns a population that
    contains extracted individuals.
    '''
    return NuclearFamilySampler(families, numOffspring, affectedParents,
        affectedOffspring, subPops, idField, fatherField, motherField).drawSample(pop)
 

def drawNuclearFamilySamples(pop, families, numOffspring, affectedParents=0,
    affectedOffspring=0, numOfSamples=1, subPops=ALL_AVAIL, idField='ind_id',
    fatherField='father_id', motherField='mother_id'):
    '''Draw ``numOfSamples`` affected sibpair samplesa from population ``pop`` and
    return a list of populations. Please refer to function
    ``drawNuclearFamilySample`` for a description of other parameters.
    '''
    return nuclearFamilySample(families, numOffspring, affectedParents,
        affectedOffspring, subPops, idField, fatherField,
        motherField).drawSamples(pop, numOfSamples)


class ThreeGenFamilySampler(PedigreeSampler):
    '''A sampler that draws three-generation families with specified pedigree
    size and number of affected individuals.
    '''
    def __init__(self, families, numOffspring, pedSize, numOfAffected=0,
        subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id'):
        '''
        families
            number of families. This can be a number or a list of numbers. In the latter
            case, specified families are drawn from each subpopulation.

        numOffspring
            number of offspring. This can be a fixed number or a range [min, max].

        pedSize
            number of individuals in the Pedigree. This can be a fixed number or
            a range [min, max].

        numAfffected
            number of affected individuals in the Pedigree. This can be a fixed number
            or a range [min, max]

        subPops
            A list of (virtual) subpopulations from which samples are drawn.
            The default value is ALL_AVAIL, which means all available
            subpopulations of a population.
        '''
        if isNumber(numOffspring):
            if numOffspring < 1:
                raise ValueError('Number of offsprings must be equal to or larger than 1.')
            self.numOffspring = numOffspring, numOffspring
        elif isSequence(numOffspring):
            if len(numOffspring) != 2:
                raise ValueError('Two boundary numbers are needed for the allowed range of number of offsprings')
            if numOffspring[0] < 1 or numOffspring[0] > numOffspring[1]:
                raise ValueError('Minimum number of offsprings must not be smaller than 1 or larger than maximum number of offsprings.')
            self.numOffspring = numOffspring
        else:
            raise ValueError('Number of offsprings should be an integer number or a range of allowed values.')
        #
        if isNumber(pedSize):
            self.pedSize = pedSize, pedSize
        elif isSequence(pedSize):
            if len(pedSize) != 2:
                raise ValueError('Two boundary numbers are needed for the range of number of individuals in a Pedigree.')
            self.pedSize = pedSize
        else:
            raise ValueError('Number of affected parents should be an integer number (<= 1) or a range of allowed values.')
        #
        if isNumber(numOfAffected):
            if numOfAffected > self.pedSize[1]:
                raise ValueError('Number of affected individuals cannot be larger than pedigree size.')
            self.numOfAffected = numOfAffected, numOfAffected
        elif isSequence(numOfAffected):
            if len(numOfAffected) != 2:
                raise ValueError('Two boundary numbers are needed for the range of number of affected individuals.')
            if numOfAffected[0] > self.pedSize[1]:
                raise ValueError('Minimum number of affected offsprings cannot be larger than number of individuals in a Pedigree.')
            self.numOfAffected = numOfAffected
        else:
            raise ValueError('Number of affected offsprings should be a proper integer nubmer or a range of allowed values.')
        #
        PedigreeSampler.__init__(self, families, subPops, idField, fatherField, motherField)

    def family(self, id):
        '''Return id, its spouse, their children, children's spouse and grandchildren'''
        father = self.pedigree.indByID(id)
        spouseID = father.spouse
        offID = [father.info('off%d' % x) for x in range(self.numOffspring[1])]
        offID = [x for x in offID if x >= 1]
        offSpouseID = [self.pedigree.indByID(x).spouse for x in offID]
        offSpouseID = [x for x in offSpouseID if x >= 1]
        grandOffID = [father.info('goff%d' % x) for x in range(self.numOffspring[1]**2)]
        grandOffID = [x for x in grandOffID if x >= 1]
        return [id, spouseID] + offID + offSpouseID + grandOffID

    def prepareSample(self, input_pop):
        # this will give us self.pop and self.pedigree
        PedigreeSampler.prepareSample(self, input_pop, isSequence(self.families))
        if isSequence(self.families) and len(self.families) != self.pop.numSubPop():
            raise ValueError('If an list of family counts is given, it should be specified for all subpopulations')
        #
        # locate all affected siblings
        minOffFields = ['off%d' % x for x in range(self.numOffspring[0])]
        offFields = ['off%d' % x for x in range(self.numOffspring[1])]
        minGrandOffFields = ['goff%d' % x for x in range(self.numOffspring[0]**2)]
        grandOffFields = ['goff%d' % x for x in range(self.numOffspring[1]**2)]
        self.pedigree.addInfoFields(['spouse'] + offFields + grandOffFields)
        # only look for wife so families will not overlap
        self.pedigree.locateRelatives(OUTBRED_SPOUSE, ['spouse'], FEMALE_ONLY)
        # look for offspring
        self.pedigree.locateRelatives(COMMON_OFFSPRING, ['spouse'] + offFields)
        # look for grand children
        self.pedigree.traceRelatives(fieldPath = [offFields, offFields], resultFields = grandOffFields)
        # check number of affected individuals and filter them out.
        def qualify(id):
            father = self.pedigree.indByID(id)
            mother = self.pedigree.indByID(father.spouse)
            offID = [father.info('off%d' % x) for x in range(self.numOffspring[1])]
            offID = [x for x in offID if x >= 0]
            offSpouseID = [self.pedigree.indByID(id).spouse for x in offID]
            grandOffID = [father.info('goff%d' % x) for x in range(self.numOffspring[1]**2)]
            grandOffID = [x for x in grandOffID if x >= 0]
            pedSize = 2 + len(offID) + len(offSpouseID) + len(grandOffID)
            if pedSize < self.pedSize[0] or pedSize > self.pedSize[1]:
                return False
            # check number of affected individuals
            numAff = sum([self.pedigree.indByID(id).affected() for id in offID + offSpouseID + grandOffID]) + father.affected() + mother.affected()
            if numAff < self.numOfAffected[0] or numAff > self.numOfAffected[1]:
                return False
            return True
        # find all families with at least minOffFields...
        if not isSequence(self.families):
            self.selectedIDs = list(filter(qualify, self.pedigree.individualsWithRelatives(['spouse'] + minOffFields + minGrandOffFields)))
        else:
            self.selectedIDs = []
            for sp in range(self.pedigree.numSubPop()):
                self.selectedIDs.append(list(filter(quality, self.pedigree.individualsWithRelatives(['spouse'] + minOffFields + minGrandOffFields, subPops=sp))))


def drawThreeGenFamilySample(pop, families, numOffspring, pedSize, numOfAffected=0,
    subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id'):
    '''Draw three-generation families from a population. Such families consist
    of grant parents, their children, spouse of these children, and grand
    children. Number of offspring, total number of individuals, and total
    number of affected individuals in a pedigree should be specified using
    parameters ``numOffspring``, ``pedSize`` and ``numOfAffected``, which can all
    be a single number, or a range ``[a, b]`` (``b`` is incldued). If a single
    ``families`` is given, Pedigrees are drawn randomly from the whole
    Population or from specified (virtual) subpopulations (parameter
    ``subPops``). Otherwise, a list of numbers should be used to specify
    numbers of families from each subpopulation, which can be all
    subpopulations if ``subPops=ALL_AVAIL`` (default), or from each of the
    specified (virtual) subpopulations. This function returns a population that
    contains extracted individuals.
    '''
    return ThreeGenFamilySampler(families, numOffspring, pedSize, numOfAffected,
        subPops, idField, fatherField, motherField).drawSample(pop)
 

def drawThreeGenFamilySamples(pop, families, numOffspring, pedSize, numOfAffected=0,
    numOfSamples=1, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id',
    motherField='mother_id'):
    '''Draw ``numOfSamples`` three-generation pedigree samples from population ``pop``
    and return a list of populations. Please refer to function
    ``drawThreeGenFamilySample`` for a description of other parameters.
    '''
    return ThreeGenFamilySampler(families, numOffspring, pedSize, numOfAffected,
        subPops, idField, fatherField, motherField).drawSamples(pop, numOfSamples)


class CombinedSampler(BaseSampler):
    '''A combined sampler accepts a list of sampler objects, draw samples and
    combine the returned sample into a single population. An id field is
    required to use this sampler, which will be used to remove extra copies of
    individuals who have been drawn by different samplers.
    '''
    def __init__(self, samplers=[], idField='ind_id'):
        '''
        samplers
            A list of samplers
        '''
        BaseSampler.__init__(self)
        self.samplers = samplers
        self.idField = idField

    def drawSample(self, pop):
        allIDs = []
        self.pop = pop.clone()
        for s in self.samplers:
            sample = s.drawSample(self.pop)
            IDs = []
            for gen in range(sample.ancestralGens() + 1):
                sample.useAncestralGen(gen)
                IDs.extend(sample.indInfo(self.idField))
            allIDs.extend(IDs)
            self.pop.removeIndividuals(IDs=IDs, idField=self.idField)
        #
        # extract these guys from the original population
        return pop.extractIndividuals(IDs = allIDs, idField = self.idField)


def drawCombinedSample(pop, samplers, idField='ind_id'):
    '''Draw different types of samples using a list of ``samplers``. A
    Population consists of all individuals from these samples will
    be returned. An ``idField`` that stores an unique ID for all individuals
    is needed to remove duplicated individuals who are drawn multiple
    numOfSamples from these samplers.
    '''
    return CombinedSampler(samplers, idField=idField).drawSample(pop)

def drawCombinedSamples(pop, samplers, numOfSamples=1, idField='ind_id'):
    '''Draw combined samples ``numOfSamples`` numOfSamples and return a list of populations.
    Please refer to function ``drawCombinedSample`` for details about
    parameters ``samplers`` and ``idField``.
    '''
    return CombinedSampler(samplers, idField=idField).drawSamples(pop, numOfSamples)

