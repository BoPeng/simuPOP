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
``times`` (number of samples to draw) and return a list of populations. Both
independent individuals and dependent individuals (pedigrees) are supported.

Independent individuals could be drawn from any population. Pedigree
information is not necessary and is usually ignored. Unique IDs are not needed
either although such IDs could help you identify samples in the parent
population.

Pedigrees could be drawn from multi-generational populations or age-structured
populations. All individuals are required to have a unique ID (usually tracked
by operator ``idTagger`` and are stored in information field ``ind_id``).
Parents of individuals are usually tracked by operator ``pedigreeTagger`` and
are stored in information fields ``father_id`` and ``mother_id``. If parental
information is tracked using operator ``parentsTagger`` and information fields
``father_idx`` and ``mother_idx``, a function ``sampling.indexToID`` can be
used to convert index based pedigree to ID based pedigree. Note that
``parentsTagger`` can not be used to track pedigrees in age-structured
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
    # Classes that can be derived to implement more complicated
    # sampling scheme
    'baseSampler',
    'randomSampler',
    'caseControlSampler',
    'pedigreeSampler',
    'affectedSibpairSampler',
    # Functions to draw samples
    'DrawRandomSample',
    'DrawRandomSamples',
    'DrawAffectedSibpairSample',
    'DrawAffectedSibpairSamples',
    'DrawCaseControlSample',
    'DrawCaseControlSamples',
    #
]

import exceptions, random

from simuPOP import AllAvail, Stat, pedigree, OutbredSpouse, CommonOffspring, FemaleOnly, \
    Affected

def isSequence(obj):
    return hasattr(obj, '__iter__')

# Sampling classes and functions

class baseSampler:
    '''
    A sampler extracts individuals from a simuPOP population and return them
    as separate populations. This base class defines the common interface of
    all sampling classes, including how samples prepared and returned.
    '''
    def __init__(self, subPops):
        '''Create a sampler with parameter ``subPops``, which will be used
        to prepare population for sampling.
        '''
        self.subPops = subPops
        self.pop = None

    def prepareSample(self, pop, rearrange):
        '''Prepare passed population object for sampling according to parameter
        ``subPops``. If samples are drawn from the whole population, a
        population will be trimmed if only selected (virtual) subpopulations
        are used. If samples are drawn separately from specified subpopulations,
        population ``pop`` will be rearranged (if ``rearrange==True``) so that
        each subpoulation corresponds to one element in parameter ``subPops``.
        '''
        if self.subPops == AllAvail:
            self.pop = pop
        else:
            self.pop = pop.extractSubPops(self.subPops, rearrange);
        return True

    def drawSample(self, pop):
        '''
        Draw and return a sample.
        '''
        raise SystemError('Please re-implement this drawSample function in the derived class.')

    def drawSamples(self, pop, times):
        '''
        Draw multiple samples and return a list of populations.
        '''
        if times < 0:
            raise ValueError("Negative number of samples are unacceptable")
        # 
        return [self.drawSample(pop) for x in range(times)]


class randomSampler(baseSampler):
    def __init__(self, size, subPops):
        baseSampler.__init__(self, subPops)
        self.size = size

    def drawSample(self, input_pop):
        '''Draw a random sample from passed population.
        '''
        if self.pop is None:
            # this will produce self.pop.
            self.prepareSample(input_pop, isSequence(self.size))
        #
        if not isSequence(self.size):
            size = self.size
            if size > self.pop.popSize():
                print 'Warning: sample size %d is greater than population size %d.' % (size, self.pop.popSize())
                size = pop.popSize()
            # randomly choose size individuals
            values = range(self.pop.popSize())
            random.shuffle(values)
            samples = values[:size]
        else:
            samples = []
            for sp in range(self.pop.numSubPop()):
                size = self.size[sp]
                if size > self.pop.subPopSize(sp):
                    print 'Warning: sample size (%d) at subpopulation %d is greater than subpopulation size %d ' \
                        % (size, sp, self.pop.subPopSize(sp))
                values = range(self.pop.subPopBegin(sp), self.pop.subPopEnd(sp))
                random.shuffle(values)
                samples.extend(values[:size])
        return self.pop.extractIndividuals(indexes = samples)


def DrawRandomSample(pop, size, subPops=AllAvail):
    '''Draw ``times`` random samples from a population. If a single ``size``
    is given, individuals are drawn randomly from the whole population or
    from specified (virtual) subpopulations (parameter ``subPops``). Otherwise,
    a list of numbers should be used to specify number of samples from each
    subpopulation, which can be all subpopulations if ``subPops=AllAvail``
    (default), or from each of the specified (virtual) subpopulations. The
    return value of this function is a list of populations.
    '''
    return randomSampler(size=size, subPops=subPops).drawSample(pop)


def DrawRandomSamples(pop, size, times=1, subPops=AllAvail):
    '''Draw ``times`` random samples from a population. Please refer to
    function ``DrawRandomSample for more details about parameters ``size``
    and ``subPops``.'''
    return randomSampler(size=size, subPops=subPops).drawSamples(pop, times=times)


class caseControlSampler(baseSampler):
    def __init__(self, cases, controls, subPops):
        baseSampler.__init__(self, subPops)
        self.cases = cases
        self.controls = controls
        if type(self.cases) != type(controls):
            raise exceptions.ValueError("Parameter cases and controls should have the same type.")
        if isSequence(self.cases) and isSequence(self.controls) and \
            len(self.cases) != len(self.controls):
            raise exceptions.ValueError("Cases and controls should have the same type")

    def prepareSample(self, input_pop):
        '''Find out indexes all affected and unaffected individuales.
        '''
        baseSampler.prepareSample(self, input_pop, isSequence(self.cases))
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
                print 'Warning: number of cases %d is greater than number of self.affectedected individuals %d.' \
                    % (self.cases, len(self.affected))
            #
            if self.controls > len(self.unaffected):
                print 'Warning: number of controls %d is greater than number of self.affectedected individuals %d.' \
                    % (self.controls, len(self.unaffected))
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
                    print 'Warning: number of cases %d is greater than number of self.affectedected individuals %d in subpopulation %d.' \
                        % (self.cases[sp], len(aff), sp)
                #
                if self.controls[sp] > len(unaff):
                    print 'Warning: number of controls %d is greater than number of self.affectedected individuals %d in subpopulation %d.' \
                        % (self.controls[sp], len(unaff), sp)

    def drawSample(self, input_pop):
        '''Draw a case control sample
        '''
        if self.pop is None:
            # this will produce self.pop and self.affected and self.unaffected
            self.prepareSample(input_pop)
        #
        if not isSequence(self.cases):
            random.shuffle(self.affected)
            random.shuffle(self.unaffected)
            case_indexes = self.affected[:self.cases]
            control_indexes = self.unaffected[:self.controls]
        else:
            case_indexes = []
            control_indexes = []
            for sp in range(self.pop.numSubPop()):
                random.shuffle(self.affected[sp])
                random.shuffle(self.unaffected[sp])
                case_indexes.extend(self.affected[:self.cases[sp]])
                control_indexes.extend(self.unaffected[:self.controls[sp]])
        return self.pop.extractIndividuals(indexes = case_indexes + control_indexes)


def DrawCaseControlSample(pop, cases, controls, subPops=AllAvail):
    '''Draw a case-control samples from a population with ``cases``
    affected and ``controls`` unaffected individuals. If single ``cases`` and
    ``controls`` are given, individuals are drawn randomly from the whole
    population or from specified (virtual) subpopulations (parameter
    ``subPops``). Otherwise, a list of numbers should be used to specify
    number of cases and controls from each subpopulation, which can be all
    subpopulations if ``subPops=AllAvail`` (default), or from each of the
    specified (virtual) subpopulations. The return value of this function is a
    list of populations.
    '''
    return caseControlSampler(cases, controls, subPops).drawSample(pop) 


def DrawCaseControlSamples(pop, cases, controls, times=1, subPops=AllAvail):
    '''Draw ``times`` case-control samples from a population with ``cases``
    affected and ``controls`` unaffected individuals. Please refer to function
    ``DrawCaseControlSample`` for a detailed descriptions of parameters.
    '''
    return caseControlSampler(cases, controls, subPops).drawSamples(pop, times) 


class pedigreeSampler(baseSampler):
    def __init__(self, families, subPops, idField='ind_id', fatherField='father_idx', motherField='mother_idx'):
        baseSampler.__init__(self, subPops)
        self.families = families
        self.idField = idField
        self.fatherField = fatherField
        self.motherField = motherField
        self.pedigree = None

    def prepareSample(self, pop, loci=[], infoFields=[], ancGen=-1):
        '''
        Prepare self.pedigree, some pedigree sampler might need additional loci and
        information fields for this sampler.
        '''
        # create self.pop
        baseSampler.prepareSample(self, pop, isSequence(self.families))
        # get self.pedigree
        self.pedigree = pedigree(self.pop, loci, infoFields,
            ancGen, self.idField, self.fatherField, self.motherField)


class affectedSibpairSampler(pedigreeSampler):
    def __init__(self, families, subPops, idField='ind_id', fatherField='father_idx', motherField='mother_idx'):
        pedigreeSampler.__init__(self, families, subPops, idField, fatherField, motherField)

    def prepareSample(self, input_pop):
        # this will give us self.pop and self.pedigree
        pedigreeSampler.prepareSample(self, input_pop, isSequence(self.families))
        if isSequence(self.families) and len(self.families) != self.pop.numSubPop():
            raise ValueError('If an list of family counts is given, it should be specified for all subpopulations')
        #
        # locate all affected siblings
        self.pedigree.addInfoFields(['spouse', 'off1', 'off2'])
        # only look for wife so families will not overlap
        self.pedigree.locateRelatives(OutbredSpouse, ['spouse'], FemaleOnly)
        # look for affected offspring
        self.pedigree.locateRelatives(CommonOffspring, ['spouse', 'off1', 'off2'], affectionStatus=Affected)
        # find all affected siblings
        if not isSequence(self.families):
            self.father_IDs = list(self.pedigree.individualsWithRelatives(['spouse', 'off1', 'off2']))
        else:
            self.father_IDs = []
            for sp in range(self.pedigree.numSubPop()):
                self.father_IDs.append(list(self.pedigree.individualsWithRelatives(['spouse', 'off1', 'off2'], subPops=sp)))

    def drawSample(self, input_pop):
        if self.pedigree is None:
            # this will give us self.pop, self.pedigree, and self.father_IDs
            self.prepareSample(input_pop)
        #
        if not isSequence(self.families):
            if self.families > len(self.father_IDs):
                print 'Warning: number of requested sibpairs %d is greater than what exists (%d).' \
                    % (self.families, len(self.father_IDs))
            #
            random.shuffle(self.father_IDs)
            selected_IDs = self.father_IDs[:self.families]
        else:
            selected_IDs = []
            for sp in range(self.pop.numSubPop()):
                if self.families[sp] > len(self.father_IDs[sp]):
                    print 'Warning: number of requested sibpairs %d is greater than what exists (%d) in subpopulation %d.' \
                        % (self.families[sp], len(self.father_IDs[sp]), sp)
                #
                random.shuffle(self.father_IDs[sp])
                selected_IDs.extend(self.father_IDs[sp][:self.families[sp]])
        # get father, spouse and their offspring
        IDs = []
        for id in selected_IDs:
            ind = self.pedigree.indByID(id)
            IDs.extend([id, ind.spouse, ind.off1, ind.off2])
        return self.pop.extractIndividuals(IDs = IDs, idField = self.idField)


def DrawAffectedSibpairSample(pop, families, subPops=AllAvail, 
    idField='ind_id', fatherField='father_id', motherField='mother_id'):
    '''Draw affected sibpair samples from a population. If a single
    ``families`` is given, affected sibpairs and their parents are drawn
    randomly from the whole population or from specified (virtual)
    subpopulations (parameter ``subPops``). Otherwise, a list of numbers should
    be used to specify number of families from each subpopulation, which can be
    all subpopulations if ``subPops=AllAvail`` (default), or from each of the
    specified (virtual) subpopulations. The return value of this function is
    a list of populations.
    '''
    return affectedSibpairSampler(families, subPops, idField, fatherField, motherField).drawSample(pop)
 

def DrawAffectedSibpairSamples(pop, families, times=1, subPops=AllAvail, 
    idField='ind_id', fatherField='father_id', motherField='mother_id'):
    '''
    Draw ``times`` affected sibpair samplesa from population ``pop``. Please
    refer to function ``DrawAffectedSibpairSample`` for a description of
    other parameters.
    '''
    return affectedSibpairSample(families, subPops, idField, fatherField, motherField).drawSamples(pop, times)

