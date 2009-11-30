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
simuPOP sampling.

This module provides functions and operators to draw samples
from a (multi-generational) simuPOP population.

"""

__all__ = [
    # Classes that can be derived to implement more complicated
    # sampling scheme
    'randomSample',
    'affectedSibpairSample',
    'caseControlSample',
    # Functions to draw samples
    'RandomSample',
    'AffectedSibpairSample',
    'CaseControlSample',
    #
]

import exceptions, operator, types

from simuOpt import simuOptions

# avoid duplicated banner message
q = simuOptions['Quiet']
simuOptions['Quiet'] = True
from simuPOP import pyOperator, pedigree, Offspring, Spouse, affectionSplitter, Stat
simuOptions['Quiet'] = q

# Ascertainment operators and functions

class _sample:
    '''
    Ascertainment/sampling refers to ways of selecting individuals from a
    population. This base class defines the common interface of all
    ascertainment operators, including how samples are saved and returned.
    Individual ascertainment operators (derived class) only need to
    write *prepareSample* and *drawSample* functions.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Create an operator that draws a certain type of samples from a
        population *times* times. '''
        self.pedigree = None

    def prepareSample(self, pop):
        '''
        This function is usually used to prepare a pedigree object so that
        samples can be drawn.
        '''
        raise SystemError('Please re-implement this prepareSample function in the derived class.')
        return True

    def drawSample(self, pop):
        '''
        Draw and return a sample, using population *pop*, and *self.pedigree*
        prepared in prepareSample.
        '''
        raise SystemError('Please re-implement this drawSample function in the derived class.')
        return True

    def drawSamples(self, pop, times):
        self.prepareSample(pop)
        if times < 0:
            raise ValueError("Negative number of samples are unacceptable")
        if not self.prepareSample(pop):
            raise RuntimeError("Failed to prepare population for sampleing")
        # 
        return [self.drawSample(pop) for x in range(times)]

    def parent(self, idx, gen):
        if self.idField == '':
            return self.pedigree.ancestor(idx, gen)
        else:
            return self.pedigree.indByID(idx, gen)

    def clone(self):
        return copy.copy(self)


class randomSample(_sample):
    '''This operator draws random individuals from a population repeatedly and
    forms a number of random samples. These samples can be put in the
    population's local namespace, or save to disk files. The function form
    of this operator returns a list of samples directly.
    '''
    def __init__(self, size, *args, **kwargs):
        '''Draw *size* random samples from a population *times* times. *size* can
        be a number or a list of numbers. In the former case, individuals are
        drawn from the whole population and the samples has only one
        subpopulation. In the latter case, a given number of individuals are
        drawn from each subpopulation and the result sample has the same number
        of subpopulation as the population from which samples are drawn. The
        samples are saved in the population's local namespace if *name* or
        *nameExpr* is given, and are saved as diskfiles if *saveAs* or
        *saveAsExpr* is given.
        '''
        _sample.__init__(self, *args, **kwargs)
        self.size = size

    def prepareSample(self, pop):
        self.pedigree = pedigree(pop, idField='', fatherField='', motherField='')
        self.pedigree.addInfoFields('sample', -1)
        return True

    def drawSample(self, pop):
        import random
        if type(self.size) not in [type(()), type([])]:
            size = self.size
            if size > pop.popSize():
                print 'Warning: sample size %d is greater than population size %d.' % (size, pop.popSize())
                size = pop.popSize()
            # randomly choose self.size individuals
            values = [0] * size + [-1] * (pop.popSize() - size)
            random.shuffle(values)
            self.pedigree.setIndInfo(values, 'sample')
        else:
            for sp in range(pop.numSubPop()):
                size = self.size[sp]
                if size > pop.subPopSize(sp):
                    print 'Warning: sample size (%d) at subpopulation %d is greater than subpopulation size %d ' \
                        % (size, sp, pop.subPopSize(sp))
                values = [sp] * size + [-1] * (pop.subPopSize(sp) - size)
                random.shuffle(values)
                self.pedigree.setIndInfo(values, 'sample', sp)
        return pop.extract(field='sample', ped=self.pedigree)


def RandomSample(pop, size, times=1, *args, **kwargs):
    'Draw random sample'
    return randomSample(size=size, *args, **kwargs).drawSamples(pop, times=times)
 
class caseControlSample(_sample):
    '''
    This operator chooses random cases and controls from a population
    repeatedly. These samples can be put in the population's local namespace,
    or save to disk files. The function form of this operator returns a list
    of samples directly.
    '''
    def __init__(self, cases, controls, *args, **kwargs):
        '''
        Draw *cases* affected and *controls* unaffected individuals from a
        population repeatedly. *cases* can be a number or a list of numbers.
        In the former case, affected individuals are drawn from the whole
        population. In the latter case, a given number of individuals are
        drawn from each subpopulation. The same hold for *controls*. The
        resulting samples have two subpopulations that hold cases and controls
        respectively. The samples are saved in the population's local namespace
        if *name* or *nameExpr* is given, and are saved as diskfiles if
        *saveAs* or *saveAsExpr* is given.
        '''
        _sample.__init__(self, *args, **kwargs)
        self.cases = cases
        self.controls = controls

    def prepareSample(self, pop):
        self.pedigree = pedigree(pop, idField='', fatherField='', motherField='')
        self.pedigree.addInfoFields('sample', -1)
        self.pedigree.setVirtualSplitter(affectionSplitter())
        return True

    def drawSample(self, pop):
        import random
        if type(self.cases) not in [type(()), type([])] and type(self.controls) not in [type(()), type([])]:
            Stat(pop, numOfAffected=True)
            allCases = pop.dvars().numOfAffected
            allControls = pop.dvars().numOfUnaffected
            #
            cases = self.cases
            if cases > allCases:
                print 'Warning: number of cases %d is greater than number of affected individuals %d.' \
                    % (cases, allCases)
                cases = allCases
            #
            controls = self.controls
            if controls > allControls:
                print 'Warning: number of controls %d is greater than number of affected individuals %d.' \
                    % (controls, allControls)
                controls = allControls
            # caseControlly choose self.size individuals
            value_cases = [0] * cases + [-1] * (allCases - cases)
            value_controls = [1] * controls + [-1] * (allControls - controls)
            random.shuffle(value_cases)
            random.shuffle(value_controls)
            # assign information fields
            idx_cases = 0
            idx_controls = 0
            for sp in range(pop.numSubPop()):
                self.pedigree.setIndInfo(value_cases[idx_cases:], 'sample', (sp, 1))
                self.pedigree.setIndInfo(value_controls[idx_controls:], 'sample', (sp, 0))
                idx_cases += self.pedigree.subPopSize((sp, 1))
                idx_controls += self.pedigree.subPopSize((sp, 0))
        else:
            if len(self.cases) != pop.numSubPop():
                raise ValueError('If an list of cases is given, it should be specified for all subpopulations')
            if len(self.controls) != pop.numSubPop():
                raise ValueError('If an list of controls is given, it should be specified for all subpopulations')
            for sp in range(pop.numSubPop()):
                allCases = self.pedigree.subPopSize((sp, 1))
                allControls = self.pedigree.subPopSize((sp, 0))
                #
                cases = self.cases[sp]
                if cases > allCases:
                    print 'Warning: number of cases %d is greater than number of affected individuals %d in subpopulation %d.' \
                        % (cases, allCases, sp)
                    cases = allCases
                #
                controls = self.controls[sp]
                if controls > allControls:
                    print 'Warning: number of controls %d is greater than number of affected individuals %d in subpopulation %d.' \
                        % (controls, allControls, sp)
                    controls = allControls
                # 
                value_cases = [0] * cases + [-1] * (allCases - cases)
                value_controls = [1] * controls + [-1] * (allControls - controls)
                random.shuffle(value_cases)
                random.shuffle(value_controls)
                # assign information fields
                self.pedigree.setIndInfo(value_cases, 'sample', (sp, 1))
                self.pedigree.setIndInfo(value_controls, 'sample', (sp, 0))
        return pop.extract(field='sample', ped=self.pedigree)


def CaseControlSample(pop, times=1, *args, **kwargs):
    'Draw case controls ample'
    s = caseControlSample(*args, **kwargs)
    return s.drawSamples(pop, times)
 


class affectedSibpairSample(_sample):
    '''
    This operator chooses affected sibpairs and their parents from a population
    repeatedly. These samples can be put in the population's local namespace,
    or save to disk files. The function form of this operator returns a list
    of samples directly.\n
    
    The population to be sampled needs to have at least one ancestral
    generation. In addition, parents of each offspring is needed so information
    fields, most likely *father_idx* and *mother_idx* should be used to track
    parents in the parental generation. An during mating operator
    *parentsTagger* is designed for such a purpose. In addition, because it is
    very unlikely for two random offspring to share parents, affected sibpairs
    can only be ascertained from populations that are generated using a mating
    scheme that produes more than one offspring at each mating event.
    '''
    def __init__(self, size, idField='', fatherField='father_idx', motherField='mother_idx', *args, **kwargs):
        '''
        Draw *size* families, including two affected siblings and their parents
        from a population repeatedly. The population to be sampled must have
        at least one ancestral generation. It should also have two information
        fields specified by parameter *infoFields* (Default to ``['father_idx',
        'mother_idx']``. Parameter *size* can be a number or a list of numbers.
        In the former case, affected sibpairs are drawn from the whole
        population. In the latter case, a given number of affected sibpairs are
        drawn from each subpopulation. In both cases, affected sibpairs in the
        resulting sample form their own subpopulations (of size two). The
        samples are saved in the population's local namespace if *name* or
        *nameExpr* is given, and are saved as diskfiles if *saveAs* or
        *saveAsExpr* is given.
        '''
        _sample.__init__(self, *args, **kwargs)
        self.size = size
        self.fields = []
        if fatherField != '':
            self.fields.append(fatherField)
        if motherField != '':
            self.fields.append(motherField)
        if idField != '':
            self.fields.append(idField)
        if len(self.fields) != 2:
            raise ValueError('Two information fields that indicate indexes of parents in the parental generation is needed')

    def prepareSample(self, pop):
        if pop.ancestralGens() < 1:
            raise ValueError('No ancestral generation if found.')
        for field in self.fields:
            if field not in pop.infoFields():
                raise ValueError('Information field %s not found in population' % field)
        #
        print self.fields
        self.pedigree = pedigree(pop, fatherField=self.fields[0], motherField=self.fields[1],
            ancGen=1)
        self.pedigree.addInfoFields(['sample', 'pedindex', 'offspring0', 'offspring1', 'spouse'], -1)
        # locate all affected siblings
        self.pedigree.locateRelatives(Offspring, ['offspring0', 'offspring1'])
        self.pedigree.locateRelatives(Spouse, ['spouse'])
        # look for affected siblings from the parental generation
        self.pedigree.useAncestralGen(1)
        parent0 = self.pedigree.infoIdx(self.fields[0])
        parent1 = self.pedigree.infoIdx(self.fields[1])
        pedindex = self.pedigree.infoIdx('pedindex')
        offspring0 = self.pedigree.infoIdx('offspring0')
        offspring1 = self.pedigree.infoIdx('offspring1')
        spouse = self.pedigree.infoIdx('spouse')
        #
        pedCount = 0
        self.validPeds = [[] for x in range(self.pedigree.numSubPop())]
        for selfIdx, ind in enumerate(self.pedigree.individuals()):
            # if this individual is used
            # or if no valid spouse
            # or if no valid first offspring
            # or if no valid second offspring
            if int(ind.info(pedindex)) != -1 \
                or int(ind.info(spouse)) == -1 \
                or int(ind.info(offspring0)) == -1 \
                or int(ind.info(offspring1)) == -1:
                continue
            # if spouse has been used
            spouseIdx = int(ind.info(spouse))
            spouseInd = self.pedigree.individual(spouseIdx)
            if int(spouseInd.info(pedindex)) != -1:
                continue
            # if the first offspring has been used, or if parents do not match, or if
            # not affected.
            offspring0Ind = self.pedigree.ancestor(int(ind.info(offspring0)), 0)
            if not offspring0Ind.affected() \
                or int(offspring0Ind.info(pedindex)) != -1 \
                or int(offspring0Ind.info(parent0)) not in [selfIdx, spouseIdx] \
                or int(offspring0Ind.info(parent1)) not in [selfIdx, spouseIdx]:
                continue
            # if the second offspring has been used, or if parents do not match, or
            # if not affected
            offspring1Ind = self.pedigree.ancestor(int(ind.info(offspring1)), 0)
            if not offspring1Ind.affected() \
                or int(offspring1Ind.info(pedindex)) != -1 \
                or int(offspring1Ind.info(parent1)) not in [selfIdx, spouseIdx] \
                or int(offspring1Ind.info(parent1)) not in [selfIdx, spouseIdx]:
                continue
            # good pedigree
            ind.setInfo(pedCount, pedindex)
            spouseInd.setInfo(pedCount, pedindex)
            offspring0Ind.setInfo(pedCount, pedindex)
            offspring1Ind.setInfo(pedCount, pedindex)
            # count the number of pedigrees
            self.validPeds[self.pedigree.subPopIndPair(selfIdx)[0]].append(pedCount)
            pedCount += 1
        return True

    def drawSample(self, pop):
        #
        import random
        pedindex = self.pedigree.infoIdx('pedindex')
        sample = self.pedigree.infoIdx('sample')
        #
        # clear information sample in case this operator is applied twice
        self.pedigree.setIndInfo([0], sample)
        #
        pedCount = sum([len(x) for x in self.validPeds])
        chosenPeds = [False] * pedCount
        #
        if type(self.size) not in [type(()), type([])]:
            size = self.size
            if size > pedCount:
                print 'Warning: number of requested sibpairs %d is greater than what exists (%d).' \
                    % (size, pedCount)
                size = pedCount
            #
            values = range(pedCount)
            random.shuffle(values)
            for v in values[:size]:
                chosenPeds[v] = True
        else:
            if len(self.size) != pop.numSubPop():
                raise ValueError('If an list of sizes is given, it should be specified for all subpopulations')
            for sp in range(pop.numSubPop()):
                allPeds = len(self.validPeds[sp])
                #
                size = self.size[sp]
                if size > allPeds:
                    print 'Warning: number of requested sibpairs %d is greater than what exists (%d) in subpopulation %d.' \
                        % (size, allPeds, sp)
                    size = allPeds
                #
                random.shuffle(self.validPeds[sp])
                for v in self.validPeds[sp][:size]:
                    chosenPeds[v] = True
        # assign genotype
        for gen in range(1, -1, -1):
            self.pedigree.useAncestralGen(gen)
            for ind in self.pedigree.individuals():
                ped = int(ind.info(pedindex))
                if ped != -1 and chosenPeds[ped]:
                    ind.setInfo(ped, sample)
                else:
                    ind.setInfo(-1, sample)
        sample = pop.extract(field='sample', ancGen=1, ped=self.pedigree)
        for gen in range(1, -1, -1):
            sample.useAncestralGen(gen)
            sample.removeSubPops([x for x in range(sample.numSubPop()) if sample.subPopSize(x) == 0])
        return sample


def AffectedSibpairSample(pop, size, times=1, *args, **kwargs):
    'Draw affected sibpair sample from ...'
    s = affectedSibpairSample(size, *args, **kwargs)
    return s.drawSamples(pop, times)
 


