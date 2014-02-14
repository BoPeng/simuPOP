#!/usr/bin/env python

#
# $File: utils.py $
# $LastChangedDate: 2014-02-05 14:38:36 -0600 (Wed, 05 Feb 2014) $
# $Rev: 4792 $
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
simuPOP demographic models

This module provides some commonly used  demographic models

"""

__all__ = [
    'migrIslandRates',
    'migrHierarchicalIslandRates',
    'migrSteppingStoneRates',
    'migr2DSteppingStoneRates',
    'ExponentialGrowthModel',
    'LinearGrowthModel',
    'InstantChangeModel',
    'MultiStageModel',
    'OutOfAfricaModel',
    'SettlementOfNewWorldModel',
    'printDemographicModel',
    'plotDemographicModel',
]

import sys
import time
import math

from simuOpt import simuOptions

from simuPOP import Population, PyEval, RandomMating, \
    ALL_AVAIL, Stat, stat, Migrator, InitSex, PyOperator

from simuPOP.utils import migrIslandRates, migrHierarchicalIslandRates, \
    migrSteppingStoneRates

from collections import OrderedDict

# the following lecture provides a good review of demographic models
#
# http://www.stats.ox.ac.uk/~mcvean/L5notes.pdf

try:
    import numpy as np
    import matplotlib.pylab as plt
    has_plotter = True
except ImportError:
    has_plotter = False

def migr2DSteppingStoneRates(r, m, n, diagonal=False, circular=False):
    '''migration rate matrix for 2D stepping stone model, with or without
    diagonal neighbors (4 or 8 neighbors for central patches). The boundaries
    are connected if circular is True. Otherwise individuals from corner and
    bounary patches will migrate to their neighbors with higher probability.
    '''
    if n < 2 and n < 2:
        return [[1]]
    rates = []
    n = int(n)
    m = int(m)
    for row in range(m):
        for col in range(n):
            if diagonal:
                neighbors = [[row-1, col], [row+1, col], [row, col-1], [row, col+1],
                    [row-1, col-1], [row-1, col+1], [row+1, col-1], [row+1, col+1]]
            else:
                neighbors = [[row-1, col], [row+1, col], [row, col-1], [row, col+1]]
            #
            if circular:
                # -1 will become n-1, n+1 will become 1
                neighbors = [(x[0] % m, x[1] % n) for x in neighbors]
            else:
                # out of boundary patches are removed
                neighbors = [(x[0], x[1]) for x in neighbors if x[0] >= 0 and x[0] < m and x[1] >= 0 and x[1] < n]
            #
            # the neighbors might overlap or cover the cell if the dimension is small
            neighbors = set(neighbors) - set([(row, col)])
            # itself
            rates.append([0]*(m*n))
            rates[-1][row * n + col] = 1. - r
            for x in neighbors:
                rates[-1][x[0] * n + x[1]] = r * 1.0 / len(neighbors)
    return rates


class BaseDemographicModel:
    '''This class is the base class for all demographic models and 
    provides common interface and utility functions for derived classes.'''
    def __init__(self, init_size=[], info_fields=[], num_gens=-1):
        '''Set attributes ``init_size``, ``info_fields``, and ``num_gens``
        to a demographic model. If size for a subpopulation for
        initial population is a list, the initial subpopulation will be split to
        multiple subpopulations. For example, ``N0=[A, [B,C]]`` is a 3-subpop
        model where the last two subpopulation will be split (and resize if needed)
        from the second subpopulation of the initial subpopulation. In addition,
        whenever a population size is allowed, a tuple of ``(size, name)`` is
        acceptable, which assigns name to the corresponding subpopulation.'''
        #
        self._raw_init_size = init_size
        self.init_size = self._extractSize(init_size)
        #
        if isinstance(info_fields, (list, tuple)):
            self.info_fields = info_fields
        else:
            self.info_fields = [info_fields]
        self.num_gens = num_gens

    def _reset(self):
        if hasattr(self, '_start_gen'):
            del self._start_gen

    def _isNamedSize(self, x):
        return isinstance(x, tuple) and len(x) == 2 and isinstance(x[1], str) and isinstance(x[0], int)

    def _extractSize(self, sz):
        # sz = 100
        if isinstance(sz, int):
            return [sz]
        elif self._isNamedSize(sz):
            return sz[0]
        res = []
        for x in sz:
            # sz = [100, 200]
            if isinstance(x, int):
                res.append(x)
            # sz = [(100, 'name')]
            elif self._isNamedSize(x):
                res.append(x[0])
            # a group 
            # sz = [[100, 200], 300]
            # sz = [[(100, 'AS'), 200], 300]
            elif isinstance(x, (tuple, list)):
                # sz = [(100, 'AS'), (200, 'ZX')]
                for y in x:
                    if isinstance(y, int):
                        res.append(y)
                    elif self._isNamedSize(y):
                        res.append(y[0])
                    else:
                        raise ValueError('Unacceptable population size: %s' % sz)
            else:
                raise ValueError('Unacceptable population size: %s' % sz)
        return res

    def _convertToNamedSize(self, sz):
        # sz = 100
        if isinstance(sz, int):
            return [(sz, '')]
        elif self._isNamedSize(sz):
            return [sz]
        res = []
        for x in sz:
            # sz = [100, 200]
            if isinstance(x, int):
                res.append((x, ''))
            # sz = [(100, 'name')]
            elif self._isNamedSize(x):
                res.append(x)
            # a group 
            # sz = [[100, 200], 300]
            # sz = [[(100, 'AS'), 200], 300]
            elif isinstance(x, (tuple, list)):
                res.append([])
                # sz = [(100, 'AS'), (200, 'ZX')]
                for y in x:
                    if isinstance(y, int):
                        res[-1].append((y, ''))
                    elif self._isNamedSize(y):
                        res[-1].append(y)
                    else:
                        raise ValueError('Unacceptable population size: %s' % sz)
            else:
                raise ValueError('Unacceptable population size: %s' % sz)
        return res

    def _fitToSize(self, pop, size):
        '''
        Fit a population to new size, split and merge population if needed
        '''
        named_size = self._convertToNamedSize(size)
        if pop.numSubPop() > 1:
            if len(named_size) == 1:
                pop.mergeSubPops()
                if named_size[0][1] != '':
                    pop.setSubPopName(named_size[0][1], 0)
            elif len(size) != pop.numSubPop():
                raise ValueError('Number of subpopulations mismatch: %d in population '
                    '%d required for ExponentialGrowthModel.' % (pop.numSubPop(), len(size)))
            elif all([self._isNamedSize(x) for x in size]):
                pop.resize([x[0] for x in named_size], propagate=True)
                for idx, (s,n) in enumerate(named_size):
                    if n != '':
                        pop.setSubPopName(n, idx)
            else:
                # this is a more complex resize method because it can resize and split
                # if a nested list is passed
                new_size = [x[0] if self._isNamedSize(x) else sum([y[0] for y in x]) for x in named_size]
                pop.resize(new_size, propagate=True)
                #
                indexes = [i for i, x in enumerate(named_size) if not self._isNamedSize(x)]
                indexes.reverse()
                for idx in indexes:
                    names = [x[1] for x in named_size[idx]]
                    pop.splitSubPop(idx, [x[0] for x in named_size[idx]],
                        names if any([x != '' for x in names]) else [])
        else:
            if len(named_size) == 1:
                pop.resize(named_size[0][0], propagate=True)
                if named_size[0][1] != '':
                    pop.setSubPopName(named_size[0][1], 0)
            else:
                if not all([self._isNamedSize(x) for x in named_size]):
                    # cannot have nested population size in this case.
                    raise ValueError('Cannot fit population with size %s to size %s' %
                        (pop.subPopSizes(), named_size))
                # split subpopulations
                pop.resize(sum([x[0] for x in named_size]), propagate=True)
                pop.splitSubPop(0, [x[0] for x in named_size])
                for idx, (s,n) in enumerate(named_size):
                    if n != '':
                        pop.setSubPopName(n, idx)

    def _expIntepolate(self, N0, NT, T, x, T0=0):
        '''x=0, ... T-1
        '''
        if x == T-1:
            return NT
        elif x >= T:
            raise ValueError('Generation number %d out of bound (0<= t < %d is expected'
                % (x, T))
        else:
            return int(math.exp(((x+1-T0)*math.log(NT) + (T-x-1)*math.log(N0))/(T-T0)))        

    def _linearIntepolate(self, N0, NT, T, x, T0=0):
        '''x=0, ... T-1
        '''
        if x == T-1:
            return NT
        elif x >= T:
            raise ValueError('Generation number %d out of bound (0<= t < %d is expected'
                % (x, T))
        else:
            return int(((x+1-T0)*NT + (T-x-1)*N0)/(T-T0))

    def __call__(self, pop):
        # the first time this function is called
        if (not hasattr(self, '_start_gen')) or self._start_gen > pop.dvars().gen:
            self._reset()
            self._start_gen = pop.dvars().gen
            # resize populations if necessary
            self._fitToSize(pop, self._raw_init_size)
        #
        self._gen = pop.dvars().gen - self._start_gen
        #
        # _reset the demographic model when it reaches the last gen
        if self._gen + 1 == self.num_gens:
            self._reset()
        

class BaseGrowthModel(BaseDemographicModel):
    '''A base model for population growth.'''
    def __init__(self, T, N0, migr=None):
        '''An population growth model that evolves a population from size ``N0``
        to ``NT`` for ``T`` generations with rate ``r``, under either an exponential,
        linear or instant population growth model. '''
        if not isinstance(T, int):
            raise ValueError('Number of generations must be an integer')
        BaseDemographicModel.__init__(self, init_size=N0,
            info_fields=[] if migr is None else 'migrate_to', num_gens=T)
        # migration
        self.migrModel = migr

    def __call__(self, pop):
        # determines generation number internally as self.gen
        BaseDemographicModel.__call__(self, pop)
        #
        if self.migrModel is not None:
            self.migrModel.apply(pop)
        

class ExponentialGrowthModel(BaseGrowthModel):
    '''A model for exponential population growth. Four parameters
    N0, T, NT, and r are allowed but only three are needed.'''
    def __init__(self, T, N0, NT=None, r=None, migr=None):
        '''An exponential population growth model that evolves a population
        from size ``N0`` to ``NT`` for ``T`` generations with rate ``r``.
        ``N0``, ``NT`` and ``r`` can be a list of population sizes or growth
        rates for multiple subpopulations. The initial population will be
        resized to ``N0`` (split if necessary). The model will automatically 
        determine ``r`` or ``NT`` if one of them is unspecified. Optionally,
        a migrator ``m`` can be provided to migrate individuals among
        subpopulations.
        '''
        BaseGrowthModel.__init__(self, T, N0, migr)
        if r is None:
            if NT is None:
                raise ValueError('Please specify ending population size NT (or growth rate r)''')
            self.NT = self._extractSize(NT)
            if len(self.NT) != len(self.init_size):
                    raise ValueError('Number of subpopulations for initial and ending generation must agree')
        elif isinstance(r, (int, float)):
            self.NT = [int(x*((1.+r)**T)) for x in self.init_size]
        elif isinstance(r, (list, tuple)):
            if len(r) != len(self.init_size):
                raise ValueError('Please specify growth rate for each subpopulation '
                    'if multiple growth rates are specified.')
            self.NT = [int(x*((1+y)**T)) for x,y in zip(self.init_size, r)]
        else:
            raise ValueError('Unacceptable growth rate (a number or a list of numbers '
                'is expected')

    def __call__(self, pop):
        BaseGrowthModel.__call__(self, pop)
        if self._gen == self.num_gens:
            return []
        else:
            return [self._expIntepolate(n0, nt, self.num_gens, self._gen)
                for (n0, nt) in zip(self.init_size, self.NT)]


class LinearGrowthModel(BaseGrowthModel):
    '''A model for linear population growth. Four parameters
    N0, T, NT, and r are allowed but only three are needed.'''
    def __init__(self, T, N0, NT=None, r=None, migr=None):
        '''An linear population growth model that evolves a population
        from size ``N0`` to ``NT`` for ``T`` generations with ``r*N0`` 
        individuals added at each generation. ``N0``, ``NT`` and ``r``
        can be a list of population sizes or growth rates for multiple 
        subpopulations. The initial population will be
        resized to ``N0`` (split if necessary). The model will automatically 
        determine ``r`` or ``NT`` if one of them is unspecified. Optionally,
        a migrator ``m`` can be provided to migrate individuals among
        subpopulations.
        '''
        BaseGrowthModel.__init__(self, T, N0, migr)
        if r is None:
            if NT is None:
                raise ValueError('Please specify ending population size NT (or growth rate r)''')
            self.NT = self._extractSize(NT)
        elif isinstance(r, (int, float)):
            self.NT = [int(x*(1+r*T)) for x in self.init_size]
        elif isinstance(r, (list, tuple)):
            if len(r) != len(self.init_size):
                raise ValueError('Please specify growth rate for each subpopulation '
                    'if multiple growth rates are specified.')
            self.NT = [int(x*(1+y*T)) for x,y in zip(self.init_size, r)]
        else:
            raise ValueError('Unacceptable growth rate (a number or a list of numbers '
                'is expected')

    def __call__(self, pop):
        BaseGrowthModel.__call__(self, pop)
        if self._gen == self.num_gens:
            return []
        else:
            return [self._linearIntepolate(n0, nt, self.num_gens, self._gen)
                for (n0, nt) in zip(self.init_size, self.NT)]    

class InstantChangeModel(BaseGrowthModel):
    '''A model for instant population growth model.'''
    def __init__(self, T, N0, G=[], NG=[], migr=None):
        '''An instant population growth model that evolves a population
        from size ``N0`` to ``NT`` for ``T`` generations with population
        size changes at generation ``G`` to ``NT``. If ``G`` is a list,
        multiple population size changes are allowed. In that case, a list
        (or a nested list) of population size should be provided to parameter
        ``NT``. This model supports population merge to and split from a single
        population.
        '''
        BaseGrowthModel.__init__(self, T, N0, migr)
        if isinstance(G, int):
            self.G = [G]
            self.NG = [NG]
        else:
            if not isinstance(NG, (tuple, list)):
                raise ValueError('Multiple sizes should be specified if multiple G is provided.')
            if len(G) != len(NG):
                raise ValueError('Please provide population size for each growth generation.')
            self.G = G
            self.NG = NG
        #
        for g in self.G:
            if g < 0 or g >= self.num_gens:
                raise ValueError('Population change generation %d exceeds total number of generations %d' \
                    % (g, self.num_gens))
                    

    def __call__(self, pop):
        BaseGrowthModel.__call__(self, pop)
        if self._gen in self.G:
            sz = self.NG[self.G.index(self._gen)]
            self._fitToSize(pop, sz)
        return pop.subPopSizes()


class MultiStageModel(BaseDemographicModel):
    '''A multi-stage demographic model connects a few demographic models 
    '''
    def __init__(self, models):
        '''An multi-stage demographic model that connects specified
        demographic models.'''
        flds = []
        gens = []
        for x in models:
            flds.extend(x.info_fields)
            gens.append(x.num_gens)
        if all([x>=0 for x in gens]):
            total_gens = sum(gens)
        else:
            total_gens = -1
        BaseDemographicModel.__init__(self, init_size=models[0].init_size,
            info_fields=list(set(flds)), num_gens=total_gens)
        #
        self.models = models
        self._model_idx = 0

    def _reset(self):
        self._model_idx = 0
        if hasattr(self, '_start_gen'):
            del self._start_gen
        for m in self.models:
            if hasattr(m, '_start_gen'):
                del m._start_gen
  
    def __call__(self, pop):
        # determines generation number internally as self.gen
        BaseDemographicModel.__call__(self, pop)
        if self._model_idx == len(self.models):
            self._reset()
            return []
        # at the end?
        if self.models[self._model_idx].num_gens == self._gen:
            self._model_idx += 1
            self._start_gen = pop.dvars().gen
        # not at the end
        sz = self.models[self._model_idx].__call__(pop)
        if sz == []:
            # this is the end, 
            self._model_idx += 1
            self._start_gen = pop.dvars().gen
            self.__call__(pop)
        return sz

class OutOfAfricaModel(MultiStageModel):
    '''A dempgrahic model for the CHB, CEU, and YRI populations, as defined in
    Gutenkunst 2009, Plos Genetics. The model is depicted in Figure 2, and the 
    default parameters are listed in Table 1 of this paper. '''
    def __init__(self, 
        T0,
        N_A=7300,
        N_AF=12300,
        N_B=2100,
        N_EU0=1000,
        r_EU=0.004,
        N_AS0=510,
        r_AS=0.0055,
        m_AF_B=0.00025,
        m_AF_EU=0.00003,
        m_AF_AS=0.000019,
        m_EU_AS=0.000096,
        T_AF=220000//25, 
        T_B=140000//25, 
        T_EU_AS=21200//25, 
        ):
        '''Counting **backward in time**, this model evolves a population for ``T0``
        generations (required parameter). The ancient population ``A`` started at
        size ``N_A`` and expanded at ``T_AF`` generations from now, to pop ``AF``
        with size ``N_AF``. Pop ``B`` split from pop ``AF`` at ``T_B`` generations
        from now, with size ``N_B``; Pop ``AF`` remains as ``N_AF`` individuals. 
        Pop ``EU`` and  ``AS`` split from pop ``B`` at ``T_EU_AS`` generations
        from now; with size ``N_EU0`` individuals and ``N_ASO`` individuals,
        respectively. Pop ``EU`` grew exponentially with rate ``r_EU``; Pop
        ``AS`` grew exponentially with rate ``r_AS``. The ``YRI``, ``CEU`` and
        ``CHB`` samples are drawn from ``AF``, ``EU`` and ``AS`` populations
        respectively.

        This model merges all subpopulations if it is applied to a population with
        multiple subpopulation.
        '''
        #
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        MultiStageModel.__init__(self, [
            InstantChangeModel(
                T=T0-T_EU_AS,
                N0=(N_A, 'Ancestral'),
                # change population size twice, one at T_AF, one at T_B
                G=[T0-T_AF, T0-T_B],
                NG=[
                    (N_AF, 'AF'), 
                    # at T_B, split to population B from subpopulation 1
                    [(N_AF, 'AF'), (N_B, 'B')]]),
            ExponentialGrowthModel(
                T=T_EU_AS,
                N0=[(N_AF, 'AF'), 
                    # split B into EU and AS at the beginning of this
                    # exponential growth stage
                    [(N_EU0, 'EU'), (N_AS0, 'AS')]],
                r=[0, r_EU, r_AS],
                migr=Migrator(rate=[
                    [0, m_AF_EU, m_AF_AS],
                    [m_EU_AS, 0, m_AF_EU],
                    [m_AF_AS, m_AF_EU, 0]
                    ])
                ),
            ]
        )

class SettlementOfNewWorldModel(MultiStageModel):
    '''A dempgrahic model for settlement of the new world of Americans, as defined
    in Gutenkunst 2009, Plos Genetics. The model is depicted in Figure 3, and the 
    default parameters are listed in Table 2 of this paper. '''
    def __init__(self,
        T0,
        N_A=7300,
        N_AF=12300,
        N_B=2100,
        N_EU0=1500,
        r_EU=0.0023,
        N_AS0=590,
        r_AS=0.0037,
        N_MX0=800,
        r_MX=0.005,
        m_AF_B=0.00025,
        m_AF_EU=0.00003,
        m_AF_AS=0.000019,
        m_EU_AS=0.0000135,
        T_AF=220000//25, 
        T_B=140000//25, 
        T_EU_AS=26400//25, 
        T_MX=21600//25,
        f_MX=0.48):
        '''Counting **backward in time**, this model evolves a population for ``T0``
        generations. The ancient population ``A`` started at size ``N_A`` and
        expanded at ``T_AF`` generations from now, to pop ``AF`` with size ``N_AF``.
        Pop ``B`` split from pop ``AF`` at ``T_B`` generations from now, with
        size ``N_B``; Pop ``AF`` remains as ``N_AF`` individuals. Pop ``EU`` and 
        ``AS`` split from pop ``B`` at ``T_EU_AS`` generations from now; with 
        size ``N_EU0`` individuals and ``N_ASO`` individuals, respectively. Pop
        ``EU`` grew exponentially with final population size ``N_EU``; Pop
        ``AS`` grew exponentially with final populaiton size ``N_AS``. Pop ``MX``
        split from pop ``AS`` at ``T_MX`` generations from now with size ``N_MX0``,
        grew exponentially to final size ``N_MX``. Migrations are allowed between
        populations with migration rates ``m_AF_B``, ``m_EU_AS``, ``m_AF_EU``,
        and ``m_AF_AS``. At the end of the evolution, the ``AF`` and ``CHB``
        populations are removed, and the ``EU`` and ``MX`` populations are merged
        with ``f_MX`` proportion for ``MX``. The Mexican American sample could
        be sampled from the last single population.

        This model merges all subpopulations if it is applied to a population with
        multiple subpopulation.
        '''
        #
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        N_EU=int(N_EU0*math.exp(r_EU*T_EU_AS))
        N_MX=int(N_MX0*math.exp(r_MX*T_MX))
        if int(N_EU/(1.-f_MX)*f_MX) <= N_MX:
            N_EU1 = N_EU
            N_MX1 = int(N_EU/(1.-f_MX)*f_MX)
        else:
            N_EU1 = int(N_MX/f_MX*(1.-f_MX))
            N_MX1 = N_MX
        #
        MultiStageModel.__init__(self, [
            InstantChangeModel(
                # leave one generation for last admixture step
                T=T0-T_EU_AS-1,
                N0=(N_A, 'Ancestral'),
                # change population size twice, one at T_AF, one at T_B
                G=[T0-T_AF, T0-T_B],
                NG=[
                    (N_AF, 'AF'), 
                    # at T_B, split to population B from subpopulation 1
                    [(N_AF, 'AF'), (N_B, 'B')]]),
            ExponentialGrowthModel(
                T=T_EU_AS - T_MX,
                N0=[(N_AF, 'AF'), 
                    # split B into EU and AS at the beginning of this
                    # exponential growth stage
                    [(N_EU0, 'EU'), (N_AS0, 'AS')]],
                r=[0, r_EU, r_AS],
                migr=Migrator(rate=[
                    [0, m_AF_EU, m_AF_AS],
                    [m_EU_AS, 0, m_AF_EU],
                    [m_AF_AS, m_AF_EU, 0]
                    ])
                ),
            ExponentialGrowthModel(T=T_MX,
                N0=[(N_AF, 'AF'),
                    # initial population size has to be calculated
                    # because we only know the final population size of
                    # EU and AS
                    (int(N_EU0*((1.+r_EU)**(T_EU_AS-T_MX))), 'EU'),
                    # split MX from AS
                    [(int(N_AS0*((1.+r_AS)**(T_EU_AS-T_MX))), 'AS'),
                        (N_MX0, 'MX')]],
                r=[0, r_EU, r_AS, r_MX],
                migr=Migrator(rate=[
                    [0, m_AF_EU, m_AF_AS],
                    [m_EU_AS, 0, m_AF_EU],
                    [m_AF_AS, m_AF_EU, 0]
                    ],
                    # the last MX population does not involve in 
                    # migration
                    subPops=[0, 1, 2],
                    toSubPops=[0, 1, 2])
                ),
            InstantChangeModel(
                T=1,
                N0=[0, N_EU1, 0, N_MX1],
                G=[0],
                NG=[(N_EU1 + N_MX1, 'MXL')]
            )]
        )


class DemographicModelReporter:
    def __init__(self):
        pass
    
    def _reportPopSize(self, pop):
        stat(pop, popSize=True)
        if 'last_size' not in pop.vars() or pop.dvars().last_size != pop.dvars().subPopSize:
            print('%d: %s' % (pop.dvars().gen, 
                ', '.join(
                    ['%d%s' % (x, ' (%s)' % y if y else '') for x, y in zip(pop.dvars().subPopSize, pop.subPopNames())])
                ))
            pop.dvars().last_size = pop.dvars().subPopSize
        return True
 
    def outputPopSize(self, model):
        pop = Population(model.init_size, infoFields=model.info_fields)
        pop.evolve(
            preOps=[
                InitSex()
            ],
            matingScheme=RandomMating(subPopSize=model),
            postOps=PyOperator(self._reportPopSize),
            gen=model.num_gens
        )

    def _recordPopSize(self, pop):
        stat(pop, popSize=True)
        gen = pop.dvars().gen
        sz = 0
        for idx, (s, n) in enumerate(zip(pop.dvars().subPopSize, pop.subPopNames())):
            if n == '':
                n = str(idx)
            if n in self.pop_regions:
                self.pop_regions[n] = np.append(self.pop_regions[n],
                    np.array([[gen, sz, gen, sz+s]]))
            else:
                self.pop_regions[n] = np.array([gen, sz, gen, sz+s], 
                    dtype=np.uint64)
            sz += s
            pop.dvars().last_size = pop.dvars().subPopSize
        return True
    #
    def plot(self, model, filename):
        if not has_plotter:
            raise RuntimeError('This function requires module numpy and matplotlib')
        self.pop_regions = OrderedDict()
        pop = Population(model.init_size, infoFields=model.info_fields)
        pop.evolve(
            preOps=[
                InitSex()
            ],
            matingScheme=RandomMating(subPopSize=model),
            postOps=PyOperator(self._recordPopSize),
            gen=model.num_gens
        )
        # 
        fig = plt.figure()
        for name, region in self.pop_regions.items():
            region = region.reshape(region.size / 4, 4)
            points = np.append(region[:, 0:2],
                region[::-1, 2:4], axis=0)
            plt.fill(points[:,0], points[:,1], label=name)
        plt.legend(loc=2)
        plt.savefig(filename)
        plt.close()


def plotDemographicModel(model, filename):
    '''Plot the specified demographic ``model`` and save figure to 
    ``filename``. This function requires python modules ``numpy`` and
    ``matplotlib``'''
    return DemographicModelReporter().plot(model, filename)

def printDemographicModel(model):
    '''Print the population size of specified demographic ``model``'''
    return DemographicModelReporter().outputPopSize(model)

if __name__ == '__main__':
    # exponential
    print('Exponential')
    printDemographicModel(ExponentialGrowthModel(10, 100, 1000))
    printDemographicModel(ExponentialGrowthModel(10, (100, 200), (1000, 2000)))
    printDemographicModel(ExponentialGrowthModel(10, (100, 200), r=0.01))
    printDemographicModel(ExponentialGrowthModel(10, (100, 200), r=(0.01, 0.2)))
    plotDemographicModel(ExponentialGrowthModel(10, (100, 200), r=(0.01, 0.2)), 'ExpDemo.png')
    # linear
    print('Linear')
    printDemographicModel(LinearGrowthModel(10, 100, 1000))
    printDemographicModel(LinearGrowthModel(10, (100, 200), (1000, 2000)))
    printDemographicModel(LinearGrowthModel(10, (100, 200), r=0.01))
    printDemographicModel(LinearGrowthModel(10, (100, 200), r=(0.1, 0.2)))
    plotDemographicModel(LinearGrowthModel(10, (100, 200), r=(0.1, 0.2)), 'LinearDemo.png')
    # instant
    print('Instat')
    printDemographicModel(InstantChangeModel(10, 100, 5, 1000))
    printDemographicModel(InstantChangeModel(10, (100, 200), 5, (1000, 2000)))
    printDemographicModel(InstantChangeModel(10, 100, [5, 8], [500, 100]))
    plotDemographicModel(InstantChangeModel(50, 100, [5, 8, 20], [[500, 200], [100, 100], [1000, 2000]]), 'InstantDemo.png')
    #
    # multi-stage model
    print('Muti-stage')
    printDemographicModel(MultiStageModel([
        InstantChangeModel(10, 100, 5, 1000),
        ExponentialGrowthModel(20, 1000, 2000)
        ]))
    plotDemographicModel(MultiStageModel([
        InstantChangeModel(10, 100, 5, 1000),
        ExponentialGrowthModel(20, 1000, 2000)
        ]), 'MultiStageDemo.png')

    # Out Of Africa Model
    print('Out of Africa')
    printDemographicModel(OutOfAfricaModel(10000))
    plotDemographicModel(OutOfAfricaModel(10000), 'OutOfAfrica.png')
    # Settlement of New World
    #
    print('Settlement of New world')
    printDemographicModel(SettlementOfNewWorldModel(10000))
    plotDemographicModel(SettlementOfNewWorldModel(10000), 'SettlementOfNewWorld.png')
     

