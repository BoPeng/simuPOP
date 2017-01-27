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

This module provides some commonly used demographic models. In addition
to several migration rate generation functions, it provides models that
encapsulate complete demographic features of one or more populations (
population growth, split, bottleneck, admixture, migration). These models
provides:

1. The model itself can be passed to parameter subPopSize of a mating
   scheme to determine the size of the next generation. More importantly,
   it performs necessary actions of population size change when needed.

2. The model provides attribute num_gens, which can be passed to parameter
   ``gens`` of ``Simulator.evolve`` or ``Population.evolve`` function.
   A demographic model can also terminate an evolutionary process by
   returnning an empty list so ``gens=model.num_gens`` is no longer required.

"""

__all__ = [
    'migrIslandRates',
    'migrHierarchicalIslandRates',
    'migrSteppingStoneRates',
    'migr2DSteppingStoneRates',
    # outcome-based models
    'InstantChangeModel',
    'ExponentialGrowthModel',
    'LinearGrowthModel',
    'AdmixtureModel',
    'MultiStageModel',
    # event-based models
    'EventBasedModel',
    'DemographicEvent',
    'ExpansionEvent',
    'ResizeEvent',
    'SplitEvent',
    'MergeEvent',
    'AdmixtureEvent',
    # existing models
    'OutOfAfricaModel',
    'SettlementOfNewWorldModel',
    'CosiModel'
]

import sys
import time
import math

from simuOpt import simuOptions

from simuPOP import Population, PyEval, RandomSelection, \
    ALL_AVAIL, Stat, stat, Migrator, InitSex, PyOperator, \
    MergeSubPops, SplitSubPops, ResizeSubPops

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

class DemographicModel:
    '''This class is the base class for all demographic models and 
    provides common interface and utility functions for derived classes. A
    demographic model is essentially a callable Python object that encapsulates
    required functions of a demographic model, to determine initial population
    size (``Population(size=model.init_size, infoFields=model.info_fields)``, 
    to determine size of offspring population during evolution (``subPopSize=model``
    of a mating scheme), and number of generations to evolve (``gen=model.num_gens``),
    although the first and last utility could be relaxed to for models that
    could be applied to populations with different sizes, and models that evolve
    indefinitely. '''
    def __init__(self, numGens=-1, initSize=[], ops=[], infoFields=[]):
        '''Set attributes ``init_size``, ``info_fields``, and ``num_gens``
        to a demographic model. The initial population will be merged or
        split to match ``initSize``. For example, ``N0=[A, [B,C]]`` is a 3-subpop
        model where the last two subpopulation will be split (and resize if needed)
        from the second subpopulation of the initial subpopulation (which should
        have two subpopulations). The population size can be an integer for fixed
        population size, None for the size of the population or subpopulation when
        the demographic model is first applied to, or a float number representing
        the proportion (can be larger than 1) of individuals for the whole or
        corresponding subpopulation. A ``None`` value will be assigned to
        attribute ``init_size`` in such a case because the initial population 
        size is determined dynamically. In addition, whenever a population size
        is allowed, a tuple of ``(size, name)`` is acceptable, which assigns 
        ``name`` to the corresponding subpopulation. ``numGens`` can be a
        non-negative number or ``-1``, which allows the demographic model to 
        be determinated by a terminator. One or more operators (e.g. a migration
        operator or a terminator) could be passed (parameter ``ops``) and will
        be applied to the population. The demographic model will return ``[]``
        (which will effectively terminate the evolutioonary process) if any of the
        operator returns ``False``. Information fields required by these operators
        should be passed to ``infoFields``. '''
        #
        self._raw_init_size = initSize
        self.init_size = self._extractSize(initSize)
        # for demographic model without fixed initial population size, set init_size to []
        if isinstance(self.init_size, int):
            self.init_size = [self.init_size]
        elif self.init_size is None or None in self.init_size or \
            any([isinstance(x, float) for x in self.init_size]):
            self.init_size = []
        #
        if isinstance(infoFields, (list, tuple)):
            self.info_fields = infoFields
        else:
            self.info_fields = [infoFields]
        if numGens is None:
            self.num_gens = -1
        else:
            self.num_gens = numGens
        if isinstance(ops, (tuple, list)):
            self.ops = list(ops)
        else:
            self.ops = [ops]
        #
        # history
        self.size_cache = {}

    def _reset(self):
        if hasattr(self, '_start_gen'):
            del self._start_gen

    def _isNamedSize(self, x):
        return isinstance(x, tuple) and len(x) == 2 and \
            isinstance(x[1], str) and self._isSize(x[0])

    def _isSize(self, x):
        if sys.version_info.major == 2:
            return isinstance(x, (int, float)) or x is None
        else:
            return isinstance(x, (int, float)) or x is None

    def _extractSize(self, sz):
        # sz = 100
        if self._isSize(sz):
            return [sz]
        elif self._isNamedSize(sz):
            return sz[0]
        res = []
        for x in sz:
            # sz = [100, 200]
            if self._isSize(x):
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
                    if self._isSize(y):
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
        if self._isSize(sz):
            return [(sz, '')]
        elif self._isNamedSize(sz):
            return [sz]
        res = []
        for x in sz:
            # sz = [100, 200]
            if self._isSize(x):
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
                    if self._isSize(y):
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
        # if size is None or size is [], return unchanged
        if not size:
            return 
        named_size = self._convertToNamedSize(size)
        if pop.numSubPop() > 1:
            # merge subpopualtions
            if len(named_size) == 1:
                pop.mergeSubPops()
                if named_size[0][1] != '':
                    pop.setSubPopName(named_size[0][1], 0)
                # resize if the type is int or float (proportion)
                if isinstance(named_size[0][0], int):
                    pop.resize(named_size[0][0])
                elif isinstance(named_size[0][0], float):
                    pop.resize(int(named_size[0][0] * pop.popSize()), propagate=True)
            elif len(size) != pop.numSubPop():
                raise ValueError('Number of subpopulations mismatch: %d in population '
                    '%d required for ExponentialGrowthModel.' % (pop.numSubPop(), len(size)))
            elif all([self._isNamedSize(x) for x in named_size]):
                # matching number of subpopulations, ...
                new_size = [x[0] for x in named_size]
                # replace None with exsiting size
                new_size = [y if x is None else x for x,y in zip(new_size, pop.subPopSizes())]
                # convert float to int
                new_size = [int(x*y) if isinstance(x, float) else x for x,y in zip(new_size, pop.subPopSizes())]
                # now resize
                pop.resize(new_size, propagate=True)
                for idx, (s,n) in enumerate(named_size):
                    if n != '':
                        pop.setSubPopName(n, idx)
            else:
                # this is a more complex resize method because it can resize and split
                # if a nested list is passed
                new_size = []
                new_names = []
                split_sizes = []
                for idx, (x, y) in enumerate(zip(named_size, pop.subPopSizes())):
                    if isinstance(x[0], int):
                        new_size.append(x[0])
                        new_names.append(x[1])
                    elif isinstance(x[0], float):
                        new_size.append(int(x[0]*y))
                        new_names.append(x[1])
                    elif x[0] is None:
                        new_size.append(y)
                        new_names.append(x[1])
                    else:  # a list
                        split_sizes.insert(0, [idx])
                        for z in x:
                            if isinstance(z[0], int):
                                split_sizes[0].append(z[0])
                            elif sys.version_info.major == 2 and isinstance(z[0], int):
                                split_sizes[0].append(z[0])
                            elif isinstance(z[0], float):
                                split_sizes[0].append(int(z[0]*y))
                            elif z[0] is None:
                                split_sizes[0].append(y)
                            else:
                                raise ValueError('Invalid size %s' % named_size)
                        new_size.append(sum(split_sizes[0][1:]))
                        new_names.append('')
                # resize and rename
                pop.resize(new_size, propagate=True)
                for idx, name in enumerate(new_names):
                    if name != '':
                        pop.setSubPopName(name, idx)
                # handle split
                indexes = [i for i, x in enumerate(named_size) if not self._isNamedSize(x)]
                indexes.reverse()
                for item in split_sizes:
                    idx = item[0]
                    new_size = item[1:]
                    names = [x[1] for x in named_size[idx]]
                    pop.splitSubPop(idx, new_size, names if any([x != '' for x in names]) else [])
        else:
            # now, if the passed population does not have any subpopulation, 
            # we can merge or split ...
            if len(named_size) == 1:
                # integer is size
                if isinstance(named_size[0][0], int):
                    pop.resize(named_size[0][0], propagate=True)
                # float is proportion
                elif isinstance(named_size[0][0], float):
                    pop.resize(int(named_size[0][0] * pop.popSize()), propagate=True)
                # None is unchanged
                if named_size[0][1] != '':
                    pop.setSubPopName(named_size[0][1], 0)
            else:
                # we need to split ...
                if not all([self._isNamedSize(x) for x in named_size]):
                    # cannot have nested population size in this case.
                    raise ValueError('Cannot fit population with size %s to size %s' %
                        (pop.subPopSizes(), named_size))
                # split subpopulations
                new_size = [x[0] for x in named_size]
                # convert None to size
                new_size = [pop.popSize() if x is None else x for x in new_size]
                # convert float to int
                new_size = [int(x*pop.popSize()) if isinstance(x, float) else x for x in new_size]
                #
                pop.resize(sum(new_size), propagate=True)
                pop.splitSubPop(0, new_size)
                for idx, (s,n) in enumerate(named_size):
                    if n != '':
                        pop.setSubPopName(n, idx)

    def _recordPopSize(self, pop):
        gen = pop.dvars().gen
        if (not hasattr(self, '_last_size')) or self._last_size != pop.subPopSizes():
            print(('%d: %s' % (gen, 
                ', '.join(
                    ['%d%s' % (x, ' (%s)' % y if y else '') for x, y in \
                        zip(pop.subPopSizes(), pop.subPopNames())])
                )))
            self._last_size = pop.subPopSizes()
        #
        if self.draw_figure:
            sz = 0
            for idx, (s, n) in enumerate(zip(pop.subPopSizes(), pop.subPopNames())):
                if n == '':
                    n = str(idx)
                if n in self.pop_base:
                    sz = max(sz, self.pop_base[n])
                self.pop_base[n] = sz
                if n in self.pop_regions:
                    self.pop_regions[n] = np.append(self.pop_regions[n],
                        np.array([[gen, sz, gen, sz+s]]))
                else:
                    self.pop_regions[n] = np.array([gen, sz, gen, sz+s], 
                        dtype=np.uint64)
                sz += s
        return True

    def plot(self, filename='', title='', initSize=[]):
        '''Evolve a haploid population using a ``RandomSelection`` mating scheme
        using the demographic model. Print population size changes duringe evolution.
        An initial population size could be specified using parameter ``initSize``
        for a demographic model with dynamic initial population size. If a filename
        is specified and if matplotlib is available, this function draws a figure
        to depict the demographic model and save it to ``filename``. An optional
        ``title`` could be specified to the figure. Note that this function can
        not be plot demographic models that works for particular mating schemes
        (e.g. genotype dependent).'''
        if not self.init_size:
            if initSize:
                self.init_size = initSize
            else:
                raise ValueError('Specific self does not have a valid initial population size')
        if filename and not has_plotter:
            raise RuntimeError('This function requires module numpy and matplotlib')
        self.draw_figure = filename and has_plotter
        self.pop_regions = OrderedDict()
        self.pop_base = OrderedDict()
        #
        self._reset()
        if title:
            print(title)
        pop = Population(self.init_size, infoFields=self.info_fields, ploidy=1)
        pop.evolve(
            matingScheme=RandomSelection(subPopSize=self),
            postOps=PyOperator(self._recordPopSize),
            gen=self.num_gens
        )
        self._recordPopSize(pop)
        # 
        if self.draw_figure:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            for name, region in list(self.pop_regions.items()):
                region = region.reshape(region.size / 4, 4)
                points = np.append(region[:, 0:2],
                    region[::-1, 2:4], axis=0)
                plt.fill(points[:,0], points[:,1], label=name, linewidth=0, edgecolor='w')
            leg = plt.legend(loc=2)
            leg.draw_frame(False)
            if title:
                plt.title(title)
            plt.savefig(filename)
            plt.close()

    def _checkSize(self, pop, param):
        gen = pop.dvars().gen
        if gen - param in self.intended_size:
            sz = self.intended_size[gen - param]
            if isinstance(sz, int):
                sz = (sz,)
            else:
                sz = tuple(sz)
            if sz != pop.subPopSizes():
                raise ValueError('Mismatch population size at generation %s (with starting gen %s): observed=%s, intended=%s' % \
                    (gen - param, param, pop.subPopSizes(), sz))
        return True
        
    def _assertSize(self, sizes, initSize=[], startGen=0):
        '''This function is provided for testing purpose.
        '''
        self.intended_size = sizes
        pop = Population(size=initSize if initSize else self.init_size,
            infoFields=self.info_fields)
        pop.dvars().gen = startGen
        pop.evolve(
            matingScheme=RandomSelection(subPopSize=self),
            postOps=PyOperator(self._checkSize, param=startGen),
            finalOps=PyOperator(self._checkSize, param=startGen),
            gen=self.num_gens
        )

    def _expIntepolate(self, N0, NT, T, x, T0=0):
        '''x=0, ... T-1
        '''
        if x == T-1:
            return NT
        elif x >= T:
            raise ValueError('Generation number %d out of bound (0<= t < %d is expected'
                % (x, T))
        else:
            return int(round(math.exp(((x+1-T0)*math.log(NT) + (T-x-1)*math.log(N0))/(T-T0))))

    def _linearIntepolate(self, N0, NT, T, x, T0=0):
        '''x=0, ... T-1
        '''
        if x == T-1:
            return NT
        elif x >= T:
            raise ValueError('Generation number %d out of bound (0<= t < %d is expected)'
                % (x, T))
        else:
            return int(round(((x+1-T0)*NT + (T-x-1)*N0)/(T-T0)))

    def _setup(self, pop):
        return True

    def _save_size(self, gen, sz):
        if self.size_cache:
            prev = [x for x in list(self.size_cache.keys()) if x < gen]
            if prev and self.size_cache[max(prev)] == sz:
                return sz
        self.size_cache[gen] = sz
        return sz

    def _cached_size(self, gen):
        if gen in self.size_cache:
            return self.size_cache[gen]
        if not self.size_cache:
            raise RuntimeError('Failed to determine size for generation {}'
                    .format(gen))
        # if we look further, keep contant size
        if gen > max(self.size_cache.keys()):
            return self.size_cache[max(self.size_cache.keys())]
        prev = [x for x in list(self.size_cache.keys()) if x < gen]
        if prev:
            return self.size_cache[max(prev)]
        else:
            raise RuntimeError('Failed to determine size for generation {}'
                    .format(gen))

    def __call__(self, pop):
        # When calling the demographic function, there are two quite separate scenarios
        #
        # 1. The demographic function is called sequentially. When a demographic model
        #   is initialized, it is considered part of its own generation zero. There can
        #   be two cases, the first one for initialization, the other one for instant
        #   population change (in the case of InstantChangeModel).
        #
        # 2. When the demographic function is called randomly, only a previously evolved
        #   generation is allowed. This is because many demographic models depend on
        #   size of previous generations to determine new population size, which makes
        #   jumping forward impossible.
        # 
        # In the case of multi-model, random access will not call fitSize for any
        # intermediate steps.
        # 
        # the first time this demographic function is called
        self._use_cached = False
        if not hasattr(self, '_start_gen'):
            self._reset()
            self._start_gen = pop.dvars().gen
            self._last_gen = pop.dvars().gen
            # resize populations if necessary
            self._fitToSize(pop, self._raw_init_size)
            # by this time, we should know the initial population size
            self.init_size = pop.subPopSizes()
            # then we can set up model if the model depends on initial
            # population size
            self._setup(pop)
        elif pop.dvars().gen != self._last_gen + 1:
            self._use_cached = True
            return self._cached_size(pop.dvars().gen)
        #
        self._gen = pop.dvars().gen - self._start_gen
        self._last_gen = pop.dvars().gen
        #
        pop.dvars()._gen = self._gen
        pop.dvars()._num_gens = self.num_gens
        for op in self.ops:
            if not op.apply(pop):
                self._reset()
                return []
        if '_expected_size' in pop.vars():
            return self._save_size(pop.dvars().gen, pop.vars().pop('_expected_size'))
        else:
            return self._save_size(pop.dvars().gen, pop.subPopSizes())

class ExponentialGrowthModel(DemographicModel):
    '''A model for exponential population growth with carry capacity'''
    def __init__(self, T=None, N0=[], NT=None, r=None, ops=[], infoFields=[]):
        '''An exponential population growth model that evolves a population from size
        ``N0`` to ``NT`` for ``T`` generations with ``r*N(t)`` individuals added
        at each generation. ``N0``, ``NT`` and ``r`` can be a list of population
        sizes or growth rates for multiple subpopulations. The initial population
        will be resized to ``N0`` (split if necessary). Zero or negative growth
        rates are allowed. The model will automatically determine ``T``, ``r``
        or ``NT`` if one of them is unspecified. If all of them are specified,
        ``NT`` is intepretted as carrying capacity of the model, namely the 
        population will keep contant after it reaches size ``NT``. Optionally,
        one or more operators (e.g. a migrator) ``ops`` can be applied to 
        population. '''
        DemographicModel.__init__(self, T, N0, ops, infoFields)
        #
        if [x is None or x == [] for x in [T, NT, r]].count(True) > 1:
            raise ValueError('Please specify at least two parameters of T, NT and r')
        self.T = T
        self.N0 = N0
        self.NT = NT
        self.r = r

    def _setup(self, pop):
        if self.r is None:
            # self.T must be known, which is num_gens
            # self.NT must also be known
            self.NT = self._extractSize(self.NT)
            if self.NT is None or None in self.NT or \
                any([isinstance(x, float) for x in self.NT]):
                raise ValueError('Relative ending population size is not allowed '
                    'for LinearGrowthModel')
        elif isinstance(self.r, (int, float)):
            # if we do not know generation number, we are in some trouble
            if self.num_gens < 0:
                # get number of generations from NT and r
                self.NT = self._extractSize(self.NT)
                if self.NT is None or None in self.NT or \
                    any([isinstance(x, float) for x in self.NT]):
                    raise ValueError('Relative ending population size is not allowed '
                        'for LinearGrowthModel')
                if len(self.NT) != len(self.init_size):
                    raise ValueError('Starting and ending population should have the '
                        'same number of subpopulations')
                # what is the number of generations
                if self.r == 0:
                    raise ValueError('Cannot reach destination size with r=0')
                T = [int(round((math.log(y) - math.log(x)) / self.r)) + 1 for (x,y) in zip(self.init_size, self.NT)]
                if max(T) < 0:
                    raise ValueError('Cannot reach destination size in this configuraiton.')
                self.num_gens = max(T + [1])
                #
            elif self.NT is None:
                # get final population size from T and r
                self.NT = [x*math.exp(self.r*self.num_gens) for x in self.init_size]
            elif None in self.NT or \
                any([isinstance(x, float) for x in self.NT]):
                raise ValueError('Relative ending population size is not allowed'
                    'for ExponentialGrowthModel')
            self.r = [self.r] * len(self.NT)
        elif isinstance(self.r, (list, tuple)):
            if len(self.r) != len(self.init_size):
                raise ValueError('Please specify growth rate for each subpopulation '
                    'if multiple growth rates are specified.')
            if self.NT is None:
                self.NT = [x*math.exp(y*self.num_gens) for x,y in zip(self.init_size, self.r)]
            elif None in self.NT or \
                any([isinstance(x, float) for x in self.NT]):
                raise ValueError('Relative ending population size is not allowed'
                    'for ExponentialGrowthModel')
        else:
            raise ValueError('Unacceptable growth rate (a number or a list of numbers '
                'is expected')

    def __call__(self, pop):
        if not DemographicModel.__call__(self, pop):
           return []
        #
        if self._use_cached:
            return self._cached_size(pop.dvars().gen)
        #
        if self._gen == self.num_gens:
            return []
        elif self.r is None:
            return self._save_size(pop.dvars().gen, 
                [self._expIntepolate(n0, nt, self.num_gens, self._gen)
                for (n0, nt) in zip(self.init_size, self.NT)])
        else:
            # with r ...
            return self._save_size(pop.dvars().gen, [min(int(round(nt)), int(round(n0 * math.exp(r * (self._gen + 1)))))
                for (n0, nt, r) in zip(self.init_size, self.NT, self.r)])


class LinearGrowthModel(DemographicModel):
    '''A model for linear population growth with carry capacity.'''
    def __init__(self, T=None, N0=[], NT=None, r=None, ops=[], infoFields=[]):
        '''An linear population growth model that evolves a population from size
        ``N0`` to ``NT`` for ``T`` generations with ``r*N0`` individuals added
        at each generation. ``N0``, ``NT`` and ``r`` can be a list of population
        sizes or growth rates for multiple subpopulations. The initial population
        will be resized to ``N0`` (split if necessary). Zero or negative growth
        rates are allowed. The model will automatically determine ``T``, ``r``
        or ``NT`` if one of them is unspecified. If all of them are specified,
        ``NT`` is intepretted as carrying capacity of the model, namely the 
        population will keep contant after it reaches size ``NT``. Optionally,
        one or more operators (e.g. a migrator) ``ops`` can be applied to 
        population. '''
        DemographicModel.__init__(self, T, N0, ops, infoFields)
        #
        if [x is None or x == [] for x in [T, NT, r]].count(True) > 1:
            raise ValueError('Please specify at least two parameters of T, NT and r')
        self.T = T
        self.N0 = N0
        self.NT = NT
        self.r = r

    def _setup(self, pop):
        if self.r is None:
            # self.T must be known, which is num_gens
            # self.NT must also be known
            self.NT = self._extractSize(self.NT)
            if self.NT is None or None in self.NT or \
                any([isinstance(x, float) for x in self.NT]):
                raise ValueError('Relative ending population size is not allowed '
                    'for LinearGrowthModel')
        elif isinstance(self.r, (int, float)):
            # if we do not know generation number, we are in some trouble
            if self.num_gens < 0:
                # get number of generations from NT and r
                self.NT = self._extractSize(self.NT)
                if self.NT is None or None in self.NT or \
                    any([isinstance(x, float) for x in self.NT]):
                    raise ValueError('Relative ending population size is not allowed '
                        'for LinearGrowthModel')
                if len(self.NT) != len(self.init_size):
                    raise ValueError('Starting and ending population should have the '
                        'same number of subpopulations')
                # what is the number of generations
                T = [int((y - x) / (x * self.r)) for (x,y) in zip(self.init_size, self.NT)]
                self.num_gens = max(T + [1])
                #
            elif self.NT is None:
                # get final population size from T and r
                self.NT = [int(x*(1+self.r*self.num_gens)) for x in self.init_size]
            elif None in self.NT or \
                any([isinstance(x, float) for x in self.NT]):
                raise ValueError('Relative ending population size is not allowed'
                    'for LinearGrowthModel')
            self.r = [self.r] * len(self.NT)
        elif isinstance(self.r, (list, tuple)):
            if len(self.r) != len(self.init_size):
                raise ValueError('Please specify growth rate for each subpopulation '
                    'if multiple growth rates are specified.')
            if self.NT is None:
                self.NT = [int(x*(1+y*self.num_gens)) for x,y in zip(self.init_size, self.r)]
            elif None in self.NT or \
                any([isinstance(x, float) for x in self.NT]):
                raise ValueError('Relative ending population size is not allowed'
                    'for LinearGrowthModel')
        else:
            raise ValueError('Unacceptable growth rate (a number or a list of numbers '
                'is expected')

    def __call__(self, pop):
        if not DemographicModel.__call__(self, pop):
            return []
        #
        if self._use_cached:
            return self._cached_size(pop.dvars().gen)
        #
        if self._gen == self.num_gens:
            return []
        elif self.r is None:
            # no r use intepolation
            return self._save_size(pop.dvars().gen, [self._linearIntepolate(n0, nt, self.num_gens, self._gen)
                for (n0, nt) in zip(self.init_size, self.NT)])
        else:
            # with r ...
            return self._save_size(pop.dvars().gen, [min(nt, int(n0 + n0 * (self._gen + 1.) * r))
                for (n0, nt, r) in zip(self.init_size, self.NT, self.r)])


class InstantChangeModel(DemographicModel):
    '''A model for instant population change (growth, resize, merge, split).'''
    def __init__(self, T=None, N0=[], G=[], NG=[], ops=[], infoFields=[], removeEmptySubPops=False):
        '''An instant population growth model that evolves a population
        from size ``N0`` to ``NT`` for ``T`` generations with population
        size changes at generation ``G`` to ``NT``. If ``G`` is a list,
        multiple population size changes are allowed. In that case, a list
        (or a nested list) of population size should be provided to parameter
        ``NT``. Both ``N0`` and ``NT`` supports fixed (an integer), dynamic
        (keep passed poulation size) and proportional (an float number) population
        size. Optionally, one or more operators (e.g. a migrator) ``ops``
        can be applied to population. Required information fields by these
        operators should be passed to parameter ``infoFields``. If ``removeEmpty``
        option is set to ``True``, empty subpopulation will be removed. This
        option can be used to remove subpopulations.'''
        DemographicModel.__init__(self, T, N0, ops, infoFields)
        if isinstance(G, int):
            self.G = [G]
            if isinstance(NG, int) :
                self.NG = [NG]
            elif isinstance(NG[0], int):
                self.NG = [NG]
            else:
                self.NG = NG
        else:
            if not isinstance(NG, (tuple, list)):
                raise ValueError('Multiple sizes should be specified if multiple G is provided.')
            if len(G) != len(NG):
                raise ValueError('Please provide population size for each growth generation.')
            self.G = sorted(G)
            self.NG = [NG[G.index(x)] for x in self.G]
        #
        for g in self.G:
            if g < 0 or g >= self.num_gens:
                raise ValueError('Population change generation %d exceeds '
                    'total number of generations %d' % (g, self.num_gens))
        self.removeEmptySubPops = removeEmptySubPops

    def __call__(self, pop):
        # this one fixes N0... (self.init_size)
        if not DemographicModel.__call__(self, pop):
            return []
        #
        if self._use_cached:
            return self._cached_size(pop.dvars().gen)
        # 
        if self._gen in self.G:
            sz = self.NG[self.G.index(self._gen)]
            self._fitToSize(pop, sz)
        if self.removeEmptySubPops:
            pop.removeSubPops([idx for idx,x in enumerate(pop.subPopSizes()) if x==0])
        sz = pop.subPopSizes()
        return self._save_size(pop.dvars().gen, sz)



class AdmixtureModel(DemographicModel):
    '''An population admixture model that mimicks the admixing of two 
    subpopulations using either an HI (hybrid isolation) or CGF
    (continuous gene flow) model. In the HI model, an admixed population
    is created instantly from two or more subpopulation with
    specified proportions. In the CGF model, a population will be reduced
    to (1-alpha) N and accept alpha*N individuals from another population
    for a number of generations. Please see Long (1990) The Stnetic Structure
    of Admixed Populations, Genetics for details. 
    
    This model is deprecated due to the introduction of event based
    implementation of admixture models (e.g. ``AdmixToNewPopAdmixture``
    and ``ContinuousGeneFlowAdmixture events``). 
    '''
    def __init__(self, T=None, N0=[], model=[], ops=[], infoFields=[]):
        '''Define a population admixture model that mixed two 
        subpopulations. The population has an initial size of ``N0``
        which contains more than one subpopulations. The demographic
        model will evolve for ``T`` generations with admixture model
        ``model``, which should be ``['HI', parent1, parent2, mu, name]`` or
        ``['CGF', receipient, doner, alpha]``. In the first case,
        a new admixed population is created with a proportion of mu, and
        1-mu individuals from parent1 and parent2 population. An optional
        ``name`` can be assigned to be new subpopulation. In the latter case,
        1-mu percent of individuals in receipient population will be replaced by
        individuals in the doner population.
        '''
        DemographicModel.__init__(self, T, N0, ops, infoFields)
        if not model or model[0] not in ('HI', 'CGF') or len(model) < 4 or model[1] == model[2] or model[3] > 1:
            raise ValueError('model should be either ("HI", parent1, paren2, mu)'
                'or ("CGF", receipient, doner, alpha)')
        self.model = model

    def HI_size(self, N1, N2, mu):
        #
        #   N1         N2
        # choose  
        #   N*mu       N*(1-mu)
        # so that
        #   N*mu < N1
        #   N*(1-mu) < N2
        # That is to say
        #   N < N1/mu
        #   N < N2/(1-mu)
        #
        if mu == 0:
            return (0, N2)
        elif mu == 1:
            return (N1, 0)
        else:
            N = min(N1/mu, N2/(1-mu))
            return (int(N*mu + 0.5), int(N*(1-mu) + 0.5))

    def __call__(self, pop):
        if not DemographicModel.__call__(self, pop):
            return []
        #
        if self._use_cached:
            return self._cached_size(pop.dvars().gen)
        #
        if self.model[0] == 'HI':
            if self._gen != 0:
                return pop.subPopSizes()
            if self.model[1] >= pop.numSubPop() or self.model[2] >= pop.numSubPop():
                raise RuntimeError('Failed to mixed populations {} and {}'.format(self.model[1], self.model[2]))
            sz1 = pop.subPopSize(self.model[1])
            sz2 = pop.subPopSize(self.model[2])
            #
            # mu from sz1, 1-mu from p2
            sz = self.HI_size(sz1, sz2, self.model[3])
            # create a new subpopulation
            pop1 = pop.clone()
            admixed_size = [0] * pop1.numSubPop()
            admixed_size[self.model[1]] = sz[0]
            admixed_size[self.model[2]] = sz[1]
            pop1.resize(admixed_size)
            if len(self.model) > 4:
                # assign a new name
                pop1.mergeSubPops(ALL_AVAIL, self.model[4])
            else:
                pop1.mergeSubPops(ALL_AVAIL)
            pop.addIndFrom(pop1)
        else:
            if self.model[1] >= pop.numSubPop() or self.model[2] >= pop.numSubPop():
                raise RuntimeError('Failed to mixed populations {} and {}'.format(self.model[1], self.model[2]))
            sz1 = pop.subPopSize(self.model[1])
            sz2 = pop.subPopSize(self.model[2])
            # requested number of individuals from sz2
            request = min(int(sz1 * (1-self.model[3]) + 0.5), sz2)
            # we need to replace requested number of individuals in sz1 
            # with individuals in sz2
            #
            # step1, enlarge doner
            sz = list(pop.subPopSizes())
            sz[self.model[1]] -= request
            sz[self.model[2]] += request
            pop.resize(sz, propagate=True)
            # step2, split doner
            pop.splitSubPop(self.model[2], [request, pop.subPopSize(self.model[2]) - request])
            # step3, merge the split indiiduals to subpopulation self.model[1]
            pop.mergeSubPops([self.model[1], self.model[2]])
        return self._save_size(pop.dvars().gen, pop.subPopSizes())


class MultiStageModel(DemographicModel):
    '''A multi-stage demographic model that connects a number of demographic
    models. '''
    def __init__(self, models, ops=[], infoFields=[]):
        '''An multi-stage demographic model that connects specified
        demographic models ``models``. It applies a model to the population
        until it reaches ``num_gens`` of the model, or if the model returns
        ``[]``. One or more operators could be specified, which will be applied
        before a demographic model is applied. Note that the last model will be
        ignored if it lasts 0 generation.'''
        flds = []
        gens = []
        for x in models:
            flds.extend(x.info_fields)
            gens.append(x.num_gens)
        if all([x>=0 for x in gens]):
            total_gens = sum(gens)
        else:
            total_gens = -1
        DemographicModel.__init__(self, numGens=total_gens,
            initSize=models[0].init_size, ops=ops, infoFields=flds+infoFields)
        #
        self.models = models
        if self.models[-1].num_gens == 0:
            raise ValueError('The last demographic model in a MultiStageModel cannot last zero generation.')
        self._model_idx = 0
        self._model_start_gen = 0

    def _reset(self):
        self._model_idx = 0
        self._model_start_gen = 0
        if hasattr(self, '_start_gen'):
            del self._start_gen
        for m in self.models:
            if hasattr(m, '_start_gen'):
                del m._start_gen
  
    def _advance(self, pop):
        self._model_idx += 1
        self._model_start_gen = pop.dvars().gen - self._start_gen
        while True:
            if self._model_idx == len(self.models):
                self._reset()
                return []
            # call and skip
            if self.models[self._model_idx].num_gens == 0:
                sz = self.models[self._model_idx].__call__(pop)
                self._model_idx += 1
                continue
            sz = self.models[self._model_idx].__call__(pop)
            if sz:
                return sz
            else:
                self._model_idx += 1
                continue

    def __call__(self, pop):
        # in a special case when the demographic model has been
        # initialized but the population has revert to a previous stage
        # so that the starting gen is after the current gen. We will
        # need to handle this case.
        # determines generation number internally as self.gen
        if not DemographicModel.__call__(self, pop):
            return []
        #
        if self._use_cached:
            return self._cached_size(pop.dvars().gen)
        #
        # There are three cases
        # 1. within the current model, a valid step is returned
        #   --> return
        # 2. within the current model, a [] is returned by a 
        #   terminator.
        #   --> proceed to the next available model, call, and return
        # 3. at the end of the current model,
        #   --> proceed to the next available model, call, and return
        # 4. at the beginning of a zero-step model
        #   --> call
        #   --> process to the next available model, call, and return
        #
        # in the middle
        if self.models[self._model_idx].num_gens < 0 or \
            self.models[self._model_idx].num_gens > self._gen - self._model_start_gen:
            sz = self.models[self._model_idx].__call__(pop)
            if not sz:
                sz = self._advance(pop)
        elif self.models[self._model_idx].num_gens == 0:
            sz = self.models[self._model_idx].__call__(pop)
            sz = self._advance(pop)
        elif self.models[self._model_idx].num_gens <= self._gen - self._model_start_gen:
            sz = self._advance(pop)
        return self._save_size(pop.dvars().gen, sz)


        
class EventBasedModel(DemographicModel):
    '''An event based demographic model in which the demographic changes are 
    triggered by demographic events such as population growth, split, join, and 
    admixture. The population size will be kept constant if no event is applied
    at a certain generation.
    '''
    def __init__(self, events=[], T=None, N0=[], ops=[], infoFields=[]):
        '''A demographic model that is driven by a list of demographic events.
        The events should be subclasses of ``DemographicEvent``, which have similar
        interface as regular operators with the exception that applicable parameters
        ``begin``, ``end``, ``step``, ``at`` are relative to the demographic model,
        not the population.
        '''
        if isinstance(events, DemographicEvent):
            events = [events]
        DemographicModel.__init__(self, numGens=T, initSize=N0,
            ops=ops + events, infoFields=infoFields)


class DemographicEvent:
    '''A demographic events that will be applied to one or more populations at
    specified generations. The interface of a DemographicEvent is very similar to
    an simuPOP operator, but the applicable parameters are handled so that
    the generations are relative to the demographic model, not the populations
    to which the event is applied.
    '''
    def __init__(self, ops=[], output='', begin=0, end=-1, step=1, at=[],
        reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[]):
        '''Create a demographic event that will be applied at specified
        generations according to applicability parameters ``reps``, ``begin``, 
        ``end``, ``step`` and ``at``. Parameter ``subPops`` is usually used
        to specify the subpopulations affected by the event. One or more simuPOP
        operators, if specified in ``ops``, will be applied when the event
        happens. Parameters ``output`` and ``infoFields`` are currently ignored. 
        '''
        if isinstance(ops, (list, tuple)):
            self.ops = ops
        else:
            self.ops = [ops]
        self.output = output
        self.begin = begin
        self.end = end
        self.step = step
        if isinstance(at, int):
            self.at = [at]
        else:
            self.at = at
        self.reps = reps
        self.subPops = subPops
        self.infoFields = infoFields
    
    def _applicable(self, pop):
        if '_gen' not in pop.vars() or '_num_gens' not in pop.vars():
            raise ValueError('Cannot apply to a population without variables _gen or _num_gens')
        #
        gen = pop.dvars()._gen
        end = pop.dvars()._num_gens - 1
        rep = pop.dvars().rep
        #
        if self.reps != ALL_AVAIL and (\
            (isinstance(self.reps, int) and self.reps != rep) or \
            rep not in self.reps):
            return False
        #
        if self.at:
            for a in self.at:
                if a >= 0:
                    if a == gen:
                        return True
                    else:
                        continue
                else:
                    if end < 0:
                        continue
                        #raise ValueError('Cannot specify negative at generation for a demographic '
                        #    'model without fixed number of generations.')
                    if end + a + 1 == gen:
                        return True
                    else:
                        continue
            return False
        #
        if end < 0:
            if self.begin < 0 or self.begin > gen:
                return False
            if ((gen - self.begin) % self.step == 0) and (end < 0 or end >= gen):
                return True
            else:
                return False
        else:
            rstart = self.begin if self.begin >= 0 else self.begin + end + 1
            rend = self.end if self.end >= 0 else self.end + end + 1
            if rstart > rend:
                return False
            return gen >= rstart and gen <= rend and (gen - rstart) % self.step == 0
        return False

    def _identifySubPops(self, pop):
        if self.subPops is ALL_AVAIL:
            return list(range(pop.numSubPop()))
        #
        ret = []
        names = pop.subPopNames()
        if isinstance(self.subPops, (str, int)):
            sps = [self.subPops]
        else:
            sps = self.subPops
        for sp in sps:
            if type(sp) == int:
                ret.append(sp)
            else:
                if sp not in names:
                    raise ValueError('Invalid subpopulation name {}'.format(sp))
                ret.append(names.index(sp))
        return ret

    def apply(self, pop):
        if self.ops and self._applicable(pop):
            for op in self.ops:
                if not op.apply(pop):
                    return False
        return True



class ResizeEvent(DemographicEvent):
    '''A demographic event that resize specified subpopulations'''
    def __init__(self, sizes=[], names=[], removeEmptySubPops=False,
        ops=[], output='', begin=0, end=-1, 
        step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[]):
        '''A demographic event that resizes given subpopulations ``subPops`` to new
        ``sizes`` (integer type), or sizes proportional to original sizes (if a float
        number is given). For example, ``sizes=[0.5, 500]`` will resize the first
        subpopulation to half of its original size, and the second subpopulation to
        size ``500``. If the new size is larger, existing individuals will be copied
        to sequentially, and repeatedly if needed. If the size of a subpopulation is
        0 and ``removeEmptySubPops`` is ``True``, empty subpopulations will be
        removed. A new set of names could be assigned to the population being resized. '''
        self.sizes = sizes
        if isinstance(names, str):
            self.names = [names]
        else:
            self.names = names
        self.removeEmptySubPops = removeEmptySubPops
        DemographicEvent.__init__(self, ops, output, begin, end, step, at, reps,
            subPops, infoFields)

    def apply(self, pop):
        if not DemographicEvent.apply(self, pop):
            return False
        #
        if not self._applicable(pop):
            return True
        # identify applicable subpopulations
        subPops = self._identifySubPops(pop)
        if isinstance(self.sizes, (tuple, list)) and len(subPops) != len(self.sizes):
            raise ValueError('If a list of sizes are specified, please specify them for all subpopulations.')
        if self.names and len(self.names) != len(subPops):
            raise ValueError('Legnth of names (%d) and subPops mismatch (%d)' % (len(self.names), len(subPops)))
        #
        if '_expected_size' in pop.vars():
            # if there was an intended size, realized it.
            pop.resize(pop.vars().pop('_expected_size'))
        # new size ...
        sz = list(pop.subPopSizes())
        for idx, sp in enumerate(subPops):
            if isinstance(self.sizes, int):
                sz[sp] = self.sizes
            elif isinstance(self.sizes, float):
                sz[sp] = int(round(sz[sp] * self.sizes))
            elif isinstance(self.sizes, (list, tuple)) and isinstance(self.sizes[idx], int):
                sz[sp] = self.sizes[idx]
            elif isinstance(self.sizes, (list, tuple)) and isinstance(self.sizes[idx], float):
                sz[sp] = int(round(self.sizes[idx] * sz[sp]))
            else:
                raise ValueError('Unrecognized sizes parameter.')
        #
        pop.resize(sz, propagate=True)
        if self.removeEmptySubPops:
            pop.removeSubPops([idx for idx,x in enumerate(pop.subPopSizes()) if x==0])
        if self.names:
            for sp,name in zip(subPops, self.names):
                pop.setSubPopName(name, sp)
        return True


class MergeEvent(DemographicEvent):
    '''A demographic event that merges one or more subpopulation to
    a single one.'''
    def __init__(self, name='', ops=[], output='', begin=0, end=-1, 
        step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[]):
        '''A demographic event that merges subpopulations into a single subpopulation.
        The merged subpopulation will have the name of the first merged subpopulation
        unless a separate ``name`` is supported.'''
        self.name = name
        DemographicEvent.__init__(self, ops, output, begin, end, step, at, reps,
            subPops, infoFields)

    def apply(self, pop):
        if not DemographicEvent.apply(self, pop):
            return False
        #
        if not self._applicable(pop):
            return True
        #
        if '_expected_size' in pop.vars():
            # if there was an intended size, realized it.
            pop.resize(pop.vars().pop('_expected_size'))
        # identify applicable subpopulations
        subPops = self._identifySubPops(pop)
        #
        pop.mergeSubPops(subPops=subPops, name=self.name)
        return True


class SplitEvent(DemographicEvent):
    '''A demographic event that splits a specified population into two or
    more subpopulations.'''
    def __init__(self, sizes=[], names=[], ops=[], output='', begin=0, end=-1, 
        step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[]):
        '''A demographic event that splits a subpopulation specified by
        ``subPops`` to two or more subpopulations, with specified ``sizes``
        and ``names``. ``sizes`` can be a list of numbers, proportions
        (e.g. ``[1., 500]`` keeps the original population and copies 500
        individuals to create a new subpupulation). Note that ``sizes``
        and ``names``, if specified, should include the source subpopulation
        as its first element.
        '''
        self.sizes = sizes
        self.names = names
        DemographicEvent.__init__(self, ops, output, begin, end, step, at, reps,
            subPops, infoFields)
    
    def apply(self, pop):
        if not DemographicEvent.apply(self, pop):
            return False
        #
        if not self._applicable(pop):
            return True
        # identify applicable subpopulations
        subPops = self._identifySubPops(pop)
        if len(subPops) != 1:
            raise ValueError('Please specify one and only one subpopulation for event CopyEvent.')
        subPop = subPops[0]
        #
        if '_expected_size' in pop.vars():
            # if there was an intended size, realized it.
            pop.resize(pop.vars().pop('_expected_size'))
        # new size ...
        sz = []
        for x in self.sizes:
            if isinstance(x, int):
                sz.append(x)
            elif isinstance(x, float):
                sz.append(int(round(x * pop.subPopSize(subPop))))
            else:
                raise ValueError('Elements of parameter sizes can only be interger (number'
                    'of individuals) or float numbers (proportions relative to the size of '
                    'subpopulation being split)')
        #
        if sum(sz) != pop.subPopSize(subPop):
            new_sz = list(pop.subPopSizes())
            new_sz[subPop] = sum(sz)
            pop.resize(new_sz, propagate=True)
        pop.splitSubPop(subPop, sz, self.names)
        return True

class ExpansionEvent(DemographicEvent):
    '''A demographic event that increase applicable population size by
    ``N*r`` (to size ``N*(1+r)``), or ``s`` (to size ``N+s``) at each applicable
    generation. The first model leads to an exponential population expansion
    model with rate ``r`` (``N(t)=N(0)*exp(r*t)``), where the second model leads to
    an linear population growth model (``N(t)=N(0)+s*t``) and this is why the
    parameter is called ``slopes``. Note that if both population
    size and ``r`` are small (e.g. ``N*r<1``),  the population might not expand
    as expected.'''
    def __init__(self, rates=[], slopes=[], capacity=[], name='', 
        ops=[], output='', begin=0, end=-1, 
        step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[]):
        '''A demographic event that expands all or specified subpopulations
        (``subPops``) exponentially by a rate of ``rates``, or linearly by a slope
        of ``slopes``, unless carray capacity (``capacity``) of the population has
        been reached. Parameter ``rates`` can be a single number or a list of rates
        for all subpopulations. Parameter ``slopes`` should be a number, or a list
        of numbers for all subpopulations. ``subPops`` can be a ``ALL_AVAIL`` or a list
        of subpopulation index or names. ``capacity`` can be empty (no limit on
        carrying capacity), or one or more numbers for each of the subpopulations.
        '''
        self._N0 = None
        self._T0 = None
        self.rates = rates
        self.slopes = slopes
        if self.rates and self.slopes:
            raise ValueError('Please specify one and only one of parameter rates and slopes.')
        self.capacity = capacity
        DemographicEvent.__init__(self, ops, output, begin, end, step, at, reps,
            subPops, infoFields)
    
    def apply(self, pop):
        if not self._applicable(pop):
            return True
        #
        if not DemographicEvent.apply(self, pop):
            return False
        #
        # identify applicable subpopulations
        subPops = self._identifySubPops(pop)
        if self.slopes and isinstance(self.slopes, (list, tuple)) and len(self.slopes) != len(subPops):
            raise ValueError('Please specify slopes for all subpopulations or '
                'each of the {} subpopulations'.format(len(subPops)))
        if self.rates and isinstance(self.rates, (list, tuple)) and len(self.rates) != len(subPops):
            raise ValueError('Please specify growth rate for all subpopulations or '
                'each of the {} subpopulations'.format(len(subPops)))
        if self.capacity and isinstance(self.capacity, (list, tuple)) and len(self.capacity) != len(subPops):
            raise ValueError('If specified, please specify carrying capacity for all '
                'subpopulations or each of the {} subpopulations'.format(len(subPops)))
        #
        if '_expected_size' in pop.vars():
            sz = pop.dvars()._expected_size
        else:
            sz = list(pop.subPopSizes())
        #
        if self._N0 is None or len(self._N0) != len(subPops):
            self._N0 = [sz[x] for x in subPops]
            self._T0 = pop.dvars().gen
        for idx, sp in enumerate(subPops):
            if self.rates:
                if type(self.rates) in [list, tuple]:
                    # if current size match, use it to do next
                    if sz[sp] == int(round(self._N0[idx] * math.exp((pop.dvars().gen - self._T0) * self.rates[idx]))):
                        sz[sp] = int(round(self._N0[idx] * math.exp((pop.dvars().gen - self._T0 + 1) * self.rates[idx])))
                    # otherwise reset ...
                    else:
                        self._N0 = [sz[x] for x in subPops]
                        self._T0 = pop.dvars().gen
                        sz[sp] = int(round(sz[sp] * math.exp(self.rates[idx])))
                else:
                    if sz[sp] == int(round(self._N0[idx] * math.exp((pop.dvars().gen - self._T0) * self.rates))):
                        sz[sp] = int(round(self._N0[idx] * math.exp((pop.dvars().gen - self._T0 + 1) * self.rates)))
                    # otherwise reset ...
                    else:
                        self._N0 = [sz[x] for x in subPops]
                        self._T0 = pop.dvars().gen
                        sz[sp] = int(round(sz[sp] * math.exp(self.rates)))
            else:
                if type(self.slopes) in [list, tuple]:
                    # if current size match, use it to do next
                    if sz[sp] == int(round(self._N0[idx] + (pop.dvars().gen - self._T0) * self.slopes[idx])):
                        sz[sp] = int(round(self._N0[idx] + (pop.dvars().gen - self._T0 + 1) * self.slopes[idx]))
                    # otherwise reset ...
                    else:
                        self._N0 = [sz[x] for x in subPops]
                        self._T0 = pop.dvars().gen
                        sz[sp] = int(round(sz[sp] + self.slopes[idx]))
                else:
                    # if current size match, use it to do next
                    if sz[sp] == int(round(self._N0[idx] + (pop.dvars().gen - self._T0) * self.slopes)):
                        sz[sp] = int(round(self._N0[idx] + (pop.dvars().gen - self._T0 + 1) * self.slopes))
                    # otherwise reset ...
                    else:
                        self._N0 = [sz[x] for x in subPops]
                        self._T0 = pop.dvars().gen
                        sz[sp] = int(round(sz[sp] + self.slopes))
        if self.capacity:
            for idx, sp in enumerate(subPops):
                if isinstance(self.capacity, (list, tuple)):
                     sz[sp] = min(sz[sp], self.capacity[idx])
                else:
                     sz[sp] = min(sz[sp], self.capacity)
        pop.dvars()._expected_size = sz 
        return True


class AdmixtureEvent(DemographicEvent):
    '''This event implements a population admixture event that mix
    individuals from specified subpopulations to either a new 
    subpopulation or an existing subpopulation.'''
    def __init__(self, sizes=[], toSubPop=None, name='',
        ops=[], output='', begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, 
        subPops=ALL_AVAIL, infoFields=[]):
        '''Create an admixed population by choosing individuals
        from all or specified subpopulations (``subPops``) and creating
        an admixed population ``toSubPop``. The admixed population will
        be appended to the population as a new subpopulation with name
        ``name`` if ``toSubPop`` is ``None`` (default), or replace an
        existing subpopulation with name or index ``toSubPop``. The admixed
        population consists of individuals from ``subPops`` according to
        specified ``sizes``. Its size is maximized to have the largest
        number of individuals from the source population when a new population
        is created, or equal to the size of the existing destination population.
        The parameter ``sizes`` should be a list of float numbers 
        between 0 and 1, and add up to 1 (e.g. ``[0.4, 0.4, 0.2]``, although
        this function ignores the last element and set it to 1 minus the 
        sum of the other numbers). Alternatively, parameter ``sizes`` can
        be a list of numbers used to explicitly specify the size of admixed
        population and number of individuals from each source subpopulation.
        In all cases, the size of source populations will be kept constant.
        '''
        DemographicEvent.__init__(self, ops, output, begin, end, step, at, reps,
            subPops, infoFields)
        if all([isinstance(x, int) for x in sizes]):
            self.numbers = sizes
            self.proportions = None
        else:
            self.numbers = None
            self.proportions = sizes
        self.subPopName = name
        self.toSubPop = toSubPop

    def apply(self, pop):
        if not self._applicable(pop):
            return True
        #
        if not DemographicEvent.apply(self, pop):
            return False
        #
        # identify applicable subpopulations
        subPops = self._identifySubPops(pop)
        #
        if '_expected_size' in pop.vars():
            # if there was an intended size, realized it.
            pop.resize(pop.vars().pop('_expected_size'))
        #
        if self.toSubPop is None:
            toSubPop = self.toSubPop
        else:
            # replacing an existing subpopulation
            if isinstance(self.toSubPop, int):
                if self.toSubPop >= pop.numSubPop():
                    raise ValueError('Subpopulation index {} out of range'.format(self.toSubPop))
                toSubPop = self.toSubPop
            else:
                if self.toSubPop not in pop.subPopNames():
                    raise ValueError('No subpopulation with name {} can be located'.format(self.toSubPop))
                toSubPop = pop.subPopNames().index(self.toSubPop)
        #
        if self.proportions and len(subPops) != len(self.proportions):
            raise ValueError('Number of subpopulations and proportions mismatch.')
        if self.numbers and len(subPops) != len(self.numbers):
            raise ValueError('Number of subpopulations and proportions mismatch.')
        #
        # determine the maximum number of individuals that 
        # can be draw from each subpopulation
        if self.numbers:
            # if specific numbers are specified
            num_migrants = []
            for sp in range(pop.numSubPop()):
                if sp in subPops:
                    num_migrants.append(self.numbers[subPops.index(sp)])
                else:
                    # not involved
                    num_migrants.append(0)
            # if to a certain subpopulation ...
            if toSubPop is not None:
                # 100, 200, 300
                #
                # pick 200, 200 to subpop 2
                #
                # we are supposed to keep num_migrants[toSubPop] individuals so the
                # number of migrant is - (N - migrants)
                #
                num_migrants[toSubPop] -= pop.subPopSize(toSubPop)
        else:
            if sum(self.proportions[:-1]) > 1:
                raise ValueError('Proportions of individual from parental populations add up more than 1.')
            if any([x < 0 or x > 1 for x in self.proportions]):
                raise ValueError('Proportion of one of the parental populations is negative or more than 1.')
            if sum(self.proportions) != 1.:
                self.proportions[-1] = 1. - sum(self.proportions[:-1])
            # 
            # create a new subpopulation
            if toSubPop is None:
                # now determine the size ... try different source subpopulation
                num_migrants = None
                for sp in range(pop.numSubPop()):
                    if sp not in subPops:
                        continue
                    idx = subPops.index(sp)
                    if self.proportions[idx] == 0:
                        continue
                    # if we use all individuals in this subpopulation
                    N = pop.subPopSize(sp) / self.proportions[idx]
                    num_migrants = [int(N*self.proportions[subPops.index(x)]) if x in subPops else 0 for x in range(pop.numSubPop())]
                    if all([x <= y for x,y in zip(num_migrants, pop.subPopSizes())]):
                        break
                if num_migrants is None:
                    raise RuntimeError('Failed to determine size of admixed subpopulation.')
            else:
                N = pop.subPopSize(toSubPop)
                # proportion of individules from each source subpopulation (migrate out)
                num_migrants = [int(N*self.proportions[subPops.index(x)]) if x in subPops else 0 for x in range(pop.numSubPop())]
                # migrate into toSubPop and replace individuals
                num_migrants[toSubPop] -= pop.subPopSize(toSubPop)
        #
        if toSubPop is None:
            # Now, we need to select specified number of individuals from the subpopulations
            indexes = []
            for sz, sp in zip(num_migrants, list(range(pop.numSubPop()))):
                indexes.extend(list(range(pop.subPopBegin(sp), pop.subPopBegin(sp) + sz)))
            sample = pop.extractIndividuals(indexes=indexes)
            sample.mergeSubPops(name=self.subPopName)
            pop.addIndFrom(sample)
        else:
            # replacing an existing subpopulation
            # adjust the size of existing subpopulations
            # copy individuals that will be in admixed population
            sz_before = list(pop.subPopSizes())
            sz_after = [x + y for x,y in zip(pop.subPopSizes(), num_migrants)]
            pop.resize(sz_after, propagate=True)
            for sp in range(pop.numSubPop() - 1, -1, -1):
                if sz_after[sp] > sz_before[sp]:
                    pop.splitSubPop(sp, [sz_before[sp], sz_after[sp] - sz_before[sp]])
                else:
                    pop.splitSubPop(sp, [sz_after[sp], 0])
            pop.mergeSubPops([x+x + 1 for x in range(pop.numSubPop()//2)], toSubPop=toSubPop+toSubPop+1)
            pop.mergeSubPops([toSubPop, toSubPop + 1])
        return True


class _OutOfAfricaModel_event(EventBasedModel):
    '''A dempgraphic model for the CHB, CEU, and YRI populations, as defined in
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
        ops=[],
        infoFields=[],
        outcome=['AF', 'EU', 'AS'],
        scale=1
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
        respectively. Additional operators could be added to ``ops``. Information
        fields required by these operators should be passed to ``infoFields``. If 
        a scaling factor ``scale`` is specified, all population sizes and
        generation numbers will be divided by a factor of ``scale``. This demographic
        model by default returns all populations (``AF``, ``EU``, ``AS``) but
        you can choose to keep only selected subpopulations using parameter
        ``outcome`` (e.g. ``outcome=['EU', 'AS']``).

        This model merges all subpopulations if it is applied to an initial 
        population with multiple subpopulation.
        '''
        #
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        #
        if isinstance(outcome, str):
            outcome = [outcome]
        #
        scale = float(scale)
        EventBasedModel.__init__(self, 
            T = int(T0/scale),
            N0 = (int(N_A/scale), 'Ancestral'),
            infoFields = 'migrate_to',
            events = [
                # resize ancestral to AF
                ResizeEvent(at= -int(T_AF/scale), subPops='Ancestral',
                    sizes=int(N_AF/scale), names='AF'),
                # split B from AF
                SplitEvent(at=-int(T_B/scale), subPops='AF',
                    sizes=(1.0, int(N_B/scale)), names=('AF', 'B')),
                # migration between AF and B
                DemographicEvent(
                    ops=Migrator(rate=[
                        [m_AF_B, 0],
                        [0, m_AF_B]]), 
                    begin=-int(T_B/scale), 
                    end=-int(T_EU_AS/scale)
                ),
                # split EU AS from B
                SplitEvent(at=-int(T_EU_AS/scale), subPops='B',
                    sizes=(int(N_EU0/scale), int(N_AS0/scale)),
                    names=('EU', 'AS')),
                # exponential growth 
                ExpansionEvent(begin=-int(T_EU_AS/scale),
                    subPops=('EU', 'AS'), 
                    rates=(r_EU*scale, r_AS*scale)),
                # migration between AF, EU and AS
                DemographicEvent(
                    ops=Migrator(rate=[
                        [0, m_AF_EU, m_AF_AS],
                        [m_AF_EU, 0, m_EU_AS],
                        [m_AF_AS, m_EU_AS, 0]
                        ]), 
                    begin=-int(T_EU_AS/scale)
                ),
                # keep only selected populations
                ResizeEvent(at=-1, 
                    sizes=[float('AF' in outcome),
                        float('EU' in outcome),
                        float('AS' in outcome)],
                    removeEmptySubPops=True)
            ]
        )


class _OutOfAfricaModel_outcome(MultiStageModel):
    '''A dempgraphic model for the CHB, CEU, and YRI populations, as defined in
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
        ops=[],
        infoFields=[],
        outcome=['AF', 'EU', 'AS'],
        scale=1
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
        respectively. Additional operators could be added to ``ops``. Information
        fields required by these operators should be passed to ``infoFields``. If 
        a scaling factor ``scale`` is specified, all population sizes and
        generation numbers will be divided by a factor of ``scale``. This demographic
        model by default returns all populations (``AF``, ``EU``, ``AS``) but
        you can choose to keep only selected subpopulations using parameter
        ``outcome`` (e.g. ``outcome=['EU', 'AS']``).

        This model merges all subpopulations if it is applied to an initial 
        population with multiple subpopulation.
        '''
        #
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        #
        if isinstance(outcome, str):
            outcome = [outcome]
        #
        final_subpops = [None, None, None]
        for (idx, name) in enumerate(['AF', 'EU', 'AS']):
            if name not in outcome:
                final_subpops[idx] = 0
            #
        if 0 in final_subpops:
            finalStage = [
                InstantChangeModel(T=1,
                    N0=final_subpops,
                    removeEmptySubPops=True)
            ]
        else:
            finalStage = [] 
        # for python 2.x and 3.x compatibility
        scale = float(scale)
        MultiStageModel.__init__(self, [
            InstantChangeModel(
                T=int((T0-T_B)/scale),
                N0=(int(N_A/scale), 'Ancestral'),
                # change population size twice, one at T_AF, one at T_B
                G=[int((T0-T_AF)/scale)],
                NG=[(int(N_AF/scale), 'AF')] 
            ),
            #
            # at T_B, split to population B from subpopulation 1
            InstantChangeModel(
                T=int((T_B - T_EU_AS)/scale),
                # change population size twice, one at T_AF, one at T_B
                N0=[None, (int(N_B/scale), 'B')],
                ops=Migrator(rate=[
                    [m_AF_B, 0],
                    [0, m_AF_B]])
                ),
            ExponentialGrowthModel(
                T=int(T_EU_AS/scale),
                N0=[None, 
                    # split B into EU and AS at the beginning of this
                    # exponential growth stage
                    [(int(N_EU0/scale), 'EU'), (int(N_AS0/scale), 'AS')]],
                r=[0, r_EU*scale, r_AS*scale],
                infoFields='migrate_to',
                ops=Migrator(rate=[
                    [0, m_AF_EU, m_AF_AS],
                    [m_AF_EU, 0, m_EU_AS],
                    [m_AF_AS, m_EU_AS, 0]
                    ])
                ),
            ] + finalStage, ops=ops, infoFields=infoFields
        )

OutOfAfricaModel = _OutOfAfricaModel_outcome

class _SettlementOfNewWorldModel_event(EventBasedModel):
    '''A dempgraphic model for settlement of the new world of Americans, as defined
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
        f_MX=0.48,
        ops=[],
        infoFields=[],
        outcome='MXL',
        scale=1
        ):
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
        with ``f_MX`` proportion for ``MX``. The Mexican American<F19> sample could
        be sampled from the last single population. Additional operators could
        be added to ``ops``. Information fields required by these operators 
        should be passed to ``infoFields``. If a scaling factor ``scale``
        is specified, all population sizes and generation numbers will be divided by
        a factor of ``scale``. This demographic model by default only returns the
        mixed Mexican America model (``outputcom='MXL'``) but you can specify any
        combination of ``AF``, ``EU``, ``AS``, ``MX`` and ``MXL``.

        This model merges all subpopulations if it is applied to an initial population
        with multiple subpopulation.
        '''
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        #
        scale = float(scale)
        if isinstance(outcome, str):
            outcome = [outcome]
        EventBasedModel.__init__(self,
            T = int(T0/scale),
            N0 = (int(N_A/scale), 'Ancestral'),
            infoFields = 'migrate_to',
            events = [
                # resize ancestral to AF
                ResizeEvent(at= -int(T_AF/scale), subPops='Ancestral',
                    sizes=int(N_AF/scale), names='AF'),
                # split B from AF
                SplitEvent(at=-int(T_B/scale), subPops='AF',
                    sizes=(1.0, int(N_B/scale)), names=('AF', 'B')),
                # migration between AF and B
                DemographicEvent(
                    ops=Migrator(rate=[
                        [m_AF_B, 0],
                        [0, m_AF_B]]), 
                    begin=-int(T_B/scale), 
                    end=-int(T_EU_AS/scale)
                ),
                # split EU AS from B
                SplitEvent(at=-int(T_EU_AS/scale), subPops='B',
                    sizes=(int(N_EU0/scale), int(N_AS0/scale)),
                    names=('EU', 'AS')),
                # exponential growth 
                ExpansionEvent(begin=-int(T_EU_AS/scale),
                    subPops=('EU', 'AS'), 
                    rates=(r_EU*scale, r_AS*scale)),
                # split MX from AS
                SplitEvent(at=-int(T_MX/scale),
                    subPops='AS',
                    sizes=(1.0, int(N_MX0/scale)),
                    names=('AS', 'MX')),
                # expansion of MX
                ExpansionEvent(begin=-int(T_MX/scale),
                    subPops='MX',
                    rates=r_MX),
                # migration between AF, EU and AS
                DemographicEvent(
                    ops=Migrator(rate=[
                        [0, m_AF_EU, m_AF_AS],
                        [m_AF_EU, 0, m_EU_AS],
                        [m_AF_AS, m_EU_AS, 0]
                        ],
                    subPops=['AF', 'EU', 'AS']), 
                    begin=-int(T_EU_AS/scale),
                ),
                # admixture
                AdmixtureEvent(
                    at=-1,
                    subPops=['EU', 'MX'],
                    sizes=(1-f_MX, f_MX),
                    name='MXL'),
                # keep only selected populations
                ResizeEvent(at=-1, 
                    sizes=[float('AF' in outcome),
                        float('EU' in outcome),
                        float('AS' in outcome),
                        float('MX' in outcome),
                        float('MXL' in outcome)],
                    removeEmptySubPops=True)
            ]
        )

class _SettlementOfNewWorldModel_outcome(MultiStageModel):
    '''A dempgraphic model for settlement of the new world of Americans, as defined
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
        f_MX=0.48,
        ops=[],
        infoFields=[],
        outcome='MXL',
        scale=1
        ):
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
        with ``f_MX`` proportion for ``MX``. The Mexican American<F19> sample could
        be sampled from the last single population. Additional operators could
        be added to ``ops``. Information fields required by these operators 
        should be passed to ``infoFields``. If a scaling factor ``scale``
        is specified, all population sizes and generation numbers will be divided by
        a factor of ``scale``. This demographic model by default only returns the
        mixed Mexican America model (``outputcom='MXL'``) but you can specify any
        combination of ``AF``, ``EU``, ``AS``, ``MX`` and ``MXL``.

        This model merges all subpopulations if it is applied to an initial population
        with multiple subpopulation.
        '''
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        # try to figure out how to mix two populations
        N_EU=int(N_EU0*math.exp(r_EU*T_EU_AS))
        N_MX=int(N_MX0*math.exp(r_MX*T_MX))
        #
        # for python 2.x and 3.x compatibility
        if isinstance(outcome, str):
            outcome = [outcome]
        if 'MXL' in outcome:
            # with admixture
            final_subpops = [None, None, None, None, None]
            for (idx, name) in enumerate(['AF', 'EU', 'AS', 'MX', 'MXL']):
                if name not in outcome:
                    final_subpops[idx] = 0
            #
            admixtureStage = [
                AdmixtureModel(T=1,
                    N0=[None, None, None, None],
                    # mixing European and Mexican population
                    model=['HI', 1, 3, 1-f_MX, 'MXL']),
                InstantChangeModel(T=1,
                    N0=final_subpops,
                    removeEmptySubPops=True)
                ]
        else:
            final_subpops = [None, None, None, None]
            for (idx, name) in enumerate(['AF', 'EU', 'AS', 'MX']):
                if name not in outcome:
                    final_subpops[idx] = 0
            admixtureStage = [
                InstantChangeModel(T=1,
                    N0=final_subpops,
                    removeEmptySubPops=True)
                ]
        #
        scale = float(scale)
        MultiStageModel.__init__(self, [
            InstantChangeModel(
                T=int((T0-T_B)/scale),
                N0=(int(N_A/scale), 'Ancestral'),
                # change population size twice, one at T_AF, one at T_B
                G=[int((T0-T_AF)/scale)],
                NG=[(int(N_AF/scale), 'AF')] 
            ),
            #
            # at T_B, split to population B from subpopulation 1
            InstantChangeModel(
                T=int((T_B - T_EU_AS)/scale),
                # change population size twice, one at T_AF, one at T_B
                N0=[None, (int(N_B/scale), 'B')],
                ops=Migrator(rate=[
                    [m_AF_B, 0],
                    [0, m_AF_B]])
                ),
            ExponentialGrowthModel(
                T=int((T_EU_AS - T_MX)/scale),
                N0=[None,
                    # split B into EU and AS at the beginning of this
                    # exponential growth stage
                    [(int(N_EU0/scale), 'EU'), (int(N_AS0/scale), 'AS')]],
                r=[0, r_EU*scale, r_AS*scale],
                infoFields='migrate_to',
                ops=Migrator(rate=[
                    [0, m_AF_EU, m_AF_AS],
                    [m_AF_EU, 0, m_EU_AS],
                    [m_AF_AS, m_EU_AS, 0]
                    ])
                ),
            ExponentialGrowthModel(T=int(T_MX/scale),
                N0=[None,
                    # initial population size has to be calculated
                    # because we only know the final population size of
                    # EU and AS
                    None,
                    # split MX from AS
                    [(None, 'AS'), (int(N_MX0//scale), 'MX')]],
                r=[0, r_EU*scale, r_AS*scale, r_MX*scale],
                infoFields='migrate_to',
                ops=Migrator(rate=[
                    [0, m_AF_EU, m_AF_AS],
                    [m_AF_EU, 0, m_EU_AS],
                    [m_AF_AS, m_EU_AS, 0]
                    ],
                    # the last MX population does not involve in 
                    # migration
                    subPops=[0, 1, 2],
                    toSubPops=[0, 1, 2])
                )
            ] + admixtureStage,
            ops=ops, infoFields=infoFields
        )


SettlementOfNewWorldModel = _SettlementOfNewWorldModel_outcome

class CosiModel(MultiStageModel):
    '''A dempgraphic model for Africa, Asia and Europe, as described in 
    Schaffner et al, Genome Research, 2005, and implemented in the coalescent
    simulator cosi.'''
    def __init__(self,
        T0,
        N_A=12500,
        N_AF=24000,
        N_OoA=7700,
        N_AF1=100000,
        N_AS1=100000,
        N_EU1=100000,
        T_AF=17000,
        T_OoA=3500,
        T_EU_AS=2000,
        T_AS_exp=400,
        T_EU_exp=350,
        T_AF_exp=200,
        F_OoA=0.085,
        F_AS=0.067,
        F_EU=0.020,
        F_AF=0.020,
        m_AF_EU=0.000032,
        m_AF_AS=0.000008,
        ops=[],
        infoFields=[],
        scale=1
        ):
        '''Counting **backward in time**, this model evolves a population for a
        total of ``T0`` generations. The ancient population ``Ancestral`` started
        at size ``N_Ancestral`` and expanded at ``T_AF`` generations from now,
        to pop ``AF`` with size ``N_AF``. The Out of Africa population split from
        the ``AF`` population at ``T_OoA`` generations ago. The ``OoA`` population
        split into two subpopulations ``AS`` and ``EU`` but keep the same size.
        At the generations of ``T_EU_exp``, ``T_AS_exp``, and ``T_AF_exp`` ago,
        three populations expanded to modern population sizes of ``N_AF1``, 
        ``N_AS1`` and ``N_EU1`` exponentially, respectively. Migrations are
        allowed between ``AF`` and ``EU`` populations
        with rate ``m_AF_EU``, and between ``AF`` and ``AS`` with rate ``m_AF_AS``.

        Four bottlenecks happens in the ``AF``, ``OoA``, ``EU`` and ``AS`` populations.
        They are supposed to happen 200 generations after population split and last
        for 200 generations. The intensity is parameterized in F, which is number
        of generations devided by twice the effective size during bottleneck.
        So the bottleneck size is 100/F. 

        This model merges all subpopulations if it is applied to a population with
        multiple subpopulation. Although parameters are configurable, we assume
        the order of events so dramatically changes of parameters might need
        to errors.  If a scaling factor ``scale`` is specified, all population
        sizes and generation numbers will be divided by, and migration rates
        will be multiplied by a factor of ``scale``.
         '''
        #
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        #
        if T_AF < T_OoA or T_OoA < T_EU_AS or T_EU_AS < T_AS_exp or T_AS_exp < T_EU_exp or T_EU_exp < T_AF_exp:
            raise ValueError('Specified parameters change the order of events to the model.')
        # for python 2.x and 3.x compatibility
        scale = float(scale)
        # by model
        N_AS = N_OoA
        N_EU = N_OoA
        r_AS = math.log(1.0*N_AS1/N_AS)/T_AS_exp
        r_EU = math.log(1.0*N_EU1/N_EU)/T_EU_exp
        r_AF = math.log(1.0*N_AF1/N_AF)/T_AF_exp
        migr = Migrator(rate=[
                    [0, m_AF_AS, m_AF_EU],
                    [m_AF_AS, 0, 0],
                    [m_AF_EU, 0, 0]
                    ])
        MultiStageModel.__init__(self, [
            InstantChangeModel(
                # constant population size before the first one to expand
                T=int((T0 - T_EU_AS)/scale),
                N0=(int(N_A/scale), 'Ancestral'),
                G=[ int((T0-T_AF)/scale), 
                    int((T0-T_OoA)/scale),
                    int((T0-T_OoA+200)/scale),
                    int((T0-T_OoA+400)/scale)],
                NG=[
                    # population size incrase to N_AF
                    (int(N_AF/scale), 'Africa'),
                    # at T_B, split to population B from subpopulation 1
                    [(int(N_AF/scale), 'Africa'), (int(N_OoA/scale), 'Out Of Africa')],
                    [int(100./F_AF/scale), int(100./F_OoA/scale)], # bottleneck
                    [int(N_AF/scale), int(N_OoA/scale)],  # recover
                    ]
                ),
            InstantChangeModel(
                # constant population size before the first one to expand
                T=int((T_EU_AS - T_AS_exp)/scale),
                N0=[int(N_AF/scale), [(int(N_OoA/scale), 'Asian'), (int(N_OoA/scale), 'Europe')]], # split
                G=[ int(200./scale), int(400./scale)],
                NG=[
                    [int(N_AF/scale), int(100./F_AS/scale), int(100./F_EU/scale)],
                    [int(N_AF/scale), int(N_OoA/scale), int(N_OoA/scale)] # recover
                    ],
                ops=migr,
                ),
            # AS expend 
            ExponentialGrowthModel(
                T=int((T_AS_exp-T_EU_exp)/scale),
                N0=[None, (None, 'Modern Asian'), None],
                r=[0, r_AS*scale, 0],
                infoFields='migrate_to',
                ops=migr),
            # EU expand
            ExponentialGrowthModel(T=int((T_EU_exp-T_AF_exp)/scale),
                N0=[None, None, (None, 'Modern Europe')],
                r=[0, r_AS*scale, r_EU*scale],
                infoFields='migrate_to',
                ops=migr),
            # AF expand
            ExponentialGrowthModel(T=int(T_AF_exp/scale),
                N0=[(None, 'Modern Africa'), None, None],
                NT=[int(N_AF1/scale), int(N_AS1/scale), int(N_EU1/scale)],
                infoFields='migrate_to',
                ops=migr),
            ]
        )


if __name__ == '__main__':
    # exponential
    ExponentialGrowthModel(10, 100, 1000).plot(title='Basic exponential growth')
    ExponentialGrowthModel(10, (100, 200), (1000, 2000)).plot(title='Subpop exp growth')
    ExponentialGrowthModel(10, (100, 200), r=0.01).plot(title='No NT exp growth')
    ExponentialGrowthModel(10, (100, 200), r=(0.01, 0.2)).plot('ExpDemo.png',
        title='Exponential population growth model')
    # linear
    LinearGrowthModel(10, 100, 1000).plot(title='Basic linear growth')
    LinearGrowthModel(10, (100, 200), (1000, 2000)).plot(title='Subpop linear growth')
    LinearGrowthModel(10, (100, 200), r=0.01).plot(title='No NT linear growth')
    LinearGrowthModel(10, (100, 200), r=(0.1, 0.2)).plot('LinearDemo.png',
        title='Linear population growth model')
    # instant
    InstantChangeModel(10, 100, 5, 1000).plot()
    InstantChangeModel(10, (100, 200), 5, (1000, 2000)).plot()
    InstantChangeModel(10, 100, [5, 8], [500, 100]).plot()
    InstantChangeModel(50, 100, [5, 8, 20], [[500, 200], [100, 100], [1000, 2000]]).plot('InstantDemo.png')
    #
    # multi-stage model
    MultiStageModel([
        InstantChangeModel(10, 100, 5, 1000),
        ExponentialGrowthModel(20, 1000, 2000)
        ]).plot('MultiStageDemo.png')
    # Out Of Africa Model
    OutOfAfricaModel(10000).plot('OutOfAfrica.png')
    # Settlement of New World
    SettlementOfNewWorldModel(10000).plot('SettlementOfNewWorld.png')
    # Cosi model 
    CosiModel(20000).plot('Cosi.png')

