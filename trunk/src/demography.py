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
    'CosiModel'
]

import sys
import time
import math

from simuOpt import simuOptions

from simuPOP import Population, PyEval, RandomSelection, \
    ALL_AVAIL, Stat, stat, Migrator, InitSex, PyOperator

from simuPOP.utils import migrIslandRates, migrHierarchicalIslandRates, \
    migrSteppingStoneRates

from collections import OrderedDict
from types import NoneType

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
    def __init__(self, initSize=[], numGens=-1, ops=[], infoFields=[]):
        '''Set attributes ``init_size``, ``info_fields``, and ``num_gens``
        to a demographic model. If size for a subpopulation for
        initial population is a list, the initial subpopulation will be split to
        multiple subpopulations. For example, ``N0=[A, [B,C]]`` is a 3-subpop
        model where the last two subpopulation will be split (and resize if needed)
        from the second subpopulation of the initial subpopulation (which should
        have two subpopulations). The population size can be an integer for
        fixed population size, None for the size of the population when the 
        demographic model is first applied to, and a float number representing
        the proportion of individuals for the whole or corresponding subpopulation.
        A ``None`` value will be assigned to attribute ``init_size`` in such a 
        case because the initial population size is determined dynamically. In
        addition, whenever a population size is allowed, a tuple of ``(size, name)`` is
        acceptable, which assigns name to the corresponding subpopulation.
        One or more operators (e.g. a migration operator or a terminator) could
        be passed (parameter ``ops``) and will be applied to the population.
        The demographic model will return ``[]`` (which will effectively terminate
        the evolutioonary process) if any of the operator returns ``False``.
        Information fields required by these operators should be passed to
        ``infoFields``. '''
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
        self.num_gens = numGens
        if isinstance(ops, (tuple, list)):
            self.ops = list(ops)
        else:
            self.ops = [ops]

    def _reset(self):
        if hasattr(self, '_start_gen'):
            del self._start_gen

    def _isNamedSize(self, x):
        return isinstance(x, tuple) and len(x) == 2 and \
            isinstance(x[1], str) and self._isSize(x[0])

    def _isSize(self, x):
        return isinstance(x, (int, float, NoneType))

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
            elif all([self._isNamedSize(x) for x in size]):
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
                for x, y in zip(named_size, pop.subPopSizes()):
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
                        total = 0
                        for z in x:
                            if isinstance(z[0], int):
                                total += z[0]
                            elif isinstance(z[0], float):
                                total += int(z[0]*y)
                            elif z[0] is None:
                                total += y
                            else:
                                raise ValueError('Invalid size %s' % named_size)
                        new_size.append(total)
                        new_names.append('')
                # resize and rename
                pop.resize(new_size, propagate=True)
                for idx, name in enumerate(new_names):
                    if name != '':
                        pop.setSubPopName(name, idx)
                # handle split
                indexes = [i for i, x in enumerate(named_size) if not self._isNamedSize(x)]
                indexes.reverse()
                for idx in indexes:
                    names = [x[1] for x in named_size[idx]]
                    # new sizes, handle None and float
                    new_size = [x[0] for x in named_size[idx]]
                    # convert None to size
                    new_size = [pop.subPopSize(idx) if x is None else x for x in new_size]
                    # convert float to int
                    new_size = [int(x*pop.subPopSize(idx)) if isinstance(x, float) else x for x in new_size]
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
        stat(pop, popSize=True)
        gen = pop.dvars().gen
        if 'last_size' not in pop.vars() or pop.dvars().last_size != pop.dvars().subPopSize:
            print('%d: %s' % (gen, 
                ', '.join(
                    ['%d%s' % (x, ' (%s)' % y if y else '') for x, y in \
                        zip(pop.dvars().subPopSize, pop.subPopNames())])
                ))
            pop.dvars().last_size = pop.dvars().subPopSize
        #
        if self.draw_figure:
            sz = 0
            for idx, (s, n) in enumerate(zip(pop.dvars().subPopSize, pop.subPopNames())):
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
                pop.dvars().last_size = pop.dvars().subPopSize
        return True

    def plot(self, filename='', title=''):
        '''Evolve a population using this demographic model. Print population
        size changes duringe evution. If a filename is specified and if matplotlib
        is available, draw a figure to depict the demographic model. An optional
        title could be specified to the figure.
        '''
        if not self.init_size:
            raise ValueError('Specific self does not have a valid initial population size')
        if filename and not has_plotter:
            raise RuntimeError('This function requires module numpy and matplotlib')
        self.draw_figure = filename and has_plotter
        self.pop_regions = OrderedDict()
        self.pop_base = OrderedDict()
        #
        self._reset()
        pop = Population(self.init_size, infoFields=self.info_fields, ploidy=1)
        pop.evolve(
            matingScheme=RandomSelection(subPopSize=self),
            postOps=PyOperator(self._recordPopSize),
            gen=self.num_gens
        )
        del pop.vars()['last_size']
        self._recordPopSize(pop)
        # 
        if self.draw_figure:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            for name, region in self.pop_regions.items():
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
        for op in self.ops:
            if not op.apply(pop):
                self._reset()
                return []
        return pop.subPopSizes()
        

class BaseGrowthModel(BaseDemographicModel):
    '''A base model for population growth.'''
    def __init__(self, T, N0, ops=[], infoFields=[]):
        '''An population growth model that evolves a population from size ``N0``
        to ``NT`` for ``T`` generations with rate ``r``, under either an exponential,
        linear or instant population growth model. '''
        if not isinstance(T, int):
            raise ValueError('Number of generations must be an integer')
        BaseDemographicModel.__init__(self, initSize=N0,
            numGens=T, ops=ops, infoFields=infoFields)
        

class ExponentialGrowthModel(BaseGrowthModel):
    '''A model for exponential population growth. Four parameters
    N0, T, NT, and r are allowed but only three are needed.'''
    def __init__(self, T, N0, NT=None, r=None, ops=[], infoFields=[]):
        '''An exponential population growth model that evolves a population
        from size ``N0`` to ``NT`` for ``T`` generations with rate ``r``.
        ``N0``, ``NT`` and ``r`` can be a list of population sizes or growth
        rates for multiple subpopulations. The initial population will be
        resized to ``N0`` (split if necessary). The model will automatically 
        determine ``r`` or ``NT`` if one of them is unspecified. Optionally,
        one or more operators (e.g. a migrator) ``ops`` can be applied to 
        population. '''
        BaseGrowthModel.__init__(self, T, N0, ops, infoFields)
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
        if not BaseGrowthModel.__call__(self, pop):
            return []
        if self._gen == self.num_gens:
            return []
        else:
            return [self._expIntepolate(n0, nt, self.num_gens, self._gen)
                for (n0, nt) in zip(self.init_size, self.NT)]


class LinearGrowthModel(BaseGrowthModel):
    '''A model for linear population growth. Four parameters
    N0, T, NT, and r are allowed but only three are needed.'''
    def __init__(self, T, N0=[], NT=None, r=None, ops=[], infoFields=[]):
        '''An linear population growth model that evolves a population
        from size ``N0`` to ``NT`` for ``T`` generations with ``r*N0`` 
        individuals added at each generation. ``N0``, ``NT`` and ``r``
        can be a list of population sizes or growth rates for multiple 
        subpopulations. The initial population will be
        resized to ``N0`` (split if necessary). The model will automatically 
        determine ``r`` or ``NT`` if one of them is unspecified. Optionally,
        one or more operators (e.g. a migrator) ``ops`` can be applied to 
        population. '''
        BaseGrowthModel.__init__(self, T, N0, ops, infoFields)
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
        if not BaseGrowthModel.__call__(self, pop):
            return []
        if self._gen == self.num_gens:
            return []
        else:
            return [self._linearIntepolate(n0, nt, self.num_gens, self._gen)
                for (n0, nt) in zip(self.init_size, self.NT)]    

class InstantChangeModel(BaseGrowthModel):
    '''A model for instant population growth model.'''
    def __init__(self, T, N0, G=[], NG=[], ops=[], infoFields=[]):
        '''An instant population growth model that evolves a population
        from size ``N0`` to ``NT`` for ``T`` generations with population
        size changes at generation ``G`` to ``NT``. If ``G`` is a list,
        multiple population size changes are allowed. In that case, a list
        (or a nested list) of population size should be provided to parameter
        ``NT``. This model supports population merge to and split from a single
        population. Optionally, one or more operators (e.g. a migrator) ``ops``
        can be applied to population. Required information fields by these
        operators should be passed to parameter ``infoFields``.'''
        BaseGrowthModel.__init__(self, T, N0, ops, infoFields)
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
        if not BaseGrowthModel.__call__(self, pop):
            return []
        if self._gen in self.G:
            sz = self.NG[self.G.index(self._gen)]
            self._fitToSize(pop, sz)
        return pop.subPopSizes()


class MultiStageModel(BaseDemographicModel):
    '''A multi-stage demographic model connects a few demographic models 
    '''
    def __init__(self, models, ops=[], infoFields=[]):
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
        BaseDemographicModel.__init__(self, initSize=models[0].init_size,
            numGens=total_gens, ops=ops, infoFields=flds+infoFields)
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
        if not BaseDemographicModel.__call__(self, pop):
            return []
        # at the end?
        if self.models[self._model_idx].num_gens == self._gen:
            self._model_idx += 1
            self._start_gen = pop.dvars().gen
        if self._model_idx == len(self.models):
            self._reset()
            return []
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
        ops=[],
        infoFields=[]
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
        fields required by these operators should be passed to ``infoFields``.

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
                infoFields='migrate_to',
                ops=Migrator(rate=[
                    [0, m_AF_EU, m_AF_AS],
                    [m_EU_AS, 0, m_AF_EU],
                    [m_AF_AS, m_AF_EU, 0]
                    ])
                ),
            ], ops=ops, infoFields=infoFields
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
        f_MX=0.48,
        ops=[],
        infoFields=[]
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
        with ``f_MX`` proportion for ``MX``. The Mexican American sample could
        be sampled from the last single population. Additional operators could
        be added to ``ops``. Information fields required by these operators 
        should be passed to ``infoFields``.

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
                infoFields='migrate_to',
                ops=Migrator(rate=[
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
                infoFields='migrate_to',
                ops=Migrator(rate=[
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
            )], ops=ops, infoFields=infoFields
        )

class CosiModel(MultiStageModel):
    '''A dempgrahic model for Africa, Asia and Europe, as described in 
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
        infoFields=[]
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
        to errors. '''
        #
        if T0 < T_AF:
            raise ValueError('Length of evolution T0=%d should be more than T_AF=%d' % (T0, T_AF))
        #
        if T_AF < T_OoA or T_OoA < T_EU_AS or T_EU_AS < T_AS_exp or T_AS_exp < T_EU_exp or T_EU_exp < T_AF_exp:
            raise ValueError('Specified parameters change the order of events to the model.')
        # by model
        N_AS = N_OoA
        N_EU = N_OoA
        r_AS = math.log(1.0*N_AS1/N_AS)/T_AS_exp
        r_EU = math.log(1.0*N_EU1/N_EU)/T_EU_exp
        r_AF = math.log(1.0*N_AF1/N_AF)/T_AF_exp
        MultiStageModel.__init__(self, [
            InstantChangeModel(
                # constant population size before the first one to expand
                T=T0-T_AS_exp,
                N0=(N_A, 'Ancestral'),
                # change population size twice, one at T_AF, one at T_B
                G=[ T0-T_AF, 
                    T0-T_OoA,
                    T0-T_OoA+200,
                    T0-T_OoA+400,
                    T0-T_EU_AS,
                    T0-T_EU_AS+200,
                    T0-T_EU_AS+400],
                NG=[
                    # population size incrase to N_AF
                    (N_AF, 'Africa'),
                    # at T_B, split to population B from subpopulation 1
                    [(N_AF, 'Africa'), (N_OoA, 'Out Of Africa')],
                    [int(100./F_AF), int(100./F_OoA)], # bottleneck
                    [N_AF, N_OoA],  # recover
                    [N_AF, [(N_OoA, 'Asian'), (N_OoA, 'Europe')]], # split
                    [N_AF, int(100./F_AS), int(100./F_EU)],
                    [N_AF, N_OoA, N_OoA] # recover
                    ]
                ),
            # AS expend 
            ExponentialGrowthModel(
                T=T_AS_exp-T_EU_exp,
                N0=[N_AF, (N_OoA, 'Modern Asian'), N_OoA],
                r=[0, r_AS, 0],
                infoFields='migrate_to',
                ops=Migrator(rate=[
                    [0, m_AF_AS, m_AF_EU],
                    [m_AF_AS, 0, 0],
                    [m_AF_EU, 0, 0]
                    ])
                ),
            # EU expand
            ExponentialGrowthModel(T=T_EU_exp-T_AF_exp,
                N0=[N_AF,
                    int(N_AS*((1.+r_AS)**(T_AS_exp-T_EU_exp))),
                    (N_OoA, 'Modern Europe')],
                r=[0, r_AS, r_EU],
                infoFields='migrate_to',
                ops=Migrator(rate=[
                    [0, m_AF_AS, m_AF_EU],
                    [m_AF_AS, 0, 0],
                    [m_AF_EU, 0, 0]
                    ])
                ),
            # AF expand
            ExponentialGrowthModel(T=T_AF_exp,
                N0=[(N_AF, 'Modern Africa'),
                    int(N_AS*((1.+r_AS)**(T_AS_exp-T_AF_exp))),
                    int(N_EU*((1.+r_EU)**(T_EU_exp-T_AF_exp)))],
                r=[r_AF, r_AS, r_EU],
                infoFields='migrate_to',
                ops=Migrator(rate=[
                    [0, m_AF_AS, m_AF_EU],
                    [m_AF_AS, 0, 0],
                    [m_AF_EU, 0, 0]
                    ])
                ),
            ]
        )


if __name__ == '__main__':
    # exponential
    ExponentialGrowthModel(10, 100, 1000).plot()
    ExponentialGrowthModel(10, (100, 200), (1000, 2000)).plot()
    ExponentialGrowthModel(10, (100, 200), r=0.01).plot()
    ExponentialGrowthModel(10, (100, 200), r=(0.01, 0.2)).plot()
    ExponentialGrowthModel(10, (100, 200), r=(0.01, 0.2)).plot('ExpDemo.png',
        title='Exponential population growth model')
    # linear
    LinearGrowthModel(10, 100, 1000).plot()
    LinearGrowthModel(10, (100, 200), (1000, 2000)).plot()
    LinearGrowthModel(10, (100, 200), r=0.01).plot()
    LinearGrowthModel(10, (100, 200), r=(0.1, 0.2)).plot()
    LinearGrowthModel(10, (100, 200), r=(0.1, 0.2)).plot('LinearDemo.png')
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
        ]).plot()
    MultiStageModel([
        InstantChangeModel(10, 100, 5, 1000),
        ExponentialGrowthModel(20, 1000, 2000)
        ]).plot('MultiStageDemo.png')

    # Out Of Africa Model
    OutOfAfricaModel(10000).plot()
    OutOfAfricaModel(10000).plot('OutOfAfrica.png')
    # Settlement of New World
    #
    SettlementOfNewWorldModel(10000).plot()
    SettlementOfNewWorldModel(10000).plot('SettlementOfNewWorld.png')
    # 
    CosiModel(20000).plot()
    CosiModel(20000).plot('Cosi.png')

