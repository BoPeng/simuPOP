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
simuPOP utilities.

This module provides some commonly used operators
and format conversion utilities.

"""

__all__ = [
    'migrIslandRates',
    'migrHierarchicalIslandRates',
    'migrSteppingStoneRates',
    'migr2DSteppingStoneRates',
    'MexicanAmerican_Model',
]

import sys
import time
import math

from simuOpt import simuOptions

from simuPOP import moduleInfo, MALE, FEMALE, StopEvolution, Population, \
    ALL_AVAIL, Stat, stat, PyEval, RandomMating, migrate, InitSex, PyOperator

from simuPOP.utils import migrIslandRates, migrHierarchicalIslandRates, \
    migrSteppingStoneRates

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

class MexicanAmerican_Model:
    def __init__(self, T0=8000, N_A=7300, T_AF=220000//25, N_AF=12300,
        T_B=140000//25, N_B=2100, 
        T_EU_AS=26400//25, 
        N_EU0=1500, N_AS0=590, 
        N_EU=16970, N_AS=29147,
        m_AF_B=0.00025,
        m_EU_AS=0.000135,
        m_AF_EU=0.00003,
        m_AF_AS=0.000019):
        self.N_A = N_A
        self.T_AF = T0 - T_AF
        self.N_AF = N_AF
        self.T_B = T0 - T_B
        self.N_B = N_B
        self.T_EU_AS = T0 - T_EU_AS
        self.N_EU0 = N_EU0
        self.N_AS0 = N_AS0
        self.N_EU = N_EU
        self.N_AS = N_AS
        self.m_AF_B = m_AF_B
        self.m_EU_AS = m_EU_AS
        self.m_AF_EU = m_AF_EU
        self.m_AF_AS = m_AF_AS
        #
        self.last_gen = T0

    def __call__(self, gen, pop=None):
        if gen < self.T_AF:
            return self.N_A
        elif gen < self.T_B:
            return self.N_AF
        elif gen == self.T_B:
            # before mating, split 
            pop.resize(self.N_AF + self.N_B, propagate=True)
            pop.splitSubPop(0, [self.N_AF, self.N_B], ['AF', 'B'])
            return pop.subPopSizes()
        elif gen <= self.T_EU_AS:
            # migration between AF and B
            migrate(pop, rate=[
                [0, self.m_AF_B],
                [self.m_AF_B, 0]])
            if gen == self.T_EU_AS:
                # split at T_EUAS
                pop.resize([self.N_AF, self.N_EU0 + self.N_AS0], propagate=True)
                pop.splitSubPop(1, [self.N_EU0, self.N_AS0], ['EU', 'AS'])
                print(pop.subPopSizes())
                return pop.subPopSizes()
            else:
                # migration will change population size, this will reset 
                # it to constant
                return [self.N_AF, self.N_B]
        elif gen < self.last_gen:
            # exponentially expand
            # migration between three populations
            migrate(pop, rate=[
                [0, self.m_AF_EU, self.m_AF_AS],
                [self.m_EU_AS, 0, self.m_AF_EU],
                [self.m_AF_AS, self.m_AF_EU, 0]
            ])
            #
            # sz[0] keep contant
            if gen == self.last_gen - 1:
                return [self.N_AF, self.N_EU, self.N_AS]
            else:
                # exponential growth with known N0, Nt, t
                x = gen - self.T_EU_AS 
                t = self.last_gen - self.T_EU_AS - 1
                return [self.N_AF, 
                    int(math.exp((x*math.log(self.N_EU) + (t-x)*math.log(self.N_EU0))/t)),
                    int(math.exp((x*math.log(self.N_AS) + (t-x)*math.log(self.N_AS0))/t))]
        else:
            # step 0 -- 
            # step 1 ---
            # ...
            # step last_gen - 1  ==> evolves last_gen generations
            # step last_gen <- stop at the beginning of this generation.
            return []
            #raise StopIteration('Ending at the beginning of generation %d' % gen)

def runDemo(model):
    def reportPopSize(pop):
        stat(pop, popSize=True)
        if 'last_size' not in pop.vars() or pop.dvars().last_size != pop.dvars().subPopSize:
            print('%d: %s' % (pop.dvars().gen, pop.dvars().subPopSize))
            pop.dvars().last_size = pop.dvars().subPopSize
        return True

    pop = Population(model(0), infoFields='migrate_to')
    return pop.evolve(
        preOps=[
            InitSex()
        ],
        matingScheme=RandomMating(subPopSize=model),
        postOps=PyOperator(reportPopSize)
    )

if __name__ == '__main__':
    print runDemo(MexicanAmerican_Model())
