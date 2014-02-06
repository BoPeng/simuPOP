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
]

import sys
import time

from simuOpt import simuOptions

from simuPOP import moduleInfo, MALE, FEMALE, Population, PointMutator, getRNG,\
    ALL_AVAIL, PyOperator, stat

from simuPOP.utils import migrIslandRates, migrHierarchicalIslandRates,
    migrSteppingStoneRates

def migr2DSteppingStoneRates(r, m, n, diagonal=False, circular=False):
    '''migration rate matrix for 2D stepping stone model, with or without
    diagonal neighbors (4 or 8 neighbors for central patches). The boundaries
    are connected if circular is True. Otherwise individuals from corner and
    bounary patches will migrate to their neighbors with higher probability.
    '''
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

