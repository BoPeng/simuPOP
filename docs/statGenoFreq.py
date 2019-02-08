#!/usr/bin/env python

#
# $File: statGenoFreq.py $
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

# This script is an example in the simuPOP user's guide. Please refer to
# the user's guide (http://simupop.sourceforge.net/manual) for a detailed
# description of this example.
#

import simuPOP as sim
pop = sim.Population(100, loci=[1, 1, 1], lociNames=['A', 'X', 'Y'],
    chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y])
sim.initGenotype(pop, freq=[0.01, 0.05, 0.94])
sim.stat(pop, genoFreq=['A', 'X']) # both loci indexes and names can be used.
print('Available genotypes on autosome:', list(pop.dvars().genoFreq[0].keys()))
for i in range(3):
    for j in range(3):
        print('%d-%d: %.3f' % (i, j, pop.dvars().genoFreq[0][(i,j)]))

print('Genotype frequency on chromosome X:\n', \
    '\n'.join(['%s: %.3f' % (x,y) for x,y in pop.dvars().genoFreq[1].items()]))

