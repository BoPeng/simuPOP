#!/usr/bin/env python

#
# $File: PyMutator.py $
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
import random
def incAllele(allele):
    return allele + random.randint(1, 5)

pop = sim.Population(size=1000, loci=[20])
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.RandomMating(),
    postOps=sim.PyMutator(func=incAllele, rates=[1e-4, 1e-3],
            loci=[2, 10]),
    gen = 1000
)
# count the average number tandem repeats at both loci
def avgAllele(pop, loc):
    ret = 0
    for ind in pop.individuals():
        ret += ind.allele(loc, 0) + ind.allele(loc, 1)
    return ret / (pop.popSize() * 2.)

print('Average number of repeats at two loci are %.2f and %.2f.' % \
    (avgAllele(pop, 2), avgAllele(pop, 10)))

