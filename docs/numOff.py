#!/usr/bin/env python

#
# $File: numOff.py $
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
def checkNumOffspring(numOffspring, ops=[]):
    '''Check the number of offspring for each family using
       information field father_idx
    '''
    pop = sim.Population(size=[30], loci=1, infoFields=['father_idx', 'mother_idx'])
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitGenotype(freq=[0.5, 0.5]),
        ],
        matingScheme=sim.RandomMating(ops=[
            sim.MendelianGenoTransmitter(),
            sim.ParentsTagger(),
            ] + ops,
            numOffspring=numOffspring),
        gen=1)
    # get the parents of each offspring
    parents = [(x, y) for x, y in zip(pop.indInfo('mother_idx'),
        pop.indInfo('father_idx'))]
    # Individuals with identical parents are considered as siblings.
    famSize = []
    lastParent = (-1, -1)
    for parent in parents:
        if parent == lastParent:
            famSize[-1] += 1
        else:
            lastParent = parent
            famSize.append(1)
    return famSize

# Case 1: produce the given number of offspring
checkNumOffspring(numOffspring=2)
# Case 2: Use a Python function
import random
def func(gen):
    return random.randint(5, 8)

checkNumOffspring(numOffspring=func)
# Case 3: A geometric distribution
checkNumOffspring(numOffspring=(sim.GEOMETRIC_DISTRIBUTION, 0.3))
# Case 4: A Possition distribution
checkNumOffspring(numOffspring=(sim.POISSON_DISTRIBUTION, 1.6))
# Case 5: A Binomial distribution
checkNumOffspring(numOffspring=(sim.BINOMIAL_DISTRIBUTION, 0.1, 10))
# Case 6: A uniform distribution
checkNumOffspring(numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 6))
# Case 7: With selection on offspring
checkNumOffspring(numOffspring=8,
    ops=[sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):0.8, (1,1):0.5})])

