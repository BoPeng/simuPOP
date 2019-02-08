#!/usr/bin/env python

#
# $File: sexMode.py $
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
def checkSexMode(ms):
    '''Check the assignment of sex to offspring'''
    pop = sim.Population(size=[40])
    pop.evolve(initOps=sim.InitSex(), matingScheme=ms, gen=1)
    # return individual sex as a string
    return ''.join(['M' if ind.sex() == sim.MALE else 'F' for ind in pop.individuals()])

# Case 1: sim.NO_SEX (all male, sim.RandomMating will not continue)
checkSexMode(sim.RandomMating(sexMode=sim.NO_SEX))
# Case 2: sim.RANDOM_SEX (sim.Male/Female with probability 0.5)
checkSexMode(sim.RandomMating(sexMode=sim.RANDOM_SEX))
# Case 3: sim.PROB_OF_MALES (Specify probability of male)
checkSexMode(sim.RandomMating(sexMode=(sim.PROB_OF_MALES, 0.8)))
# Case 4: sim.NUM_OF_MALES (Specify number of male in each family)
checkSexMode(sim.RandomMating(numOffspring=3, sexMode=(sim.NUM_OF_MALES, 1)))
# Case 5: sim.NUM_OF_FEMALES (Specify number of female in each family)
checkSexMode(sim.RandomMating(
    numOffspring=(sim.UNIFORM_DISTRIBUTION, 4, 6),
    sexMode=(sim.NUM_OF_FEMALES, 2))
)
# Case 6: sim.SEQUENCE_OF_SEX
checkSexMode(sim.RandomMating(
    numOffspring=4, sexMode=(sim.SEQUENCE_OF_SEX, sim.MALE, sim.FEMALE))
)
# Case 7: sim.GLOBAL_SEQUENCE_OF_SEX
checkSexMode(sim.RandomMating(
    numOffspring=3, sexMode=(sim.GLOBAL_SEQUENCE_OF_SEX, sim.MALE, sim.FEMALE))
)
# Case 8: A generator function
def sexFunc():
    i = 0
    while True:
        i += 1
        if i % 2 == 0:
            yield sim.MALE
        else:
            yield sim.FEMALE

checkSexMode(sim.RandomMating(numOffspring=3, sexMode=sexFunc))

