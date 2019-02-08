#!/usr/bin/env python

#
# $File: dynamicNumOff.py $
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

class RandomNumOff:
    # a demographic model
    def __init__(self):
        self.numOff = []
    
    def getNumOff(self):
        # return the pre-simulated number of offspring as a generator function
        for item in self.numOff:
            yield item
    
    def __call__(self, pop):
        # define __call__ so that a RandomNumOff object is callable.
        #
        # Each male produce from 1 to 3 offspring. For large population, get the
        # number of males instead of checking the sex of each individual
        self.numOff = [random.randint(1, 3) for ind in pop.individuals() if ind.sex() == sim.MALE]
        # return the total population size
        print('{} mating events with number of offspring {}'.format(len(self.numOff), self.numOff))
        return sum(self.numOff)


pop = sim.Population(10)

# create a demogranic model
numOffModel = RandomNumOff()

pop.evolve(
    preOps=sim.InitSex(),
    matingScheme=sim.RandomMating(
        # the model will be called before mating to deteremine
        # family and population size
        subPopSize=numOffModel,
        # the getNumOff function (generator) returns number of offspring
        # for each mating event
        numOffspring=numOffModel.getNumOff
    ),
    gen=3
)


