#!/usr/bin/env python

#
# $File: extract.py $
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
pop = sim.Population(size=[200, 200], loci=[5, 5], infoFields='age')
sim.initGenotype(pop, genotype=range(10))
sim.initInfo(pop, lambda: random.randint(0,75), infoFields='age')
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[20, 60]))
# remove individuals
pop.removeIndividuals(indexes=range(0, 300, 10))
print(pop.subPopSizes())
# remove individuals using IDs
pop.setIndInfo([1, 2, 3, 4], field='age')
pop.removeIndividuals(IDs=[2, 4], idField='age')
# remove indiviuals using a filter function
sim.initSex(pop)
pop.removeIndividuals(filter=lambda ind: ind.sex() == sim.MALE)
print([pop.individual(x).sex() for x in range(8)])
#
# remove subpopulation
pop.removeSubPops(1)
print(pop.subPopSizes())
# remove virtual subpopulation (people with age between 20 and 60)
pop.removeSubPops([(0, 1)])
print(pop.subPopSizes())
# extract another virtual subpopulation (people with age greater than 60)
pop1 = pop.extractSubPops([(0,2)])
sim.dump(pop1, structure=False, max=10)

