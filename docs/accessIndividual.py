#!/usr/bin/env python

#
# $File: accessIndividual.py $
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
# create a sim.population with two generations. The current generation has values
# 0-9 at information field x, the parental generation has values 10-19.
pop = sim.Population(size=[5, 5], loci=[2, 3], infoFields='x', ancGen=1)
pop.setIndInfo(range(10, 20), 'x')
pop1 = pop.clone()
pop1.setIndInfo(range(10), 'x')
pop.push(pop1)
#
ind = pop.individual(5)       # using absolute index
ind.x
ind.x       # the same as ind.x
# use a for loop, and relative index
for idx in range(pop.subPopSize(1)):
    print(pop.individual(idx, 1).x)

# It is usually easier to use an iterator
for ind in pop.individuals(1):
    print(ind.x)

# Access individuals in VSPs
pop.setVirtualSplitter(sim.InfoSplitter(cutoff=[3, 7, 17], field='x'))
for ind in pop.individuals([1, 1]):
    print(ind.x)

# Access all individuals in all ancestral generations
print([ind.x for ind in pop.allIndividuals()])
# or only specified subpopulations or ancestral generations
print([ind.x for ind in pop.allIndividuals(subPops=[(0,2), (1,3)], ancGens=1)])

# Access individuals in ancetral generations
pop.ancestor(5, 1).x        # absolute index
pop.ancestor(0, 1, 1).x     # relative index
# Or make ancestral generation the current generation and use 'individual'
pop.useAncestralGen(1)
pop.individual(5).x         # absolute index
pop.individual(0, 1).x      # relative index
# 'ancestor' can still access the 'present' (generation 0) generation
pop.ancestor(5, 0).x
# access individual by ID
pop.addInfoFields('ind_id')
sim.tagID(pop)
[int(ind.ind_id) for ind in pop.individuals()]
# access individual by ID. Note that individual 12 is in the parental generation
pop.indByID(12).x

