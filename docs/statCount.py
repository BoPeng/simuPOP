#!/usr/bin/env python

#
# $File: statCount.py $
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
pop = sim.Population(10000, loci=1)
pop.setVirtualSplitter(sim.CombinedSplitter(
    [sim.SexSplitter(), sim.AffectionSplitter()]))
sim.initSex(pop)
sim.initGenotype(pop, freq=[0.2, 0.8])
sim.maPenetrance(pop, loci=0, penetrance=[0.1, 0.2, 0.5])
# Count sim.population size
sim.stat(pop, popSize=True, subPops=[(0, 0), (0, 2)])
# popSize is the size of two VSPs, does not equal to total sim.population size.
# Because two VSPs overlap (all males and all unaffected), popSize can be
# greater than real sim.population size.
print(pop.dvars().subPopSize, pop.dvars().popSize)
# print popSize of each virtual subpopulation.
sim.stat(pop, popSize=True, subPops=[(0, 0), (0, 2)], vars='popSize_sp')
# Note the two ways to access variable in (virtual) subpopulations.
print(pop.dvars((0,0)).popSize, pop.dvars().subPop[(0,2)]['popSize'])
# Count number of male (should be the same as the size of VSP (0,0).
sim.stat(pop, numOfMales=True)
print(pop.dvars().numOfMales)
# Count the number of affected and unaffected male individual
sim.stat(pop, numOfMales=True, subPops=[(0, 2), (0, 3)], vars='numOfMales_sp')
print(pop.dvars((0,2)).numOfMales, pop.dvars((0,3)).numOfMales)
# or number of affected male and females
sim.stat(pop, numOfAffected=True, subPops=[(0, 0), (0, 1)], vars='numOfAffected_sp')
print(pop.dvars((0,0)).numOfAffected, pop.dvars((0,1)).numOfAffected)
# These can also be done using a sim.ProductSplitter...
pop.setVirtualSplitter(sim.ProductSplitter(
    [sim.SexSplitter(), sim.AffectionSplitter()]))
sim.stat(pop, popSize=True, subPops=[(0, x) for x in range(4)])
# counts for male unaffected, male affected, female unaffected and female affected
print(pop.dvars().subPopSize)

