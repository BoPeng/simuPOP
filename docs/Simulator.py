#!/usr/bin/env python

#
# $File: Simulator.py $
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
pop = sim.Population(100, loci=10)
# five copies of the same population
simu = sim.Simulator(pop, rep=5)
simu.numRep()
# evolve for ten generations and save the populations
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.3, 0.7])
    ],
    matingScheme=sim.RandomMating(),
    finalOps=sim.SavePopulation('!"pop%d.pop"%rep'),
    gen=10
)
# load the population and create another Simulator
simu = sim.Simulator([sim.loadPopulation('pop%d.pop' % x) for x in range(5)])
# continue to evolve
simu.evolve(
    matingScheme=sim.RandomMating(),
    gen=10
)
# print out allele frequency
for pop in simu.populations():
    sim.stat(pop, alleleFreq=0)
    print('%.2f' % pop.dvars().alleleFreq[0][0])

# get a population
pop = simu.extract(0)
simu.numRep()

