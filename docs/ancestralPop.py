#!/usr/bin/env python

#
# $File: ancestralPop.py $
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
pop = sim.Population(500, loci=1, ancGen=2)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    matingScheme = sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=0, begin=-3),
        sim.PyEval(r"'%.3f\n' % alleleFreq[0][0]", begin=-3)
    ],
    gen = 20
)
# information
pop.ancestralGens()
pop.popSize(ancGen=1)
pop.setVirtualSplitter(sim.SexSplitter())
# number of males in the current and parental generation
pop.subPopSize((0,0)), pop.subPopSize((0,0), ancGen=1)
# start from current generation
for i in range(pop.ancestralGens(), -1, -1):
  pop.useAncestralGen(i)
  sim.stat(pop, alleleFreq=0)
  print('%d   %.3f' % (i, pop.dvars().alleleFreq[0][0]))

# restore to the current generation  
pop.useAncestralGen(0)  

