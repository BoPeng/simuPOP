#!/usr/bin/env python

#
# $File: DiscardIf.py $
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
pop = sim.Population(size=500, loci=1)
pop.setVirtualSplitter(sim.ProductSplitter([
    sim.AffectionSplitter(),
    sim.RangeSplitter([[0,500], [500, 1000]]),
    ])
)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5]),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            sim.MendelianGenoTransmitter(),
            sim.MaPenetrance(loci=0, penetrance=[0, 0.01, 0.1]),
            sim.DiscardIf(True, subPops=[
                (0, 'Unaffected, Range [0, 500)'),
                (0, 'Affected, Range [500, 1000)')])
        ],
        subPopSize=1000,
    ),
    gen = 1
)
sim.stat(pop, numOfAffected=True)
print(pop.dvars().numOfAffected, pop.dvars().numOfUnaffected)

