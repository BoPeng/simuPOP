#!/usr/bin/env python

#
# $File: IfElse.py $
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
simu = sim.Simulator(
    sim.Population(size=1000, loci=1),
    rep=4)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5]),
        sim.PyExec('below40, above60 = 0, 0')
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=0),
        sim.IfElse('alleleFreq[0][1] < 0.4',
            sim.PyExec('below40 += 1')),
        sim.IfElse('alleleFreq[0][1] > 0.6',
            sim.PyExec('above60 += 1')),
        sim.IfElse('len(alleleFreq[0]) == 1',
            sim.PyExec('stoppedAt = gen')),
        sim.TerminateIf('len(alleleFreq[0]) == 1')
    ]
)
for pop in simu.populations():
    print('Overall: %4d, below 40%%: %4d, above 60%%: %4d' % \
        (pop.dvars().stoppedAt, pop.dvars().below40, pop.dvars().above60))


