#!/usr/bin/env python

#
# $File: vspDuringMatingSelector.py $
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
pop = sim.Population(size=[5000, 5000], loci=1, infoFields='fitness')
pop.setVirtualSplitter(sim.SexSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[.5, .5])
    ],
    matingScheme=sim.RandomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.MaSelector(loci=0, fitness=[1, 1, 0.98], subPops=[(0,0), (1,1)]),
        sim.MaSelector(loci=0, fitness=[1, 0.99, 0.98], subPops=[(0,1), (1,0)]),
        ]),
    postOps=[
        sim.Stat(alleleFreq=[0], subPops=[(sim.ALL_AVAIL, sim.ALL_AVAIL)],
            vars='alleleFreq_sp', step=50),
        sim.PyEval(r"'%.4f\t%.4f\t%.4f\t%.4f\n' % "
            "tuple([subPop[x]['alleleFreq'][0][1] for x in ((0,0),(0,1),(1,0),(1,1))])",
            step=50)
    ],
    gen=151
)

