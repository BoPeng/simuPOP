#!/usr/bin/env python

#
# $File: HeteroMatingWeight.py $
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
pop = sim.Population(size=[1000], loci=2,
    infoFields='mark')
pop.setVirtualSplitter(sim.RangeSplitter([[0, 500], [200, 1000]]))

pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.HeteroMating([
        sim.RandomMating(subPops=0, weight=-0.5,
            ops=[sim.InfoExec('mark=0'), sim.MendelianGenoTransmitter()]),
        sim.RandomMating(subPops=[(0, 0)], weight=2,
            ops=[sim.InfoExec('mark=1'), sim.MendelianGenoTransmitter()]),
        sim.RandomMating(subPops=[(0, 1)], weight=3,
            ops=[sim.InfoExec('mark=2'), sim.MendelianGenoTransmitter()])
    ]),
    gen = 10
)
marks = list(pop.indInfo('mark'))
marks.count(0.)
marks.count(1.)
marks.count(2.)

