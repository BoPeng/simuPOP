#!/usr/bin/env python

#
# $File: InheritTagger.py $
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
pop = sim.Population(size=[1000]*10, loci=1, infoFields='x')
# tag the first individual of each subpopulation.
for sp in range(pop.numSubPop()):
    pop.individual(0, sp).x = 1

pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.RandomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.InheritTagger(mode=sim.MAXIMUM, infoFields='x'),
    ]),
    postOps=[
        sim.Stat(sumOfInfo='x', vars=['sumOfInfo_sp']),
        sim.PyEval(r'", ".join(["%3d" % subPop[i]["sumOfInfo"]["x"] for i in range(10)])+"\n"'),
    ],
    gen = 5
)

