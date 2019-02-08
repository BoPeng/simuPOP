#!/usr/bin/env python

#
# $File: SavePopulation.py $
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
simu = sim.Simulator(sim.Population(100, loci=2),
    rep=5)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.2, 0.8])
    ],
    matingScheme=sim.RandomMating(),
    postOps=sim.SavePopulation(output="!'snapshot_%d_%d.pop' % (rep, gen)",
            step = 10),
    gen = 50
)

