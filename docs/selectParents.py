#!/usr/bin/env python

#
# $File: selectParents.py $
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
pop = sim.Population(4000, loci=1, infoFields='fitness')
simu = sim.Simulator(pop, rep=3)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    preOps=sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):0.98, (1,1):0.97}),
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=0, step=10),
        sim.PyEval("'Gen:%3d ' % gen", reps=0, step=10),
        sim.PyEval(r"'%.3f\t' % alleleFreq[0][1]", step=10),
        sim.PyOutput('\n', reps=-1, step=10)
    ],
    gen = 50
)

