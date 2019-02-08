#!/usr/bin/env python

#
# $File: MlSelector.py $
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
pop = sim.Population(size=10000, loci=2, infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[.5, .5])
    ],
    preOps=[
        sim.MlSelector([
            sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}),
            sim.MapSelector(loci=1, fitness={(0,0):1, (0,1):0.9, (1,1):.8}),
            ], mode = sim.ADDITIVE, reps=0),
        sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}, reps=1),
        sim.MapSelector(loci=1, fitness={(0,0):1, (0,1):0.9, (1,1):.8}, reps=2)
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
         sim.Stat(alleleFreq=[0,1]),
         sim.PyEval(r"'REP %d:\t%.3f\t%.3f\t' % (rep, alleleFreq[0][1], alleleFreq[1][1])"),
         sim.PyOutput('\n', reps=-1),
    ],
    gen = 5
)

