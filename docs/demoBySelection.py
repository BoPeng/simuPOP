#!/usr/bin/env python

#
# $File: demoBySelection.py $
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
def demo(pop):
    return int(pop.popSize() * 1.05)

pop = sim.Population(size=10000, loci=1)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.7, 0.3])
    ],
    preOps=[
        sim.Stat(popSize=True),
        sim.PyEval(r'"%d %s --> " % (gen, subPopSize)'),
        sim.ResizeSubPops(0, proportions=[0.5], at=2),
        sim.MaPenetrance(loci=0, penetrance=[0.01, 0.2, 0.6], begin=4),
        sim.DiscardIf('ind.affected()', exposeInd='ind', begin=4),
        sim.Stat(popSize=True),
        sim.PyEval(r'"%s --> " % subPopSize'),
    ],
    matingScheme=sim.RandomMating(subPopSize=demo),
    postOps=[
        sim.Stat(popSize=True),
        sim.PyEval(r'"%s\n" % subPopSize')
    ],
    gen = 6
)

