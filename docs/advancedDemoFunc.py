#!/usr/bin/env python

#
# $File: advancedDemoFunc.py $
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
def demo(gen, pop):
    if gen < 2:
        return 1000 + 100 * gen
    if gen == 2:
        # this happens right before mating at generation 2
        size = pop.popSize()
        pop.splitSubPop(0, [size // 2, size - size//2]) 
    # for generation two and later
    return [x + 50 * gen for x in pop.subPopSizes()]

pop = sim.Population(1000)
pop.evolve(
    preOps=[
        sim.Stat(popSize=True),
        sim.PyEval(r'"Gen %d:\t%s (before mating)\t" % (gen, subPopSize)')
    ],
    matingScheme=sim.RandomSelection(subPopSize=demo),
    postOps=[
        sim.Stat(popSize=True),
        sim.PyEval(r'"%s (after mating)\n" % subPopSize')
    ],
    gen = 5
)

