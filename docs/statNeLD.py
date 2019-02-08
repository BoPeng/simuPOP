#!/usr/bin/env python

#
# $File: statNeLD.py $
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
pop = sim.Population([2000], loci=[1]*50)
pop.setVirtualSplitter(sim.RangeSplitter([0, 500]))
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.005]*4 + [0.015]*2 + [0.25, 0.7]),
    ],
    preOps=[
        sim.Stat(effectiveSize=sim.ALL_AVAIL, subPops=[(0,0)], 
            vars='Ne_LD', step=20),
        sim.PyEval(r'"LD Ne (gen %d): %.1f (%.1f - %.1f)'
            r', %.1f (%.1f - %.1f, adjusted)\n" % '
            'tuple([gen] + Ne_LD[0.] + Ne_LD[0.02])',
            step=20)
    ],
    matingScheme=sim.RandomMating(),
    gen = 101
)

