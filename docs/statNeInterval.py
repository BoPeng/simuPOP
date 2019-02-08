#!/usr/bin/env python

#
# $File: statNeInterval.py $
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
        sim.InitGenotype(freq=[0.3, 0.7]),
        sim.Stat(effectiveSize=range(50), subPops=[(0,0)],
            vars='Ne_temporal_base'),
    ],
    preOps=[
        sim.Stat(effectiveSize=range(50), subPops=[(0,0)], 
            vars='Ne_waples89_P1', step=20),
        sim.Stat(effectiveSize=range(50), subPops=[(0,0)], step=20,
            suffix='_i', vars=['Ne_temporal_base', 'Ne_waples89_P1']),
        sim.PyEval(r'"Waples Ne (till %d): %.1f (%.1f - %.1f), '
            r'(interval) %.1f (%.1f - %.1f)\n" % '
            'tuple([gen] + Ne_waples89_P1 + Ne_waples89_P1_i)',
            step=20)
    ],
    matingScheme=sim.RandomMating(),
    gen = 101
)

