#!/usr/bin/env python

#
# $File: statNeTemporal.py $
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit https://github.com/BoPeng/simuPOP
# for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (Bo.Peng@bcm.edu)
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
# the user's guide (https://github.com/BoPeng/simuPOP/manual) for a detailed
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
            vars=['Ne_waples89_P1', 'Ne_tempoFS_P1'], step=20),
        sim.PyEval(r'"Waples Ne: %.1f (%.1f - %.1f), TempoFS: '
            r'%.1f (%.1f - %.1f), at generation %d\n" % '
            'tuple(Ne_waples89_P1 + Ne_tempoFS_P1 + [gen])', step=20)
    ],
    matingScheme=sim.RandomMating(),
    gen = 101
)
