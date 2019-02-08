#!/usr/bin/env python

#
# $File: statNeDemographic.py $
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

import simuOpt
simuOpt.setOptions(alleleType='lineage', quiet=True)
import simuPOP as sim
pop = sim.Population([2000], loci=[1]*3,
    chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y])
pop.setVirtualSplitter(sim.SexSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.3, 0.7]),
    ],
    preOps=[
        sim.Stat(effectiveSize=range(3), subPops=[0, (0,0), (0,1)],
            vars='Ne_demo_base_sp'),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(effectiveSize=range(3), subPops=[0, (0,0), (0,1)],
            vars='Ne_demo_sp'),
        sim.PyEval(r'"Demographic Ne: %.1f (auto), %.1f (X), %.1f (Y), '
            r'Males: %.1f, %.1f, %.1f, Females: %.1f, %.1f, %.1f\n"'
            '% tuple([subPop[0]["Ne_demo"][x] for x in (0, 1, 2)] + '
            '[subPop[(0,0)]["Ne_demo"][x] for x in (0, 1, 2)] + '
            '[subPop[(0,1)]["Ne_demo"][x] for x in (0, 1, 2)])')
    ],
    gen = 5
)

