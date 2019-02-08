#!/usr/bin/env python

#
# $File: MaSelectorHaploid.py $
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
pop = sim.Population(size=10000, ploidy=1, loci=[1,1], infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[.5, .5])
    ],
    # fitness values for AB, Ab, aB and ab
    preOps=sim.MaSelector(loci=[0,1], fitness=[1, 1, 1, 0.95]),
    matingScheme=sim.RandomSelection(),
    postOps=[
        sim.Stat(haploFreq=[0, 1], step=25),
        sim.PyEval(r"'%.3f\t%.3f\t%.3f\t%.3f\n' % (haploFreq[(0,1)][(0,0)],"
                "haploFreq[(0,1)][(0,1)], haploFreq[(0,1)][(1,0)],"
                "haploFreq[(0,1)][(1,1)])", step=25)
    ],
    gen = 100
)

