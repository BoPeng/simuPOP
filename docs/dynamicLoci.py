#!/usr/bin/env python

#
# $File: dynamicLoci.py $
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
pop = sim.Population(100, loci=[10], infoFields='fitness')

def mostPopular(pop):
    sim.stat(pop, alleleFreq=sim.ALL_AVAIL)
    freq = [pop.dvars().alleleFreq[x][1] for x in range(pop.totNumLoci())]
    max_freq = max(freq)
    pop.dvars().selLoci = (freq.index(max_freq), max_freq)
    return [freq.index(max_freq)]

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.6, 0.4]),
    ],
    preOps=[
        sim.MaSelector(fitness=[1, 0.9, 0.8], loci=mostPopular),
        sim.PyEval(r"'gen=%d, select against %d with frequency %.2f\n' % (gen, selLoci[0], selLoci[1])"),
    ],
    matingScheme=sim.RandomMating(),
    gen=10,
)

