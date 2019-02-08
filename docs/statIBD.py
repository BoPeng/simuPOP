#!/usr/bin/env python

#
# $File: statIBD.py $
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
simuOpt.setOptions(alleleType='lineage')
import simuPOP as sim
pop = sim.Population([500], loci=[1]*100)
pop.evolve(
    initOps=[
        sim.InitLineage(),
        sim.InitSex(),
        sim.InitGenotype(freq=[0.2]*5),
    ],
    preOps=[
        sim.Stat(inbreeding=sim.ALL_AVAIL, popSize=True, step=10),
        sim.PyEval(r'"gen %d: IBD freq %.4f, IBS freq %.4f, est: %.4f\n" % '
            '(gen, sum(IBD_freq.values()) /len(IBD_freq), '
            ' sum(IBS_freq.values()) /len(IBS_freq), '
            ' 1 - (1-1/(2.*popSize))**gen)', step=10)
    ],
    matingScheme=sim.RandomMating(),
    gen = 100
)

