#!/usr/bin/env python

#
# $File: statHWE.py $
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
pop = sim.Population([1000], loci=1)
pop.setVirtualSplitter(sim.ProportionSplitter([0.4, 0.4, 0.2]))
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(genotype=[0,0], subPops=[(0,0)]),
        sim.InitGenotype(genotype=[0,1], subPops=[(0,1)]),
        sim.InitGenotype(genotype=[1,1], subPops=[(0,2)]),
    ],
    preOps=[
        sim.Stat(HWE=0, genoFreq=0),
        sim.PyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
            'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])'),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(HWE=0, genoFreq=0),
        sim.PyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
            'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])'),
    ],
    gen = 1
)

