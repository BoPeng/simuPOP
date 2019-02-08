#!/usr/bin/env python

#
# $File: randomSeed.py $
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
import random
def simulate():
    pop = sim.Population(1000, loci=10, infoFields='age')
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitGenotype(freq=[0.5, 0.5]),
            sim.InitInfo(lambda: random.randint(0, 10), infoFields='age')
        ],
        matingScheme=sim.RandomMating(),
        finalOps=sim.Stat(alleleFreq=0),
        gen=100
    )
    return pop.dvars().alleleFreq[0][0]

seed = sim.getRNG().seed()
random.seed(seed)
print('%.4f' % simulate())
# will yield different result
print('%.4f' % simulate())
sim.setRNG(seed=seed)
random.seed(seed)
# will yield identical result because the same seeds are used
print('%.4f' % simulate())

