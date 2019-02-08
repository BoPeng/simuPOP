#!/usr/bin/env python

#
# $File: PySelector.py $
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
pop = sim.Population(size=2000, loci=[1]*2, infoFields=['fitness', 'smoking'])
s1 = .02
s2 = .03
# the second parameter gen can be used for varying selection pressure
def sel(geno, smoking):
    #     BB  Bb   bb
    # AA  1   1    1
    # Aa  1   1-s1 1-s2
    # aa  1   1    1-s2
    #
    # geno is (A1 A2 B1 B2)
    if geno[0] + geno[1] == 1 and geno[2] + geno[3] == 1:
        v = 1 - s1  # case of AaBb
    elif geno[2] + geno[3] == 2:
        v = 1 - s2  # case of ??bb
    else:                
        v = 1       # other cases
    if smoking:
        return v * 0.9
    else:
        return v

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[.5, .5])
    ],
    preOps=sim.PySelector(loci=[0, 1], func=sel),
    matingScheme=sim.RandomMating(),
    postOps=[
        # set smoking status randomly
        sim.InitInfo(lambda : random.randint(0,1), infoFields='smoking'),
        sim.Stat(alleleFreq=[0, 1], step=20),
        sim.PyEval(r"'%.4f\t%.4f\n' % (alleleFreq[0][1], alleleFreq[1][1])", step=20)
    ],
    gen = 50
)

