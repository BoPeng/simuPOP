#!/usr/bin/env python

#
# $File: PyQuanTrait.py $
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
pop = sim.Population(size=5000, loci=2, infoFields=['qtrait1', 'qtrait2', 'age'])
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[40]))
def qtrait(geno, age):
    'Return two traits that depends on genotype and age'
    return random.normalvariate(age * sum(geno), 10), random.randint(0, 10*sum(geno))

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.2, 0.8]),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        # use random age for simplicity
        sim.InitInfo(lambda:random.randint(20, 75), infoFields='age'),
        sim.PyQuanTrait(loci=(0,1), func=qtrait, infoFields=['qtrait1', 'qtrait2']),
        sim.Stat(meanOfInfo=['qtrait1'], subPops=[(0, sim.ALL_AVAIL)],
            vars='meanOfInfo_sp'),
        sim.PyEval(r"'Mean of trait1: %.3f (age < 40), %.3f (age >=40)\n' % "
            "(subPop[(0,0)]['meanOfInfo']['qtrait1'], subPop[(0,1)]['meanOfInfo']['qtrait1'])"),
    ],
    gen = 5
)


