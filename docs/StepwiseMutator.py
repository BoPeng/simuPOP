#!/usr/bin/env python

#
# $File: StepwiseMutator.py $
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
pop = sim.Population(size=1000, loci=[1, 1])
pop.evolve(
    # all start from allele 50
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq= [0]*50 + [1])
    ],
    matingScheme=sim.RandomMating(),
    preOps=[
        sim.StepwiseMutator(rates=1e-3, loci=0),
        sim.StepwiseMutator(rates=1e-3, incProb=0.6, loci=1,
            mutStep=(sim.GEOMETRIC_DISTRIBUTION, 0.2)),
    ],
    gen=100
)
# count the average number tandem repeats at both loci
cnt0 = cnt1 = 0
for ind in pop.individuals():
    cnt0 += ind.allele(0, 0) + ind.allele(0, 1)
    cnt1 += ind.allele(1, 0) + ind.allele(1, 1)

print('Average number of repeats at two loci are %.2f and %.2f.' % \
    (cnt0/2000., cnt1/2000.))

