#!/usr/bin/env python

#
# $File: PyParentsChooser.py $
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
from random import randint
def randomChooser(pop, subPop):
    males = []
    females = []
    # identify males and females in each social rank
    for rank in range(3):
        males.append([x for x in pop.individuals(subPop) \
            if x.sex() == sim.MALE and x.rank == rank])
        females.append([x for x in pop.individuals(subPop) \
            if x.sex() == sim.FEMALE and x.rank == rank])
    #
    while True:
        # choose a rank randomly
        rank = int(pop.individual(randint(0, pop.subPopSize(subPop) - 1), subPop).rank)
        yield males[rank][randint(0, len(males[rank]) - 1)], \
            females[rank][randint(0, len(females[rank]) - 1)]

def setRank(rank):
    'The rank of offspring can increase or drop to zero randomly'
    # only use rank of the father
    return (rank[0] + randint(-1, 1)) % 3

pop = sim.Population(size=[1000, 2000], loci=1, infoFields='rank')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitInfo(lambda : randint(0, 2), infoFields='rank')
    ],
    matingScheme=sim.HomoMating(
        sim.PyParentsChooser(randomChooser),
        sim.OffspringGenerator(ops=[
            sim.MendelianGenoTransmitter(),
            sim.PyTagger(setRank),
            ])
    ),
    gen = 5
)
