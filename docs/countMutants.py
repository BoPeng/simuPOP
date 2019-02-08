#!/usr/bin/env python

#
# $File: countMutants.py $
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
from collections import defaultdict
# count number of mutants at each locus
counter = defaultdict(int)
def countMutants(mutants):
    global counter
    for line in mutants.split('\n'):
        # a trailing \n will lead to an empty string
        if not line:  
            continue
        (gen, loc, ploidy, a1, a2, id) = line.split('\t')
        counter[int(loc)] += 1

pop = sim.Population([5000]*3, loci=[2,1,1], infoFields='ind_id',
    chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5]),
        sim.IdTagger(),
    ],
    preOps=[
        sim.KAlleleMutator(rates=[0.001] + [0.01]*3,
            loci=range(4), k=100, output=countMutants),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            sim.IdTagger(),
            sim.MendelianGenoTransmitter()
        ]),
    gen = 10
)
print(counter.items())

