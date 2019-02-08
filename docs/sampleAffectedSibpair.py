#!/usr/bin/env python

#
# $File: sampleAffectedSibpair.py $
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
from simuPOP.sampling import indexToID
pop = sim.Population(size=15, loci=5, infoFields=['father_idx', 'mother_idx'], ancGen=2)
pop.evolve(
    preOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.7, 0.3]),
    ],
    matingScheme=sim.RandomMating(numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
        ops=[sim.MendelianGenoTransmitter(), sim.ParentsTagger()]),
    postOps=sim.MaPenetrance(loci=3, penetrance=(0.1, 0.4, 0.7)),
    gen = 5
)
indexToID(pop, reset=True)
# three information fields were added
print(pop.infoFields())
# save this population for future use
pop.save('log/pedigree.pop')

from simuPOP.sampling import drawAffectedSibpairSample
pop = sim.loadPopulation('log/pedigree.pop')
sample = drawAffectedSibpairSample(pop, families=2)

