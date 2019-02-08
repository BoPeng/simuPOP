#!/usr/bin/env python

#
# $File: monogamous.py $
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
pop = sim.Population(20, infoFields=['father_idx', 'mother_idx'])
pop.evolve(
    initOps=sim.InitSex(sex=(sim.MALE, sim.FEMALE)),
    matingScheme=sim.MonogamousMating(
        numOffspring=2,
        sexMode=(sim.NUM_OF_MALES, 1),
        ops=[
            sim.MendelianGenoTransmitter(),
            sim.ParentsTagger(),
        ],
    ),
    gen = 5
)
[ind.sex() for ind in pop.individuals()]
[int(ind.father_idx) for ind in pop.individuals()]
[int(ind.mother_idx) for ind in pop.individuals()]
# count the number of distinct parents
len(set(pop.indInfo('father_idx')))
len(set(pop.indInfo('mother_idx')))

