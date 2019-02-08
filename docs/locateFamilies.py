#!/usr/bin/env python

#
# $File: locateFamilies.py $
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
pop = sim.Population(1000, ancGen=-1, infoFields=['ind_id', 'father_id', 'mother_id'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
    ],
    matingScheme=sim.RandomMating(
        numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
        ops=[
            sim.MendelianGenoTransmitter(),
            sim.IdTagger(),
            sim.PedigreeTagger()
        ],
    ),
    gen = 19
)
# we now have the complete pedigree of 20 generations
pop.asPedigree()
# total number of individuals should be 20 * 1000
# how many families do we have?
fam = pop.identifyFamilies()
len(fam)
# but how many families with more than 1 individual?
# The rest of them must be in the initial generation
len([x for x in fam if x > 1])
# let us look backward. allAnc are the ancestors who have offspring in the
# last generation. You can see this is a small number compared the number of
# ancestors.
allAnc = pop.identifyAncestors()
len(allAnc)

