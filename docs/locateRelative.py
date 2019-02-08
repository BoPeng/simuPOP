#!/usr/bin/env python

#
# $File: locateRelative.py $
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
pop = sim.Population(1000, ancGen=2, infoFields=['ind_id', 'father_id', 'mother_id'])
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
    gen = 5
)
ped = sim.Pedigree(pop)
offFields = ['off%d' % x for x in range(4)]
grandOffFields = ['grandOff%d' % x for x in range(5)]
ped.addInfoFields(['spouse'] + offFields + grandOffFields)
# only look spouse for fathers...
ped.locateRelatives(sim.OUTBRED_SPOUSE, ['spouse'], sex=sim.FEMALE_ONLY)
ped.locateRelatives(sim.COMMON_OFFSPRING, ['spouse'] + offFields)
# trace offspring of offspring
ped.traceRelatives([offFields, offFields], resultFields=grandOffFields)
# 
IDs = ped.individualsWithRelatives(grandOffFields)
# check on ID.
grandFather = IDs[0]
grandMother = ped.indByID(grandFather).spouse
# some ID might be invalid.
children = [ped.indByID(grandFather).info(x) for x in offFields]
childrenSpouse = [ped.indByID(x).spouse for x in children if x >= 1]
childrenParents = [ped.indByID(x).father_id for x in children if x >= 1] \
    + [ped.indByID(x).mother_id for x in children if x >= 1]
grandChildren = [ped.indByID(grandFather).info(x) for x in grandOffFields]
grandChildrenParents = [ped.indByID(x).father_id for x in grandChildren if x >= 1] \
    + [ped.indByID(x).mother_id for x in grandChildren if x >= 1]

def idString(IDs):
    uniqueIDs = list(set(IDs))
    uniqueIDs.sort()
    return ', '.join(['%d' % x for x in uniqueIDs if x >= 1])

print('''GrandParents: %d, %d
Children: %s
Spouses of children: %s
Parents of children: %s
GrandChildren: %s
Parents of grandChildren: %s ''' % \
(grandFather, grandMother, idString(children), idString(childrenSpouse),
    idString(childrenParents), idString(grandChildren), idString(grandChildrenParents)))

# let us look at the structure of this complete pedigree using another method
famSz = ped.identifyFamilies()
# it is amazing that there is a huge family that connects almost everyone
len(famSz), max(famSz)
# if we only look at the last two generations, things are much better
ped.addInfoFields('ped_id')
famSz = ped.identifyFamilies(pedField='ped_id', ancGens=[0,1])
len(famSz), max(famSz)

