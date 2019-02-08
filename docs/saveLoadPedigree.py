#!/usr/bin/env python

#
# $File: saveLoadPedigree.py $
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
pop = sim.Population(4, loci=1, infoFields=['ind_id', 'father_id', 'mother_id'],
    ancGen=-1)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
        sim.InitGenotype(freq=[0.5, 0.5]),
        sim.PedigreeTagger(output='>>pedigree.ped', outputLoci=0)
    ],
    matingScheme=sim.RandomMating(
        ops=[
            sim.MendelianGenoTransmitter(),
            sim.IdTagger(),
            sim.PedigreeTagger(output='>>pedigree.ped', outputLoci=0)
        ],
    ),
    gen = 2
)
#
print(open('pedigree.ped').read())
pop.asPedigree()
pop.save('pedigree1.ped', loci=0)
print(open('pedigree1.ped').read())
# 
ped = sim.loadPedigree('pedigree1.ped')
sim.dump(ped, ancGens=range(3))

