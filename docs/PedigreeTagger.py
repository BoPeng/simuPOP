#!/usr/bin/env python

#
# $File: PedigreeTagger.py $
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
pop = sim.Population(100, infoFields=['ind_id', 'father_id', 'mother_id'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
        sim.PedigreeTagger(output='>>pedigree.txt'),
    ],
    matingScheme=sim.RandomMating(ops=[
        sim.IdTagger(),
        sim.PedigreeTagger(output='>>pedigree.txt'),
        sim.MendelianGenoTransmitter()]
    ),
    gen = 100
)
ped = open('pedigree.txt')
lines = ped.readlines()
ped.close()
# first few lines, saved by the first PedigreeTagger
print(''.join(lines[:3]))
# last several lines, saved by the second PedigreeTagger
print(''.join(lines[-3:]))
# load this file
ped = sim.loadPedigree('pedigree.txt')
# should have 100 ancestral generations (plus one present generation)
ped.ancestralGens()

