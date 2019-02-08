#!/usr/bin/env python

#
# $File: sampleNuclearFamily.py $
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
from simuPOP.sampling import drawNuclearFamilySample
pop = sim.loadPopulation('log/pedigree.pop')
sample = drawNuclearFamilySample(pop, families=2, numOffspring=(2,4),
    affectedParents=(1,2), affectedOffspring=(1, 3))
# try to separate two families?
sample.asPedigree()
#= sim.Pedigree(sample, loci=sim.ALL_AVAIL, infoFields=sim.ALL_AVAIL)
sample.addInfoFields('ped_id')
# return size of families
sz = sample.identifyFamilies(pedField='ped_id')
print(sz)
ped1 = sample.extractIndividuals(IDs=0, idField='ped_id')
# print the ID of all individuals in the first pedigree
print([ind.ind_id for ind in ped1.allIndividuals()])

