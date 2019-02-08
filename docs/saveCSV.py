#!/usr/bin/env python

#
# $File: saveCSV.py $
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
from simuPOP.utils import saveCSV
pop = sim.Population(size=[10], loci=[2, 3],
    lociNames=['r11', 'r12', 'r21', 'r22', 'r23'],
    alleleNames=['A', 'B'], infoFields='age')
sim.initSex(pop)
sim.initInfo(pop, [2, 3, 4], infoFields='age')
sim.initGenotype(pop, freq=[0.4, 0.6])
sim.maPenetrance(pop, loci=0, penetrance=(0.2, 0.2, 0.4))
# no filename so output to standard output
saveCSV(pop, infoFields='age')
# change affection code and how to output genotype
saveCSV(pop, infoFields='age', affectionFormatter={True: 1, False: 2},
    genoFormatter={(0,0):'AA', (0,1):'AB', (1,0):'AB', (1,1):'BB'})
# save to a file
saveCSV(pop, filename='pop.csv', infoFields='age', affectionFormatter={True: 1, False: 2},
    genoFormatter=lambda geno: (geno[0] + 1, geno[1] + 1), sep=' ')
print(open('pop.csv').read())

