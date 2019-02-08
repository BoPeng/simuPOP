#!/usr/bin/env python

#
# $File: importData.py $
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
def importData(filename):
    'Read data from ``filename`` and create a population'
    data = open(filename)
    header = data.readline()
    fields = header.split(',')
    # columns 1, 3, 5, ..., without trailing '_1'
    names = [fields[x].strip()[:-2] for x in range(1, len(fields), 2)]
    popSize = 0
    alleleNames = set()
    for line in data.readlines():
        # get all allele names
        alleleNames |= set([x.strip() for x in line.split(',')[1:]])
        popSize += 1
    # create a population
    alleleNames = list(alleleNames)
    pop = sim.Population(size=popSize, loci=len(names), lociNames=names,
        alleleNames=alleleNames)
    # start from beginning of the file again
    data.seek(0)
    # discard the first line
    data.readline()
    for ind, line in zip(pop.individuals(), data.readlines()):
        fields = [x.strip() for x in line.split(',')]
        sex = sim.MALE if fields[0] == '1' else sim.FEMALE
        ploidy0 = [alleleNames.index(fields[x]) for x in range(1, len(fields), 2)]
        ploidy1 = [alleleNames.index(fields[x]) for x in range(2, len(fields), 2)]
        ind.setGenotype(ploidy0, 0)
        ind.setGenotype(ploidy1, 1)
        ind.setSex(sex)
    # close the file
    data.close()
    return pop

from simuPOP.utils import saveCSV
pop = sim.Population(size=[10], loci=[3, 2], lociNames=['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
    alleleNames=['A', 'B'])
sim.initSex(pop)
sim.initGenotype(pop, freq=[0.5, 0.5])
# output sex but not affection status.
saveCSV(pop, filename='sample.csv', affectionFormatter=None,
    sexFormatter={sim.MALE:1, sim.FEMALE:2})
# have a look at the file
print(open('sample.csv').read())
pop1 = importData('sample.csv')
sim.dump(pop1)

