#!/usr/bin/env python

#
# $File: geneticContribution.py $
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

import simuOpt
simuOpt.setOptions(alleleType='lineage', quiet=True)
import simuPOP as sim
pop = sim.Population(1000, loci=[10]*4)

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.25]*4),
        sim.InitLineage(range(1000), mode=sim.PER_INDIVIDUAL),
    ],
    matingScheme=sim.RandomMating(ops=sim.Recombinator(rates=0.001)),
    gen = 100
)
# average number of 'contributors'
num_contributors = [len(set(ind.lineage())) for ind in pop.individuals()]
print('Average number of contributors is %.2f' % (sum(num_contributors) / float(pop.popSize())))
# percent of genetic information from each ancestor (baseline is 1/1000)
lineage = pop.lineage()
lin_perc = [lineage.count(x)/float(len(lineage)) for x in range(1000)]
# how many of ancestors do not have any allele left?
print('Number of ancestors with no allele left: %d' % lin_perc.count(0.))
# top five contributors
lin_perc.sort()
lin_perc.reverse()
print('Top contributors (started with 0.001): %.5f %.5f %.5f' % (lin_perc[0], lin_perc[1], lin_perc[2]))

