#!/usr/bin/env python

#
# $File: InitGenotype.py $
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
pop = sim.Population(size=[2000, 3000], loci=[5, 7])
# by allele frequency
def printFreq(pop, loci):
    sim.stat(pop, alleleFreq=loci)
    print(', '.join(['{:.3f}'.format(pop.dvars().alleleFreq[x][0]) for x in loci]))

sim.initGenotype(pop, freq=[.4, .6])
sim.dump(pop, max=6, structure=False)
printFreq(pop, range(5))
# by proportion
sim.initGenotype(pop, prop=[0.4, 0.6])
printFreq(pop, range(5))
# by haplotype frequency
sim.initGenotype(pop, freq=[.4, .6], haplotypes=[[1, 1, 0, 1], [0, 0, 1]])
sim.dump(pop, max=6, structure=False)
printFreq(pop, range(5))
# by haplotype proportion
sim.initGenotype(pop, prop=[0.4, 0.6], haplotypes=[[1, 1, 0], [0, 0, 1, 1]])
printFreq(pop, range(5))
# by genotype
pop = sim.Population(size=[2, 3], loci=[5, 7])
sim.initGenotype(pop, genotype=[1]*5 + [2]*7 + [3]*5 +[4]*7)
sim.dump(pop, structure=False)
# 
# use virtual subpopulation
pop = sim.Population(size=[2000, 3000], loci=[5, 7])
pop.setVirtualSplitter(sim.SexSplitter())
sim.initSex(pop)
sim.initGenotype(pop, genotype=range(10), loci=range(5))
# initialize all males
sim.initGenotype(pop, genotype=[2]*7, loci=range(5, 12),
    subPops=[(0, 0), (1, 0)])
# assign genotype by proportions
pop.setVirtualSplitter(sim.ProportionSplitter([0.4, 0.6]))
sim.initGenotype(pop, freq=[0.2, 0.8], subPops=[(0,0)])
sim.initGenotype(pop, freq=[0.5, 0.5], subPops=[(0,1)])
#
# initialize by random allele frequency
import random
sim.initGenotype(pop, freq=lambda : random.random())
printFreq(pop, range(5))
# initialize with loci specific frequency. here
# lambda loc: 0.01*loc is equivalent to 
# lambda loc: [0.01*loc, 1-0.01*loc]
sim.initGenotype(pop,
    freq=lambda loc: 0.01*loc)
printFreq(pop, range(5))
# initialize with VSP-specific frequency
sim.initGenotype(pop,
    freq=lambda vsp: [[0.2, 0.8], [0.5, 0.5]][vsp[1]],
    subPops=[(0, 0), (0, 1)])


