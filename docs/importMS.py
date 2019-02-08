#!/usr/bin/env python

#
# $File: importMS.py $
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
from simuPOP.utils import importPopulation, export
pop = sim.Population([20,20], loci=[10, 10])
# simulate a population but mutate only a subset of loci
pop.evolve(
    preOps=[
        sim.InitSex(),
        sim.SNPMutator(u=0.1, v=0.01, loci=range(5, 17))
    ],
    matingScheme=sim.RandomMating(),
    gen=100
)
# export first chromosome, all individuals
export(pop, format='ms', output='ms.txt')
# export first chromosome, subpops as replicates
export(pop, format='ms', output='ms_subPop.txt', splitBy='subPop')
# export all chromosomes, but limit to all males in subPop 1
pop.setVirtualSplitter(sim.SexSplitter())
export(pop, format='ms', output='ms_chrom.txt', splitBy='chrom', subPops=[(1,0)])
# 
print(open('ms_chrom.txt').read())
# import as haploid sequence
pop = importPopulation(format='ms', filename='ms.txt')
# import as diploid 
pop = importPopulation(format='ms', filename='ms.txt', ploidy=2)
# import as a single chromosome
pop = importPopulation(format='ms', filename='ms_subPop.txt', mergeBy='subPop')

