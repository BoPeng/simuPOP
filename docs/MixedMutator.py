#!/usr/bin/env python

#
# $File: MixedMutator.py $
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
pop = sim.Population(5000, loci=[1, 1])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(genotype=[50, 50])
    ],
    preOps=[
        # the first locus uses a pure stepwise mutation model
        sim.StepwiseMutator(rates=0.001, loci=0),
        # the second locus uses a mixed model
        sim.MixedMutator(rates=0.001, loci=1, mutators=[        
            sim.KAlleleMutator(rates=1, k=100),
            sim.StepwiseMutator(rates=1)
        ], prob=[0.1, 0.9])],
    matingScheme=sim.RandomMating(),
    gen = 20
)
# what alleles are there?
geno0 = []
geno1 = []
for ind in pop.individuals():
    geno0.extend([ind.allele(0, 0), ind.allele(0, 1)])
    geno1.extend([ind.allele(1, 0), ind.allele(1, 1)])

print('Locus 0 has alleles', ', '.join([str(x) for x in set(geno0)]))
print('Locus 1 has alleles', ', '.join([str(x) for x in set(geno1)]))

