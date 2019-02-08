#!/usr/bin/env python

#
# $File: statAssociation.py $
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
from simuPOP.utils import *
from simuPOP.sampling import drawCaseControlSample
def assoTest(pop):
    'Draw case-control sample and apply association tests'
    sample = drawCaseControlSample(pop, cases=500, controls=500)
    sim.stat(sample, association=(0, 2), vars=['Allele_ChiSq_p', 'Geno_ChiSq_p', 'Armitage_p'])
    print('Allele test: %.2e, %.2e, Geno test: %.2e, %.2e, Trend test: %.2e, %.2e' \
        % (sample.dvars().Allele_ChiSq_p[0], sample.dvars().Allele_ChiSq_p[2],
        sample.dvars().Geno_ChiSq_p[0], sample.dvars().Geno_ChiSq_p[2],
        sample.dvars().Armitage_p[0], sample.dvars().Armitage_p[2]))
    return True

pop = sim.Population(size=100000, loci=3)
pop.setVirtualSplitter(sim.ProportionSplitter([0.5, 0.5]))
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(genotype=[0]*3, subPops=[(0,0)]),
        sim.InitGenotype(genotype=[1]*3, subPops=[(0,1)]),
    ],
    matingScheme=sim.RandomMating(ops=sim.Recombinator(loci=[0, 1], rates=[0.01, 0.005])),
    postOps=[
        sim.MaPenetrance(loci=1, penetrance=[0.1, 0.2, 0.4]),
        sim.PyOperator(func=assoTest, step=20),
    ],
    gen = 100
)

