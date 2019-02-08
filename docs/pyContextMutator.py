#!/usr/bin/env python

#
# $File: pyContextMutator.py $
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
import random
pop = sim.Population(5000, loci=[3, 3])
def contextMut(allele, context):
    if context == [0, 0]:
        if allele == 0 and random.random() < 0.1:
            return 1
    elif context == [1, 1]:
        if allele == 0:
            return 1
    # do not mutate
    return allele

pop.evolve(
    # initialize locus by 0, 0, 0, 1, 0, 1
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(genotype=[1, 1], loci=[3, 5])
    ],
    preOps=[
        sim.PyMutator(func=contextMut, context=1,
            loci=[1, 4],  rates=0.01
        ),
        #sim.SNPMutator(u=0.01, v= 0.01, loci=[1, 4]),
        sim.Stat(alleleFreq=[1, 4], step=5),
        sim.PyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    matingScheme=sim.RandomMating(),
    gen = 20
)

