#!/usr/bin/env python

#
# $File: controlledOffGenerator.py $
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit https://github.com/BoPeng/simuPOP
# for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (Bo.Peng@bcm.edu)
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
# the user's guide (https://github.com/BoPeng/simuPOP/manual) for a detailed
# description of this example.
#

import simuPOP as sim
def traj(gen):
    return [0.5 + gen * 0.01]

pop = sim.Population(1000, loci=[10]*2)
# evolve the sim.Population while keeping allele frequency 0.5
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    matingScheme=sim.HomoMating(sim.RandomParentChooser(),
        sim.ControlledOffspringGenerator(loci=5,
            alleles=[0], freqFunc=traj,
            ops = sim.SelfingGenoTransmitter())),
    postOps=[
        sim.Stat(alleleFreq=[5, 15]),
        sim.PyEval(r'"%.2f\t%.2f\n" % (alleleFreq[5][0], alleleFreq[15][0])')
    ],
    gen = 5
)
