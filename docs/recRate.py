#!/usr/bin/env python

#
# $File: recRate.py $
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
simu = sim.Simulator(sim.Population(size=[1000], loci=[100]),
    rep=2)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(genotype=[0]*100 + [1]*100)
    ],
    matingScheme=sim.RandomMating(ops = [
        sim.Recombinator(rates=0.01, reps=0),
        sim.Recombinator(rates=[0.01]*10, loci=range(50, 60), reps=1),
    ]),
    postOps=[
        sim.Stat(LD=[[40, 55], [60, 70]]),
        sim.PyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[60][70])'),
        sim.PyOutput('\n', reps=-1)
    ],
    gen = 5
)
