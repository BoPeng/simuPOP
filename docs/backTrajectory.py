#!/usr/bin/env python

#
# $File: backTrajectory.py $
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
from simuPOP.utils import Trajectory, simulateBackwardTrajectory
from math import exp
def Nt(gen):
    'An exponential sim.Population growth demographic model.'
    return int((5000) * exp(.00115 * gen))

def fitness(gen, sp):
    'Constant positive selection pressure.'
    return [1, 1.01, 1.02]

# simulate a trajectory backward in time, from generation 1000
traj = simulateBackwardTrajectory(N=Nt, fitness=fitness, nLoci=2,
     endGen=1000, endFreq=[0.1, 0.2])
# matplotlib syntax
#traj.plot('log/backTrajectory.png', set_ylim_top=0.3, set_ylim_bottom=0,
#        plot_c_loc=['r', 'b'], set_title_label='Simulated Trajectory (backward-time)')

print('Trajectory simulated with length %s ' % len(traj.traj))
pop = sim.Population(size=Nt(0), loci=[1]*2)
# save Trajectory function in the sim.population's local namespace
# so that the sim.PyEval operator can access it.
pop.dvars().traj = traj.func()
pop.evolve(
    initOps=[sim.InitSex()],
    preOps=traj.mutators(loci=[0, 1]),
    matingScheme=sim.ControlledRandomMating(loci=[0, 1], alleles=[1, 1],
        subPopSize=Nt, freqFunc=traj.func()),
    postOps=[
        sim.Stat(alleleFreq=[0, 1], begin=500, step=100),
        sim.PyEval(r"'%4d: %.3f (exp: %.3f), %.3f (exp: %.3f)\n' % (gen, alleleFreq[0][1],"
            "traj(gen)[0], alleleFreq[1][1], traj(gen)[1])",
            begin=500, step=100)
    ],
    gen=1001  # evolve 1001 generations to reach the end of generation 1000
)

