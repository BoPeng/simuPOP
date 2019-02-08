#!/usr/bin/env python

#
# $File: forwardTrajectory.py $
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
simuOpt.setOptions(quiet=True)
import simuPOP as sim
from simuPOP.utils import Trajectory, simulateForwardTrajectory

traj = simulateForwardTrajectory(N=[2000, 4000], fitness=[1, 0.99, 0.98],
    beginGen=0, endGen=100, beginFreq=[0.2, 0.3],
    endFreq=[[0.1, 0.11], [0.2, 0.21]])
# 
#traj.plot('log/forwardTrajectory.png', set_ylim_top=0.5,
#    plot_c_sp=['r', 'b'], set_title_label='Simulated Trajectory (forward-time)')
pop = sim.Population(size=[2000, 4000], loci=10, infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.8, 0.2], subPops=0),
        sim.InitGenotype(freq=[0.7, 0.3], subPops=1),
        sim.PyOutput('Sp0: loc2\tloc5\tSp1: loc2\tloc5\n'),
    ],
    matingScheme=sim.ControlledRandomMating(
        ops=[sim.Recombinator(rates=0.01)],
        loci=5, alleles=1, freqFunc=traj.func()),
    postOps=[
        sim.Stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=20),
        sim.PyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
            "subPop[0]['alleleFreq'][5][1], subPop[1]['alleleFreq'][2][1],"
            "subPop[1]['alleleFreq'][5][1])", step=20)
    ],
    gen = 101
)

