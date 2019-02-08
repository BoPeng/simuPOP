#!/usr/bin/env python

#
# $File: debug.py $
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
# redirect system stderr
import sys
debugOutput = open('debug.txt', 'w')
old_stderr = sys.stderr
sys.stderr = debugOutput
# start simulation
simu = sim.Simulator(sim.Population(100, loci=1), rep=5)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.1, 0.9])
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=0),
        sim.IfElse('alleleNum[0][0] == 0',
            ifOps=[
                # the is None part makes the function return True
                sim.PyOperator(lambda : sim.turnOnDebug("DBG_MUTATOR") is None),
                sim.PointMutator(loci=0, allele=0, inds=0),
            ],
            elseOps=sim.PyOperator(lambda : sim.turnOffDebug("DBG_MUTATOR") is None)),
    ],
    gen = 100
)
# replace standard stdandard error
sys.stderr = old_stderr
debugOutput.close()
print(''.join(open('debug.txt').readlines()[:5]))

