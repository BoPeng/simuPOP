#!/usr/bin/env python

#
# $File: simuGen.py $
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
simu = sim.Simulator(sim.Population(50, loci=[10], ploidy=1),
    rep=3)
simu.evolve(gen = 5)
simu.dvars(0).gen
simu.evolve(
    initOps=[sim.InitGenotype(freq=[0.5, 0.5])],
    matingScheme=sim.RandomSelection(),
    postOps=[
        sim.Stat(alleleFreq=5),
        sim.IfElse('alleleNum[5][0] == 0',
            sim.PyEval(r"'Allele 0 is lost in rep %d at gen %d\n' % (rep, gen)")),
        sim.IfElse('alleleNum[5][0] == 50',
            sim.PyEval(r"'Allele 0 is fixed in rep %d at gen %d\n' % (rep, gen)")),
        sim.TerminateIf('len(alleleNum[5]) == 1'),
    ],
)
[simu.dvars(x).gen for x in range(3)]

