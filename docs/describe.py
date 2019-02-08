#!/usr/bin/env python

#
# $File: describe.py $
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

def outputstat(pop):
    'Calculate and output statistics, ignored'
    return True

# describe this evolutionary process
print(sim.describeEvolProcess(
    initOps=[
        sim.InitSex(),
        sim.InitInfo(lambda: random.randint(0, 75), infoFields='age'),
        sim.InitGenotype(freq=[0.5, 0.5]),
        sim.IdTagger(),
        sim.PyOutput('Prevalence of disease in each age group:\n'),
    ],
    preOps=sim.InfoExec('age += 1'),
    matingScheme=sim.HeteroMating([
        sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
        sim.RandomMating(ops=[
            sim.IdTagger(),
            sim.Recombinator(intensity=1e-4)
        ], subPops=[(0,1)]),
    ]),
    postOps=[
        sim.MaPenetrance(loci=0, penetrance=[0.01, 0.1, 0.3]),
        sim.PyOperator(func=outputstat)
    ],
    gen = 100,
    numRep = 3
))     

