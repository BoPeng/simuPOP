#!/usr/bin/env python

#
# $File: migrateVSP.py $
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
pop = sim.Population(size=[1000]*2, infoFields='migrate_to')
pop.setVirtualSplitter(sim.SexSplitter())
pop.evolve(
    # 500 males and 500 females
    initOps=sim.InitSex(sex=[sim.MALE, sim.FEMALE]),
    preOps=[
        sim.Migrator(rate=[
            [0, 0.10],
            [0, 0.05],
            ],
            mode = sim.BY_PROPORTION,
            subPops=[(0, 0), (0, 1)]),
        sim.Stat(popSize=True, numOfMales=True, vars='numOfMales_sp'),
        sim.PyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMales'], subPopSize[0], "
            "subPop[1]['numOfMales'], subPopSize[1])"),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(popSize=True, numOfMales=True, vars='numOfMales_sp'),
        sim.PyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMales'], subPopSize[0], "
            "subPop[1]['numOfMales'], subPopSize[1])"),
    ],
    gen = 2
)   

