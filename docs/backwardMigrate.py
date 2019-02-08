#!/usr/bin/env python

#
# $File: backwardMigrate.py $
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
sim.turnOnDebug('DBG_MIGRATOR')
pop = sim.Population(size=[10000, 5000, 8000], infoFields=['migrate_to', 'migrate_from'])
def originOfInds(pop):
    print('Observed backward migration matrix at generation {}'.format(pop.dvars().gen))
    for sp in range(pop.numSubPop()): 
        # get source subpop for all individuals in subpopulation i
        origins = pop.indInfo('migrate_from', sp)
        spSize = pop.subPopSize(sp)
        B_sp = [origins.count(j) * 1.0 /spSize for j in range(pop.numSubPop())]
        print('    ' + ', '.join(['{:.3f}'.format(x) for x in B_sp]))
    return True

pop.evolve(
    initOps=sim.InitSex(),
    preOps=
        # mark the source subpopulation of each individual
        [sim.InitInfo(i, subPops=i, infoFields='migrate_from') for i in range(3)] + [
        # perform migration
        sim.BackwardMigrator(rate=[
            [0, 0.04, 0.02],
            [0.05, 0, 0.02],
            [0.02, 0.01, 0]
        ]),
        # calculate and print observed backward migration matrix 
        sim.PyOperator(func=originOfInds),
        # calculate population size
        sim.Stat(popSize=True),
        # and print it
        sim.PyEval(r'"Pop size after migration: {}\n".format(", ".join([str(x) for x in subPopSize]))'),
        ], 
    matingScheme=sim.RandomMating(),
    gen = 5
)        

