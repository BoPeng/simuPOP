#!/usr/bin/env python

#
# $File: peneSelector.py $
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
pop = sim.Population(size=2000, loci=1, infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[.5, .5])
    ],
    preOps=[
        sim.MaPenetrance(loci=0, penetrance=[0.01, 0.1, 0.2]),
        sim.Stat(numOfAffected=True, step=25, vars='propOfAffected'),
        sim.PyEval(r"'Percent of affected: %.3f\t' % propOfAffected", step=50),
        sim.InfoExec('fitness = not ind.affected()', exposeInd='ind')
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=0),
        sim.PyEval(r"'%.4f\n' % alleleFreq[0][1]", step=50)
    ],
    gen=151
)
