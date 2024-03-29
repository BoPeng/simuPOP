#!/usr/bin/env python

#
# $File: PyTagger.py $
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
import random
def randomMove(x, y):
    '''Pass parental information fields to offspring'''
    # shift right with high concentration of alleles...
    off_x = random.normalvariate((x[0]+x[1])/2., 0.1)
    off_y = random.normalvariate((y[0]+y[1])/2., 0.1)
    return off_x, off_y

pop = sim.Population(1000, loci=[1], infoFields=['x', 'y'])
pop.setVirtualSplitter(sim.GenotypeSplitter(loci=0, alleles=[[0, 0], [0,1], [1, 1]]))
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5]),
        sim.InitInfo(random.random, infoFields=['x', 'y'])
    ],
    matingScheme=sim.RandomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.PyTagger(func=randomMove),
    ]),
    postOps=[
        sim.Stat(minOfInfo='x', maxOfInfo='x'),
        sim.PyEval(r"'Range of x: %.2f, %.2f\n' % (minOfInfo['x'], maxOfInfo['x'])")
    ],
    gen = 5
)
