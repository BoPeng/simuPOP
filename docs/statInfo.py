#!/usr/bin/env python

#
# $File: statInfo.py $
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
import random
pop = sim.Population([500], infoFields='anc')
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sim.InfoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
#
pop.evolve(
    initOps=[
        sim.InitSex(),
        # anc is 0 or 1
        sim.InitInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    matingScheme=sim.RandomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.InheritTagger(mode=sim.MEAN, infoFields='anc')
    ]),
    postOps=[
        sim.Stat(popSize=True, meanOfInfo='anc', varOfInfo='anc',
            subPops=[(0, sim.ALL_AVAIL)]),
        sim.PyEval(r"'Anc: %.2f (%.2f), #inds: %s\n' %" + \
            "(meanOfInfo['anc'], varOfInfo['anc'], " + \
            "', '.join(['%4d' % x for x in subPopSize]))")
    ],
    gen = 5,
)

