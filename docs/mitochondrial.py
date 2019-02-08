#!/usr/bin/env python

#
# $File: mitochondrial.py $
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
pop = sim.Population(1000, loci=[5]*4,
    # one autosome, two sex chromosomes, and one mitochondrial chromosomes
    chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y, sim.MITOCHONDRIAL],
    infoFields=['fitness'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.25]*4)
    ],
    preOps=[
        sim.MapSelector(loci=17, fitness={(0,): 1, (1,): 1, (2,): 1, (3,): 0.4})
    ],
    matingScheme=sim.RandomMating(ops= [
        sim.Recombinator(rates=0.1),
        sim.MitochondrialGenoTransmitter(),
    ]),
    postOps=[
        sim.Stat(alleleFreq=17, step=10),
        sim.PyEval(r'"%.2f %.2f %.2f %.2f\n" % (alleleNum[17][0],'
            'alleleNum[17][1], alleleNum[17][2], alleleNum[17][3])', step=10),
    ],
    gen = 100
)

