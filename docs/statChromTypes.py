#!/usr/bin/env python

#
# $File: statChromTypes.py $
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
    chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y, sim.MITOCHONDRIAL])
pop.setVirtualSplitter(sim.SexSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(haplotypes=[ [0, 1, 2, 0, 1]*4, [2, 1, 0, 2, 3]*4 ],
            prop=[0.4, 0.6]),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            sim.MendelianGenoTransmitter(),
            sim.MitochondrialGenoTransmitter()]),
    preOps=[
        sim.Stat(neutrality=range(5)),
        sim.Stat(neutrality=range(5, 10), suffix='_X'),
        sim.Stat(neutrality=range(10, 15), suffix='_Y'),
        sim.Stat(neutrality=range(15, 20), suffix='_mt'),
        sim.PyEval(r'"%.3f %.3f %.3f %.3f\n" % (Pi, Pi_X, Pi_Y, Pi_mt)'),
    ],
    gen = 2
)

