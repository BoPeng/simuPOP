#!/usr/bin/env python

#
# $File: applicableGen.py $
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
pop = sim.Population(1000, loci=[20])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.8, 0.2])
    ],
    preOps=[
        sim.PyEval(r"'At the beginning of gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1])
    ],
    matingScheme = sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=0, begin=80, step=10),
        sim.PyEval(r"'At the end of gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        sim.PyEval(r"'At the end of gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1])
    ],
    finalOps=sim.SavePopulation(output='sample.pop'),
    gen=100
)

