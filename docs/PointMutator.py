#!/usr/bin/env python

#
# $File: PointMutator.py $
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
pop = sim.Population(1000, loci=1, infoFields='fitness')
pop.evolve(
    initOps=sim.PyOutput('Introducing alleles at generation'),
    preOps=sim.MaSelector(loci=0, wildtype=0, fitness=[1, 1.05, 1.1]),
    matingScheme=sim.RandomSelection(),
    postOps=[
        sim.Stat(alleleFreq=0),
        sim.IfElse('alleleNum[0][1] == 0', ifOps=[
            sim.PyEval(r"' %d' % gen"),
            sim.PointMutator(inds=0, loci=0, allele=1),
        ]),
        sim.IfElse('alleleFreq[0][1] > 0.05', ifOps=[
            sim.PyEval(r"'.\nTerminate at generation %d at allele freq %.3f.\n'" +
                " % (gen, alleleFreq[0][1])"),
            sim.TerminateIf('True'),
        ])
    ],
)

