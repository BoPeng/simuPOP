#!/usr/bin/env python

#
# $File: demoTerminate.py $
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
import simuPOP.demography as demo

model = demo.MultiStageModel([
    demo.InstantChangeModel(N0=1000, 
        ops=[
            sim.Stat(alleleFreq=sim.ALL_AVAIL, numOfSegSites=sim.ALL_AVAIL),
            # terminate if the average allele frequency of segregating sites
            # are more than 0.1 
            sim.TerminateIf('sum([x[1] for x in alleleFreq.values() if '
                'x[1] != 0])/(1 if numOfSegSites==0 else numOfSegSites) > 0.1')
        ]
    ),
    demo.ExponentialGrowthModel(N0=[0.5, 0.5], r=0.01, NT=[2000, 5000])
    ]
)

pop = sim.Population(size=model.init_size, loci=100)
pop.evolve(
    initOps=sim.InitSex(),
    preOps=sim.SNPMutator(u=0.001, v=0.001),
    matingScheme=sim.RandomMating(subPopSize=model),
    postOps=[
        sim.Stat(alleleFreq=sim.ALL_AVAIL, numOfSegSites=sim.ALL_AVAIL,
            popSize=True, step=50),
        sim.PyEval(r'"%d: %s, %.3f\n" % (gen, subPopSize, sum([x[1] for x '
            'in alleleFreq.values() if x[1] != 0])/(1 if numOfSegSites == 0 '
            'else numOfSegSites))', step=50)
    ],
)


