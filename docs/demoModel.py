#!/usr/bin/env python

#
# $File: demoModel.py $
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
from simuPOP.demography import *
model = MultiStageModel([
    InstantChangeModel(T=200, 
        # start with an ancestral population of size 1000
        N0=(1000, 'Ancestral'),
        # change population size at 50 and 60
        G=[50, 60], 
        # change to population size 200 and back to 1000
        NG=[(200, 'bottleneck'), (1000, 'Post-Bottleneck')]),
    ExponentialGrowthModel(
        T=50, 
        # split the population into two subpopulations
        N0=[(400, 'P1'), (600, 'P2')],
        # expand to size 4000 and 5000 respectively
        NT=[4000, 5000])]
    )
#
# model.init_size returns the initial population size
# migrate_to is required for migration
pop = sim.Population(size=model.init_size, loci=1,
    infoFields=model.info_fields)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    matingScheme=sim.RandomMating(subPopSize=model),
    finalOps=
        sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
    gen=model.num_gens
)
# print out population size and frequency
for idx, name in enumerate(pop.subPopNames()):
    print('%s (%d): %.4f' % (name, pop.subPopSize(name), 
        pop.dvars(idx).alleleFreq[0][0]))

# get a visual presentation of the demographic model
model.plot('log/demoModel.png',
    title='A bottleneck + exponential growth demographic model')

