#!/usr/bin/env python

#
# $File: ageOfMutants.py $
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

import simuOpt
simuOpt.setOptions(alleleType='lineage', quiet=True)
import simuPOP as sim
pop = sim.Population(size=10000, loci=[10]*10, infoFields='ind_id')
# just to make sure IDs starts from 1
sim.IdTagger().reset(1)
pop.evolve(
    initOps = [
        sim.InitSex(),
        sim.InitGenotype(freq=[0.2, 0.3, 0.4, 0.1]),
        sim.IdTagger(),
        sim.InitLineage(mode=sim.FROM_INFO),
    ],
    # an extremely high mutation rate, just for demonstration
    preOps = sim.AcgtMutator(rate=0.01, model='JC69'),
    matingScheme=sim.RandomMating(
        ops=[
            sim.IdTagger(),
            sim.MendelianGenoTransmitter(),
        ]
    ),
    gen = 10
)
lin = pop.lineage()
# Number of alleles from each generation
for gen in range(10):
    id_start = gen*10000 + 1
    id_end = (gen+1)*10000
    num_mut = len([x for x in lin if x >= id_start and x <= id_end])
    print('Gen %d: %5.2f %%' % (gen, num_mut / (2*10000*100.) * 100))


