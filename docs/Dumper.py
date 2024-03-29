#!/usr/bin/env python

#
# $File: Dumper.py $
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
pop = sim.Population(size=[10, 10], loci=[20, 30], infoFields='gen',
    ancGen=-1)
sim.initSex(pop)
pop.setVirtualSplitter(sim.SexSplitter())
pop1 = pop.clone()
sim.initGenotype(pop, freq=[0]*20 + [0.1]*10)
pop.setIndInfo(1, 'gen')
sim.initGenotype(pop1, freq=[0]*50 + [0.1]*10)
pop1.setIndInfo(2, 'gen')
pop.push(pop1)
sim.dump(pop, width=3, loci=[5, 6, 30], subPops=([0, 0], [1, 1]),
    max=10, structure=False)
# list all male individuals in all subpopulations
sim.dump(pop, width=3, loci=[5, 6, 30], subPops=[(sim.ALL_AVAIL, 0)],
    max=10, structure=False)
