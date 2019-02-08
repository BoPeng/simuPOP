#!/usr/bin/env python

#
# $File: virtualSubPop.py $
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
pop = sim.Population(10, loci=[2, 3], infoFields='Sex')
sim.initSex(pop)
pop.setVirtualSplitter(sim.SexSplitter())
# initialize male and females with different genotypes. 
sim.initGenotype(pop, genotype=[0]*5, subPops=[(0, 0)])
sim.initGenotype(pop, genotype=[1]*5, subPops=[(0, 1)])
# set Sex information field to 0 for all males, and 1 for all females
pop.setIndInfo([sim.MALE], 'Sex', [0, 0])
pop.setIndInfo([sim.FEMALE], 'Sex', [0, 1])
# Print individual genotypes, followed by values at information field Sex
sim.dump(pop, structure=False)

