#!/usr/bin/env python

#
# $File: InfoEval.py $
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
pop = sim.Population(20, loci=1, infoFields='a')
pop.setVirtualSplitter(sim.InfoSplitter('a', cutoff=[3]))
sim.initGenotype(pop, freq=[0.2, 0.8])
pop.setIndInfo([random.randint(2, 5) for x in range(20)], 'a')
sim.infoEval(pop, 'a', subPops=[(0, 0)]);print(' ')
sim.infoEval(pop, 'ind.allele(0, 0)', exposeInd='ind');print(' ')
# use sim.population variables
pop.dvars().b = 5
sim.infoEval(pop, '"%d " % (a+b)');print(' ')

