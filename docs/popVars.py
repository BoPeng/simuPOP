#!/usr/bin/env python

#
# $File: popVars.py $
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
from pprint import pprint
pop = sim.Population(100, loci=2)
sim.initGenotype(pop, freq=[0.3, 0.7])
print(pop.vars())    # No variable now
pop.dvars().myVar = 21
print(pop.vars())
sim.stat(pop, popSize=1, alleleFreq=0)
# pprint prints in a less messy format
pprint(pop.vars())
# print number of allele 1 at locus 0
print(pop.vars()['alleleNum'][0][1])
# use the dvars() function to access dictionary keys as attributes
print(pop.dvars().alleleNum[0][1])
print(pop.dvars().alleleFreq[0])

