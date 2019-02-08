#!/usr/bin/env python

#
# $File: MapPenetrance.py $
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
pop = sim.Population(size=2000, loci=2)
sim.initGenotype(pop, freq=[.2, .8])
sim.mapPenetrance(pop, loci=0,
    penetrance={(0,0):0, (0,1):.2, (1,1):.3})
sim.stat(pop, genoFreq=0, numOfAffected=1, vars='genoNum')
# number of affected individuals
pop.dvars().numOfAffected
# which should be roughly (#01 + #10) * 0.2 + #11 * 0.3
(pop.dvars().genoNum[0][(0,1)] + pop.dvars().genoNum[0][(1,0)]) * 0.2 \
+ pop.dvars().genoNum[0][(1,1)] * 0.3

