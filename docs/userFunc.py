#!/usr/bin/env python

#
# $File: userFunc.py $
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
pop = sim.Population(1000, loci=1, infoFields='smoking')
sim.initInfo(pop, lambda:random.randint(0,1), infoFields='smoking')
sim.initGenotype(pop, freq=[0.3, 0.7])

# a penetrance function that depends on smoking
def func(geno, smoking):
    if smoking:
        return (geno[0]+geno[1])*0.4
    else:
        return (geno[0]+geno[1])*0.1

sim.pyPenetrance(pop, loci=0, func=func)
sim.stat(pop, numOfAffected=True)
print(pop.dvars().numOfAffected)


