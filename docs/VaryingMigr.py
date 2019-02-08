#!/usr/bin/env python

#
# $File: VaryingMigr.py $
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

from simuPOP.utils import migrIslandRates
import random

def demo(pop):
  # this function randomly split populations
  numSP = pop.numSubPop()
  if random.random() > 0.3:
      pop.splitSubPop(random.randint(0, numSP-1), [0.5, 0.5])
  return pop.subPopSizes()

def migr(pop):
  numSP = pop.numSubPop()
  sim.migrate(pop, migrIslandRates(0.01, numSP))
  return True

pop = sim.Population(10000, infoFields='migrate_to')
pop.evolve(
    initOps=sim.InitSex(),
    preOps=[
        sim.PyOperator(func=migr),
        sim.Stat(popSize=True),
        sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    matingScheme=sim.RandomMating(subPopSize=demo),
    gen = 5
)

