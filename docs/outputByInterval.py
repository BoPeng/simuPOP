#!/usr/bin/env python

#
# $File: outputByInterval.py $
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
import time
pop = sim.Population(1000, loci=10)
pop.dvars().init_time = time.time()
pop.dvars().last_time = time.time()
exec('import time', pop.vars(), pop.vars())
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.IfElse('time.time() - last_time > 5', [
            sim.PyEval(r'"Gen: %d\n" % gen'),
            sim.PyExec('last_time = time.time()')
            ]),
        sim.TerminateIf('time.time() - init_time > 20')
    ]
)
        

