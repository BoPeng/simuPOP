#!/usr/bin/env python

#
# $File: outputFunc.py $
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
import logging
# logging to a file simulation.log, with detailed debug information
logging.basicConfig(
    filename='simulation.log',
    level=logging.DEBUG,
    format='%(levelname)s: %(message)s',
    filemode='w'
)
formatter = logging.Formatter('%(message)s')
logger = logging.getLogger('')
pop = sim.Population(size=1000, loci=2)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(genotype=[1, 2, 2, 1])
    ],
    matingScheme = sim.RandomMating(ops=sim.Recombinator(rates=0.01)),
    postOps=[
        sim.Stat(LD=[0, 1]),
        sim.PyEval(r"'LD: %d, %.2f' % (gen, LD[0][1])", step=20,
            output=logger.info),   # send LD to console and a logfile
        sim.PyEval(r"'R2: %d, %.2f' % (gen, R2[0][1])", step=20,
            output=logger.debug),  # send R2 only to a logfile
    ],
    gen=100
)
print(open('simulation.log').read())

