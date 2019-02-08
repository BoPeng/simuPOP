#!/usr/bin/env python

#
# $File: matingSchemeByRepAndGen.py $
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
simu = sim.Simulator(sim.Population(1000, loci=[10]), rep=3)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    matingScheme=sim.ConditionalMating('rep == 0', 
        # the first replicate use standard random mating
        sim.RandomMating(),
        sim.ConditionalMating('rep == 1 and gen >= 5',
            # the second replicate produces more males for the first 5 generations
            sim.RandomMating(),
            # the last replicate produces more males all the time
            sim.RandomMating(sexMode=(sim.PROB_OF_MALES, 0.7))
            )
        ),
    postOps=[
        sim.Stat(numOfMales=True),
        sim.PyEval("'gen=%d' % gen", reps=0),
        sim.PyEval(r"'\t%d' % numOfMales"),
        sim.PyOutput('\n', reps=-1)
    ],        
    gen=10
)

