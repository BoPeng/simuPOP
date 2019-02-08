#!/usr/bin/env python

#
# $File: PyPenetrance.py $
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
pop = sim.Population(size=2000, loci=[1]*2, infoFields=['p', 'smoking'])
pop.setVirtualSplitter(sim.InfoSplitter(field='smoking', values=[0,1]))
# the second parameter gen can be used for varying selection pressure
def penet(geno, smoking):
    #     BB     Bb      bb
    # AA  0.01   0.01    0.01
    # Aa  0.01   0.03    0.03
    # aa  0.01   0.03    0.05
    #
    # geno is (A1 A2 B1 B2)
    if geno[0] + geno[1] == 1 and geno[2] + geno[3] != 0:
        v = 0.03   # case of AaBb
    elif geno[0] + geno[1] == 2 and geno[2] + geno[3] == 1:
        v = 0.03   # case of aaBb
    elif geno[0] + geno[1] ==2 and geno[2] + geno[3] == 2:
        v = 0.05   # case of aabb
    else:                
        v = 0.01   # other cases
    if smoking:
        return v * 2
    else:
        return v

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[.5, .5]),
        sim.PyOutput('Calculate prevalence in smoker and non-smokers\n'),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        # set smoking status randomly
        sim.InitInfo(lambda : random.randint(0,1), infoFields='smoking'),
        # assign affection status
        sim.PyPenetrance(loci=[0, 1], func=penet),
        sim.Stat(numOfAffected=True, subPops=[(0, sim.ALL_AVAIL)], 
            vars='propOfAffected_sp', step=20),
        sim.PyEval(r"'Non-smoker: %.2f%%\tSmoker: %.2f%%\n' % "
            "(subPop[(0,0)]['propOfAffected']*100, subPop[(0,1)]['propOfAffected']*100)",
            step=20)
    ],
    gen = 50
)


