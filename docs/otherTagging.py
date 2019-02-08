#!/usr/bin/env python

#
# $File: otherTagging.py $
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
pop = sim.Population(1000, loci=[1], infoFields=['aff', 'numOfAff'])
# define virtual subpopulations by affection sim.status
pop.setVirtualSplitter(sim.AffectionSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5]),
    ],
    preOps=[
        # get affection sim.status for parents
        sim.MaPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4]),
        # set 'aff' of parents
        sim.InfoExec('aff = ind.affected()', exposeInd='ind'),
    ],
        # get number of affected parents for each offspring and store in numOfAff
    matingScheme=sim.RandomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.SummaryTagger(mode=sim.SUMMATION, infoFields=['aff', 'numOfAff'])]),
    postOps=[
        # get affection sim.status for offspring
        sim.MaPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4]),
        # calculate mean 'numOfAff' of offspring, for unaffected and affected subpopulations.
        sim.Stat(meanOfInfo='numOfAff', subPops=[(0,0), (0,1)], vars=['meanOfInfo_sp']),
        # print mean number of affected parents for unaffected and affected offspring.
        sim.PyEval(r"'Mean number of affected parents: %.2f (unaff), %.2f (aff)\n' % "
            "(subPop[(0,0)]['meanOfInfo']['numOfAff'], subPop[(0,1)]['meanOfInfo']['numOfAff'])")
    ],
    gen = 5
)


