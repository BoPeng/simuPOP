#!/usr/bin/env python

#
# $File: PyMlSelector.py $
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

import simuOpt
simuOpt.setOptions(quiet=True, alleleType='mutant')
import simuPOP as sim
import random
pop = sim.Population(size=2000, loci=[10000], infoFields=['fitness'])

class GammaDistributedFitness:
    def __init__(self, alpha, beta):
        self.coefMap = {}
        self.alpha = alpha
        self.beta = beta
     
    def __call__(self, loc, alleles):
        # because s is assigned for each locus, we need to make sure the
        # same s is used for fitness of genotypes 01 (1-s) and 11 (1-2s)
        # at each locus
        if loc in self.coefMap:
            s = self.coefMap[loc]
        else:
            s = random.gammavariate(self.alpha, self.beta)
            self.coefMap[loc] = s
        #
        if 0 in alleles:
            return 1. - s
        else:
            return 1. - 2.*s

pop.evolve(
    initOps=sim.InitSex(),
    preOps=[
        sim.AcgtMutator(rate=[0.00001], model='JC69'),
        sim.PyMlSelector(GammaDistributedFitness(0.23, 0.185),
            output='>>sel.txt'),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(numOfSegSites=sim.ALL_AVAIL, step=50),
        sim.PyEval(r"'Gen: %2d #seg sites: %d\n' % (gen, numOfSegSites)",
            step=50)
    ],
    gen = 201
)
print(''.join(open('sel.txt').readlines()[:5]))

