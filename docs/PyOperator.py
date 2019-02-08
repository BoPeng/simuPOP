#!/usr/bin/env python

#
# $File: PyOperator.py $
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
def dynaMutator(pop, param):
    '''This mutator mutates commom loci with low mutation rate and rare
    loci with high mutation rate, as an attempt to raise allele frequency
    of rare loci to an higher level.'''
    # unpack parameter
    (cutoff, mu1, mu2) = param;
    sim.stat(pop, alleleFreq=range(pop.totNumLoci()))
    for i in range(pop.totNumLoci()):
        # Get the frequency of allele 1 (disease allele)
        if pop.dvars().alleleFreq[i][1] < cutoff:
            sim.kAlleleMutate(pop, k=2, rates=mu1, loci=[i])
        else:
            sim.kAlleleMutate(pop, k=2, rates=mu2, loci=[i])
    return True

pop = sim.Population(size=10000, loci=[2, 3])
pop.evolve(
    initOps=[ 
        sim.InitSex(),
        sim.InitGenotype(freq=[.99, .01], loci=[0, 2, 4]),
        sim.InitGenotype(freq=[.8, .2], loci=[1, 3])
    ],
    preOps=sim.PyOperator(func=dynaMutator, param=(.2, 1e-2, 1e-5)),
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=range(5), step=10),
        sim.PyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)                

