#!/usr/bin/env python

#
# $File: reichEvolve.py $
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


def simulate(model, N0, N1, G0, G1, spec, s, mu, k):
    '''Evolve a sim.Population using given demographic model
    and observe the evolution of its allelic spectrum.
    model: type of demographic model.
    N0, N1, G0, G1: parameters of demographic model.
    spec: initial allelic spectrum, should be a list of allele
        frequencies for each allele.
    s: selection pressure.
    mu: mutation rate.
    k: k for the k-allele model
    '''
    demo_func = demo_model(model, N0, N1, G0, G1)
    pop = sim.Population(size=demo_func(0), loci=1, infoFields='fitness')
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitGenotype(freq=spec, loci=0)
        ],
        matingScheme=sim.RandomMating(subPopSize=demo_func),
        postOps=[
            sim.KAlleleMutator(k=k, rates=mu),
            sim.MaSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=[0], step=100),
            sim.PyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
                step=100),
        ],
        gen = G0 + G1
    )

simulate('instant', 1000, 10000, 500, 500, [0.9]+[0.02]*5, 0.01, 1e-4, 200)

