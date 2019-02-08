#!/usr/bin/env python

#
# $File: infiniteSites.py $
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
simuOpt.setOptions(alleleType='long')
import simuPOP as sim

def infSitesMutate(pop, param):
    '''Apply an infinite mutation model'''
    (startPos, endPos, rate) = param
    # for each individual
    for ind in pop.individuals():
        # for each homologous copy of chromosomes
        for p in range(2):
            # using a geometric distribution to determine
            # the first mutation location
            loc = sim.getRNG().randGeometric(rate)
            # if a mutation happens, record the mutated location
            if startPos + loc < endPos:
                try:
                    # find the first non-zero location
                    idx = ind.genotype(p).index(0)
                    # record mutation here
                    ind.setAllele(startPos + loc, idx, ploidy=p)
                except:
                    raise
                    print('Warning: more than %d mutations have accumulated' % pop.totNumLoci())
                    pass
    return True

pop = sim.Population(size=[2000], loci=[100])
pop.evolve(
    initOps=sim.InitSex(),
    preOps=[
        # mutate in a 10Mb region at rate 1e-8
        sim.PyOperator(func=infSitesMutate, param=(1, 10000000, 1e-8)),
    ],
    matingScheme=sim.RandomMating(),
    gen = 100
)
# now, we get a sim.Population. Let us have a look at the 'alleles'.
# print the first five mutation locations
print(pop.individual(0).genotype()[:5])
# how many alleles are there (does not count 0)?
print(len(set(pop.genotype())) - 1)
# Allele count a simple count of alleles.
cnt = {}
for allele in pop.genotype():
    if allele == 0:
        continue
    if allele in cnt:
        cnt[allele] += 1
    else:
        cnt[allele] = 1

# highest allele frequency?
print(max(cnt.values()) *0.5 / pop.popSize())

