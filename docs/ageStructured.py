#!/usr/bin/env python

#
# $File: ageStructured.py $
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
N = 10000
pop = sim.Population(N, loci=1, infoFields=['age', 'ind_id', 'father_id', 'mother_id'])
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[20, 50, 75]))
def demoModel(gen, pop):
    '''A demographic model that keep a constant supply of new individuals'''
    # number of individuals that will die
    sim.stat(pop, popSize=True, subPops=[(0,3)])
    # individuals that will be kept, plus some new guys.
    return pop.popSize() - pop.dvars().popSize + N // 75

def pene(geno, age, ind):
    'Define an age-dependent penetrance function'
    # this disease does not occur in children
    if age < 16:
        return 0
    # if an individual is already affected, keep so
    if ind.affected():
        return 1
    # the probability of getting disease increases with age
    return (0., 0.001*age, 0.001*age)[sum(geno)]

def outputstat(pop):
    'Calculate and output statistics'
    sim.stat(pop, popSize=True, numOfAffected=True,
        subPops=[(0, sim.ALL_AVAIL)],
        vars=['popSize_sp', 'propOfAffected_sp'])
    for sp in range(3):
        print('%s: %.3f%% (size %d)' % (pop.subPopName((0,sp)),
            pop.dvars((0,sp)).propOfAffected * 100.,
            pop.dvars((0,sp)).popSize))
    #
    return True


pop.evolve(
    initOps=[
        sim.InitSex(),
        # random assign age
        sim.InitInfo(lambda: random.randint(0, 75), infoFields='age'),
        # random genotype
        sim.InitGenotype(freq=[0.5, 0.5]),
        # assign an unique ID to everyone.
        sim.IdTagger(),
        sim.PyOutput('Prevalence of disease in each age group:\n'),
    ],
    # increase the age of everyone by 1 before mating.
    preOps=sim.InfoExec('age += 1'),
    matingScheme=sim.HeteroMating([
        # all individuals with age < 75 will be kept. Note that
        # CloneMating will keep individual sex, affection status and all
        # information fields (by default).
        sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
        # only individuals with age between 20 and 50 will mate and produce
        # offspring. The age of offspring will be zero.
        sim.RandomMating(ops=[
            sim.IdTagger(),                   # give new born an ID
            sim.PedigreeTagger(),             # track parents of each individual
            sim.MendelianGenoTransmitter(),   # transmit genotype
        ],
        numOffspring=(sim.UNIFORM_DISTRIBUTION, 1, 3),
        subPops=[(0,1)]),],
        subPopSize=demoModel),
    # number of individuals?
    postOps=[
        sim.PyPenetrance(func=pene, loci=0),
        sim.PyOperator(func=outputstat, step=20)
    ],
    gen = 200
)

# draw two Pedigrees from the last age-structured population
from simuPOP import sampling
sample = sampling.drawNuclearFamilySample(pop, families=2, numOffspring=(2,3),
    affectedParents=(1,2), affectedOffspring=(1,3))
sim.dump(sample)


