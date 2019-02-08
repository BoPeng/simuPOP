#!/usr/bin/env python

#
# $File: pedigreeMating.py $
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
# create a population without any genotype
from simuPOP.utils import migrSteppingStoneRates
ped = sim.Population(size=[1000]*5, ancGen=-1, 
    infoFields=['ind_id', 'father_id', 'mother_id', 'migrate_to'])
ped.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
    ],
    preOps=sim.Migrator(rate=migrSteppingStoneRates(0.1, 5)),
    matingScheme=sim.RandomMating(
        numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
        ops=[
            # we do not even need a genotype transmitter...
            sim.IdTagger(),
            sim.PedigreeTagger(),
        ]),
    gen=100
)
# convert itself to a pedigree object
ped.asPedigree()
# we should have 100 ancestral generations
N = ped.ancestralGens()
# We should have 101 * 1000 * 5 individuals, but how many actually
# contribute genotype to the last generation?
anc = ped.identifyAncestors()
len(anc)
# remove individuals who do not contribute genotype to the last generation
allIDs = [x.ind_id for x in ped.allIndividuals()]
removedIDs = list(set(allIDs) - set(anc))
ped.removeIndividuals(IDs=removedIDs)
# now create a top most population, but we do not need all of them
# so we record only used individuals
IDs = [x.ind_id for x in ped.allIndividuals(ancGens=N)]
sex = [x.sex() for x in ped.allIndividuals(ancGens=N)]
# create a population, this time with genotype. Note that we do not need
# populaton structure because PedigreeMating disregard population structure.
pop = sim.Population(size=len(IDs), loci=1000, infoFields='ind_id')
# manually initialize ID and sex
sim.initInfo(pop, IDs, infoFields='ind_id')
sim.initSex(pop, sex=sex)
pop.evolve(
    initOps=sim.InitGenotype(freq=[0.4, 0.6]),
    # we do not need migration, or set number of offspring,
    # or demographic model, but we do need a genotype transmitter
    matingScheme=sim.PedigreeMating(ped, 
        ops=sim.MendelianGenoTransmitter()),
    gen=100
)
# let us compare the pedigree and the population object
print(ped.indInfo('ind_id')[:5])
print(pop.indInfo('ind_id')[:5])
print([ped.individual(x).sex() for x in range(5)])
print([pop.individual(x).sex() for x in range(5)])
print(ped.subPopSizes())
print(pop.subPopSizes())

