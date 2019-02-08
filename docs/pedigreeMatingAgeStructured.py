#!/usr/bin/env python

#
# $File: pedigreeMatingAgeStructured.py $
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
pop = sim.Population(N, infoFields=['age', 'ind_id', 'father_id', 'mother_id'])
# we simulate age 0, 1, 2, 3 
pop.setVirtualSplitter(sim.InfoSplitter(field='age', values=[0, 1, 2, 3]))
pop.evolve(
    initOps=[
        sim.InitSex(),
        # random assign age
        sim.InitInfo(lambda: random.randint(0, 3), infoFields='age'),
        # random genotype
        sim.InitGenotype(freq=[0.5, 0.5]),
        # assign an unique ID to everyone.
        sim.IdTagger(),
    ],
    # increase the age of everyone by 1 before mating.
    preOps=sim.InfoExec('age += 1'),
    matingScheme=sim.HeteroMating([
        # age 1, 2 will be copied
        sim.CloneMating(
            ops=[
                # This will set offspring ID
                sim.CloneGenoTransmitter(),
                # new ID for offspring in order to track pedigree
                sim.IdTagger(),
                # both offspring and parental IDs will be the same
                sim.PedigreeTagger(output='>>structured.ped'),
            ],
            subPops=[(0,1), (0,2)],
            weight=-1
        ),
        # age 2 produce offspring
        sim.RandomMating(
            ops=[
                # new ID for offspring
                sim.IdTagger(),
                # record complete pedigree
                sim.PedigreeTagger(output='>>structured.ped'),
                sim.MendelianGenoTransmitter(),   # transmit genotype
            ],
            subPops=[(0,2)]
        )]
    ),
    gen=20
)

# use a pedigree object recovered from a file saved by operator PedigreeTagger
ped = sim.loadPedigree('structured.ped')
# create a top most population, but we do not need all of them
# so we record only used individuals
IDs = [x.ind_id for x in ped.allIndividuals(ancGens=ped.ancestralGens())]
sex = [x.sex() for x in ped.allIndividuals(ancGens=ped.ancestralGens())]
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
        ops=sim.IfElse(lambda mom: mom is None,
                sim.CloneGenoTransmitter(),
                sim.MendelianGenoTransmitter())
    ),
    gen=100
)
# 
print(pop.indInfo('ind_id')[:5])
print([pop.individual(x).sex() for x in range(5)])
# The pedigree object does not have population structure
print(pop.subPopSizes())

