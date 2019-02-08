#!/usr/bin/env python

#
# $File: mtDNA_evolve.py $
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

def alleleCount(pop):
    summary = [0]* 6
    for ind in pop.individuals():
        geno = ind.genotype(ploidy=0)
        summary[geno[0] + geno[2] + geno[4] + geno[6] + geno[8]] += 1
    print('%d %s' % (pop.dvars().gen, summary))
    return True

pop = sim.Population(1000, loci=[2]*5, chromTypes=[sim.CUSTOMIZED]*5)
pop.evolve(
    # every one has miDNAs 10, 00, 00, 00, 00
    initOps=[
        sim.InitGenotype(haplotypes=[[1]+[0]*9]),
    ],
    # random select cells for cytoplasmic segregation
    matingScheme=sim.RandomSelection(ops= [
        sim.MitochondrialGenoTransmitter(),
    ]),
    postOps=sim.PyOperator(func=alleleCount, step=10),
    gen = 51
)

