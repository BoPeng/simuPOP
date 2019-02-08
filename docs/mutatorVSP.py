#!/usr/bin/env python

#
# $File: mutatorVSP.py $
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
def fragileX(geno):
    '''A disease model where an individual has increased risk of 
    affected if the number of tandem repeats exceed 75.
    '''
    # Alleles A1, A2.
    maxRep = max(geno)
    if maxRep < 50:
        return 0
    else:
        # individuals with allele >= 70 will surely be affected
        return min(1, (maxRep - 50)*0.05)

def avgAllele(pop):
    'Get average allele by affection sim.status.'
    sim.stat(pop, alleleFreq=(0,1), subPops=[(0,0), (0,1)],
        numOfAffected=True, vars=['alleleNum', 'alleleNum_sp'])
    avg = []
    for alleleNum in [\
            pop.dvars((0,0)).alleleNum[0],  # first locus, unaffected
            pop.dvars((0,1)).alleleNum[0],  # first locus, affected
            pop.dvars().alleleNum[1],       # second locus, overall
        ]:
        alleleSum = numAllele = 0
        for idx,cnt in enumerate(alleleNum):
            alleleSum += idx * cnt
            numAllele += cnt
        if numAllele == 0:
            avg.append(0)
        else:
            avg.append(alleleSum * 1.0 /numAllele)
    # unaffected, affected, loc2
    pop.dvars().avgAllele = avg
    return True

pop = sim.Population(10000, loci=[1, 1])
pop.setVirtualSplitter(sim.AffectionSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(genotype=[50, 50])
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        # determine affection sim.status for each offspring (duringMating)
        sim.PyPenetrance(func=fragileX, loci=0),
        # unaffected offspring, mutation rate is high to save some time
        sim.StepwiseMutator(rates=1e-3, loci=1),
        # unaffected offspring, mutation rate is high to save some time
        sim.StepwiseMutator(rates=1e-3, loci=0, subPops=[(0, 0)]),
        # affected offspring have high probability of mutating upward
        sim.StepwiseMutator(rates=1e-2, loci=0, subPops=[(0, 1)],
           incProb=0.7, mutStep=3),
        # number of affected
        sim.PyOperator(func=avgAllele, step=20),
        sim.PyEval(r"'Gen: %3d #Aff: %d AvgRepeat: %.2f (unaff), %.2f (aff), %.2f (unrelated)\n'"
            + " % (gen, numOfAffected, avgAllele[0], avgAllele[1], avgAllele[2])",
            step=20),
    ],
    gen = 101
)

