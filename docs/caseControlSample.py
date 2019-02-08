#!/usr/bin/env python

#
# $File: caseControlSample.py $
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
from simuPOP.sampling import drawCaseControlSamples
pop = sim.Population([10000], loci=5)
sim.initGenotype(pop, freq=[0.2, 0.8])
sim.maPenetrance(pop, loci=2, penetrance=[0.11, 0.15, 0.20])
# draw multiple case control sample
samples = drawCaseControlSamples(pop, cases=500, controls=500, numOfSamples=5)
for sample in samples:
    sim.stat(sample, association=range(5))
    print(', '.join(['%.6f' % sample.dvars().Allele_ChiSq_p[x] for x in range(5)]))


