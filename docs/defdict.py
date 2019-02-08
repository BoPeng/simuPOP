#!/usr/bin/env python

#
# $File: defdict.py $
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
pop = sim.Population([100]*2, loci=1)
sim.initGenotype(pop, freq=[0, 0.2, 0.8], subPops=0)
sim.initGenotype(pop, freq=[0.2, 0.8], subPops=1)
sim.stat(pop, alleleFreq=0, vars=['alleleFreq_sp'])
for sp in range(2):
    print('Subpop %d (with %d alleles): ' % (sp, len(pop.dvars(sp).alleleFreq[0])))
    for a in range(3):
        print('%.2f ' % pop.dvars(sp).alleleFreq[0][a])


