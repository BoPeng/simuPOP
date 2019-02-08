#!/usr/bin/env python

#
# $File: importExport.py $
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
from simuPOP.utils import importPopulation, export
pop = sim.Population([2,4], loci=5, lociNames=['a1', 'a2', 'a3', 'a4', 'a5'],
    infoFields='BMI')
sim.initGenotype(pop, freq=[0.3, 0.5, 0.2])
sim.initSex(pop)
sim.initInfo(pop, [20, 30, 40, 50, 30, 25], infoFields='BMI')
export(pop, format='fstat', output='fstat.txt')
print(open('fstat.txt').read())
export(pop, format='structure', phenotype='BMI', output='stru.txt')
print(open('stru.txt').read())
pop1 = importPopulation(format='fstat', filename='fstat.txt')
sim.dump(pop1)

