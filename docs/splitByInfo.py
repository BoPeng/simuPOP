#!/usr/bin/env python

#
# $File: splitByInfo.py $
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit https://github.com/BoPeng/simuPOP
# for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (Bo.Peng@bcm.edu)
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
# the user's guide (https://github.com/BoPeng/simuPOP/manual) for a detailed
# description of this example.
#

import simuPOP as sim
import random
pop = sim.Population([1000]*3, subPopNames=['a', 'b', 'c'], infoFields='x')
pop.setIndInfo([random.randint(0, 3) for x in range(1000)], 'x')
print(pop.subPopSizes())
print(pop.subPopNames())
sim.splitSubPops(pop, subPops=[0, 2], infoFields=['x'])
print(pop.subPopSizes())
print(pop.subPopNames())
