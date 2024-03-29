#!/usr/bin/env python

#
# $File: chromType.py $
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
pop = sim.Population(size=6, ploidy=2, loci=[3, 3, 3, 2, 2, 4, 4],
    chromTypes=[sim.AUTOSOME]*2 + [sim.CHROMOSOME_X, sim.CHROMOSOME_Y, sim.MITOCHONDRIAL]
        + [sim.CUSTOMIZED]*2)
sim.initGenotype(pop, freq=[0.3, 0.7])
sim.dump(pop, structure=False) # does not display genotypic structure information
