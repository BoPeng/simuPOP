#!/usr/bin/env python

#
# $File: reichDemo.py $
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
import math
def demo_model(model, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    model: linear or exponential
    N0:   Initial sim.population size.
    N1:   Ending sim.population size.
    G0:   Length of burn-in stage.
    G1:   Length of sim.population expansion stage.
    '''
    def ins_expansion(gen):
        if gen < G0:
            return N0
        else:
            return N1
    rate = (math.log(N1) - math.log(N0))/G1
    def exp_expansion(gen):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    if model == 'instant':
        return ins_expansion
    elif model == 'exponential':
        return exp_expansion

# when needed, create a demographic function as follows
demo_func = demo_model('exponential', 1000, 100000, 500, 500)
# sim.population size at generation 700
print(demo_func(700))

