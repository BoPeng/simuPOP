#!/usr/bin/env python

#
# $File: sexSpecificRec.py $
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

from simuPOP import *
class sexSpecificRecombinator(PyOperator):
    def __init__(self, intensity=0, rates=0, loci=[], convMode=NO_CONVERSION,
            maleIntensity=0, maleRates=0, maleLoci=[], maleConvMode=NO_CONVERSION,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.Recombinator = Recombinator(rates, intensity, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = Recombinator(maleRates, maleIntensity,
            maleLoci, maleConvMode)
        #
        PyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)
    #
    def transmitGenotype(self, pop, off, dad, mom):
        # Form the first homologous copy of offspring.
        self.Recombinator.transmitGenotype(mom, off, 0)
        # Form the second homologous copy of offspring.
        self.maleRecombinator.transmitGenotype(dad, off, 1)
        return True

pop = Population(10, loci=[15]*2, infoFields=['father_idx', 'mother_idx'])
pop.evolve(
    initOps=[
        InitSex(),
        InitGenotype(freq=[0.4] + [0.2]*3)
    ],
    matingScheme=RandomMating(ops=[
        sexSpecificRecombinator(rates=0.1, maleRates=0),
        ParentsTagger()
    ]),
    postOps=Dumper(structure=False),
    gen = 2
)

