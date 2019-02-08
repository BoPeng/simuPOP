#!/usr/bin/env python

#
# $File: advancedVSP.py $
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
pop = sim.Population(size=[2000, 4000], loci=[30], infoFields='x')
# assign random information fields
sim.initSex(pop)
sim.initInfo(pop, lambda: random.randint(0, 3), infoFields='x')
#
# 1, use a combined splitter
pop.setVirtualSplitter(sim.CombinedSplitter(splitters = [
    sim.SexSplitter(),
    sim.InfoSplitter(field='x', values=[0, 1, 2, 3])
]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.MALE
pop.subPopSize([1, 4])    # individuals in sp 1 with value 2 at field x
#
# use a product splitter that defines additional VSPs by sex and info
pop.setVirtualSplitter(sim.ProductSplitter(splitters = [
    sim.SexSplitter(names=['M', 'F']),  # give a new set of names
    sim.InfoSplitter(field='x', values=[0, 1, 2, 3])
]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.MALE with value 1 in sp 0
pop.subPopSize([1, 5])    # sim.FEMALE with value 1 in sp 1
#
# use a combined splitter to join VSPs defined by a
# product splitter
pop.setVirtualSplitter(sim.CombinedSplitter([
    sim.ProductSplitter([
        sim.SexSplitter(),
        sim.InfoSplitter(field='x', values=[0, 1, 2, 3])])],
    vspMap = [[0,1,2], [4,5,6], [7]],
    names = ['Male x<=3', 'Female x<=3', 'Female x=4']))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.MALE with value 0, 1, 2 at field x
pop.subPopSize([1, 1])    # sim.FEMALE with value 0, 1 or 2 at field x

