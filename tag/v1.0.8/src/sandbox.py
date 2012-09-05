#!/usr/bin/env python

#
# $File: simuPOP.py $
# $LastChangedDate: 2010-05-26 11:09:40 -0700 (Wed, 26 May 2010) $
# $Rev: 3569 $
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

"""
simuPOP sandbox. This module contains classes (operators and mating schemes)
and functions that are either experimental (and have a chance to be formally
added to simuPOP), or designed for special purposes. Because of the temporary
nature of these classes, it is recommended that modules that use these classes
include Python version of them so that they can be used when the sandbox
version is no longer available.
"""

__all__ = []

__long__all__ = [
    'revertFixedSites',
    'RevertFixedSites',
    'MutSpaceMutator',
    'MutSpaceSelector',
    'MutSpaceRecombinator',
]

# get options
from simuOpt import simuOptions
import os, sys

#
# These classes are available in binary modules but are not exposed
# in simuPOP. Existing functions are only available in long module.
#
if simuOptions['Optimized']:
    if simuOptions['AlleleType'] == 'short':
        pass
    elif simuOptions['AlleleType'] == 'long':
        from simuPOP_laop import sb_RevertFixedSites as RevertFixedSites
        from simuPOP_laop import sb_MutSpaceMutator as MutSpaceMutator
        from simuPOP_laop import sb_MutSpaceSelector as MutSpaceSelector
        from simuPOP_laop import sb_MutSpaceRecombinator as MutSpaceRecombinator
        __all__.extend(__long__all__)
    elif simuOptions['AlleleType'] == 'binary':
        pass
    elif simuOptions['AlleleType'] == 'mutant':
        pass
    else:
        pass
else:
    if simuOptions['AlleleType'] == 'short':
        pass
    elif simuOptions['AlleleType'] == 'long':
        from simuPOP_la import sb_RevertFixedSites as RevertFixedSites
        from simuPOP_la import sb_MutSpaceMutator as MutSpaceMutator
        from simuPOP_la import sb_MutSpaceSelector as MutSpaceSelector
        from simuPOP_la import sb_MutSpaceRecombinator as MutSpaceRecombinator
        __all__.extend(__long__all__)
    elif simuOptions['AlleleType'] == 'binary':
        pass
    elif simuOptions['AlleleType'] == 'mutant':
        pass
    else:
        pass

def revertFixedSites(pop):
    '''Apply operator ``RevertFixedSites`` to ``pop``'''
    RevertFixedSites().apply(pop)
