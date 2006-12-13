#!/usr/bin/env python
#
# Purpose:
#    testing the mpi version of simuPOP
#
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision: 475 $
# $LastChangedDate: 2006-10-04 01:00:38 -0500 (Wed, 04 Oct 2006) $
# 


import simuOpt
simuOpt.setOptions(quiet=True, optimized=False, mpi=True)

from simuPOP import *
import unittest, os, sys, exceptions, time

class TestMPI(unittest.TestCase):

    def TestMPIStart(self):
        '''Start mpi'''
        testMPI()


if __name__ == '__main__':
    unittest.main()
