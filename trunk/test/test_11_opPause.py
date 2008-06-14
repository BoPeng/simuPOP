#!/usr/bin/env python
#
# Purpose:
#  testing pause operator
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys, exceptions

#NOTE:
#
# These tests will not be automatically run
# since they need user-interaction.
#
# To run the test:
# test_11_opPause.py TestPause.interactiveTestPauseAtGen
# test_11_opPause.py TestPause.interactiveTestExitToShell
#

class TestPause(unittest.TestCase):

  def interactiveTestPauseAtGen(self):
    'Testing resume to simulation'
    simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
      randomMating(), rep=5)
    print "\n\nUSER INTERACTION: Please press q\n\n"
    self.assertRaises( exceptions.SystemError, simu.evolve,
      ops=[ pause(at=[10]),
            # should quite, can not reach generation 12
            terminateIf("True", at=[12] ) ] )

  def interactiveTestExitToShell(self):
    'Testing exit to a shell'
    simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
      randomMating(), rep=5)
    print "\n\nUSER INTERACTION: Please press s and then Ctrl-D"
    print "Please check the existence of variable pop\n\n"
    simu.evolve(
      ops=[ pause(at=[10]) ], end=12)
    print "\n\nUSER INTERACTION: Please press s and then Ctrl-D"
    print "Please check the existence of variable tmpPop\n\n"
    simu.evolve(
      ops=[ pause(at=[20], popName='tmpPop') ], end=25)

if __name__ == '__main__':
  unittest.main()
