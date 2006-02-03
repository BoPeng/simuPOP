#!/usr/bin/env python

############################################################################
#    Copyright (C) 2004 by Bo Peng                                         
#    bpeng@rice.edu                                                        
#                                                                          
#    $LastChangedDate$          
#    $Rev$                       
#                                                                          
#    This program is free software; you can redistribute it and/or modify  
#    it under the terms of the GNU General Public License as published by  
#    the Free Software Foundation; either version 2 of the License, or     
#    (at your option) any later version.                                   
#                                                                          
#    This program is distributed in the hope that it will be useful,       
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
#    GNU General Public License for more details.                          
#                                                                          
#    You should have received a copy of the GNU General Public License     
#    along with this program; if not, write to the                         
#    Free Software Foundation, Inc.,                                       
#    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
############################################################################

#
# Purpose:
#  testing pause operator
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# Dec, 2005
#

from simuPOP import *
from simuUtil import *

import unittest

class TestPause(unittest.TestCase):
  def testPauseAtGen(self):
    simu = simulator( population(size=10, ploidy=2, loci=[2, 3]),
      randomMating(), rep=5)
    print "\n\nUSER INTERACTION: Please press q\n\n"
    self.assertRaises( exceptions.SystemError, simu.evolve, 
      ops=[ pause(at=[10]), 
            # should quite, can not reach generation 12
            terminateIf("True", at=[12] ) ] )
      
  def testExitToShell(self):
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
