#!/usr/bin/env python

#
# $File: run_tests.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit http://simupop.sourceforge.net
# for details.
#
# Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
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


# This script finds all test scripts in the current directory and runs them.
# 
# Usage:
#   run_tests.py
#       Run tests for 'short', 'long', and 'binary' modules and summarize
#       outputs
#
#   run_tests.py [ short | long | binary ]
#       Run the tests for specified modules.
#
import dircache, re, unittest, sys, os

from simuOpt import setOptions
if len(sys.argv) > 1:
    setOptions(alleleType = sys.argv[1])

def importTests():
    tests = unittest.TestSuite()

    # Find all files with names like test_MODULE.py.  For each such file,
    # import it as a module.  If it defines a function called "tests", call
    # that function and add its result to the test suite.
    for file in dircache.listdir('.'):
        match = re.match("^(test_(.*))\\.py$", file)
        if match:
            m = match.group(1)
            print "Adding test cases in ", m
            module = __import__(m)
            tests.addTest( unittest.defaultTestLoader.loadTestsFromModule( module ) )
    return tests

if __name__ == '__main__':
    if len(sys.argv) == 1:
        from subprocess import Popen
        for allele in ['binary', 'short', 'long']:
            print '\n\n===== Testing %s module =====\n\n' % allele
            p = Popen('python %s %s' % (sys.argv[0], allele), bufsize=0, shell=True)
            os.waitpid(p.pid, 0)
    else:
        test_runner = unittest.TextTestRunner(verbosity=2)
        test_runner.run(importTests())
