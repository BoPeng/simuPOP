#!/usr/bin/env python
#
# This script finds all test scripts in the current directory and runs
# them.
#
# NOTE: to run test on different versions of libraries,
#   setenv SIMUALLELETYPE short
#                         long
#                         binary
#
# NOTE: you can stop testing test_20_rpy with an arbitrary commandline
#   argument.
#
import dircache, re, unittest, sys, os

# Find all files with names like test_MODULE.py.  For each such file,
# import it as a module.  If it defines a function called "tests", call
# that function and add its result to the test suite.

tests = unittest.TestSuite()

for file in dircache.listdir('.'):
  match = re.match("^(test_(.*))\\.py$", file)
  if match:
    m = match.group(1)
    if len(sys.argv) > 1 and m == 'test_20_rpy':
      continue
    print "Adding test cases in ", m
    module = __import__(m)
    tests.addTest( unittest.defaultTestLoader.loadTestsFromModule( module ) )

# 3. RUN TESTS
#
# Set verbosity=2 to get more output.  Sadly it's not possible to get
# the output of failed tests until all the tests are completed.

test_runner = unittest.TextTestRunner(verbosity=2)
test_runner.run(tests)
