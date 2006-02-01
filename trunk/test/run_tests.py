#!/usr/bin/env python
#
# This script finds all test scripts in the current directory and runs
# them.
#
# Test scripts should be called test_MODULE.py (for unit tests, MODULE
# should be the module being tested.  Each test script should define a
# function called "tests", which returns an object belonging to the
# unittest.TestSuite class.  See [PyUnit] for details of the Python Unit
# Test Framework.

import dircache, re, unittest, sys
import simuOpt

# if no option, test standard libraries
if len(sys.argv) > 1:
  for arg in sys.argv[1:]:
    if arg in ['standard', 'short', 'long', 'binary']:
      simuOpt.setOptions(AlleleType = sys.argv[1])
    elif srg == 'optimized':
      simuOpt.setOptions(optimized=True)

# Find all files with names like test_MODULE.py.  For each such file,
# import it as a module.  If it defines a function called "tests", call
# that function and add its result to the test suite.

tests = unittest.TestSuite()

for file in dircache.listdir('.'):
    match = re.match("^(test_(.*))\\.py$", file)
    if match:
        module = __import__(match.group(1))
        if hasattr(module, 'tests'):
            tests.addTest(getattr(module, 'tests')())


# 3. RUN TESTS
#
# Set verbosity=2 to get more output.  Sadly it's not possible to get
# the output of failed tests until all the tests are completed.

test_runner = unittest.TextTestRunner(verbosity=2)

test_runner.run(tests)


