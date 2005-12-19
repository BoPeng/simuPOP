#!/usr/bin/env python
#             Perforce Defect Tracking Integration Project
#              <http://www.ravenbrook.com/project/p4dti/>
#
#                    RUN_TESTS.PY -- RUN ALL TESTS
#
#             Gareth Rees, Ravenbrook Limited, 2001-03-14
#
#
# 1. INTRODUCTION
#
# This script finds all test scripts in the current directory and runs
# them.
#
# Test scripts should be called test_MODULE.py (for unit tests, MODULE
# should be the module being tested.  Each test script should define a
# function called "tests", which returns an object belonging to the
# unittest.TestSuite class.  See [PyUnit] for details of the Python Unit
# Test Framework.

import dircache
import re
import unittest


# 2. FIND TESTS
#
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


# A. REFERENCES
#
# [PyUnit] "PyUnit - a unit testing framework for Python"; Steve
# Purcell; <http://pyunit.sourceforge.net/>.
#
#
# B. DOCUMENT HISTORY
#
# 2001-03-14 GDR Created.
#
# 2002-10-25 NB Add verbosity control to get more output from the unit
# test framework.
#
#
# C. COPYRIGHT AND LICENSE
#
# This file is copyright (c) 2001 Perforce Software, Inc.  All rights
# reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1.  Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
# 2.  Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDERS AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
#
