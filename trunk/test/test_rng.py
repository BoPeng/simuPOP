#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of random number selector
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# June, 2004.
#
# Important:  set path to simuPOP.py
#
#  if you do not set PYTHONPATH, you can uncomment
#  the following two lines and specify path to
#  simuPOP.py here.
#

from simuPOP import *
print listAllRNG()
setRNG()
setRNG("ranlux389")
setRNG("random64-libc5")

# help(RNG)

# get the random number generator.
# currently, there is one global RNG used in the system
rng = rng()
print rng.name()
# return 0,...,10
for n in range(1,10):
  print rng.randBinomial(10, .7)
# return 0,1,2,3,4
  
for n in range(1,10):
  print rng.randUniform01()
  
