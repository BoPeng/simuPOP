#!/usr/bin/env python
# test the getParam function of simuUtil

from simuUtil import *
from types import *

def usage():
  return "usage"

#getParam(arg, longarg,name,prompt,default,allowedTypes)
print getParam(arg='d:',prompt='enter distance: ',
        allowedTypes=[IntType,LongType,FloatType])

