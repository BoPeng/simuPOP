#!/usr/bin/env python
#
# Purpose:
#  testing other features of simuPOP
#  Some of them are advanced.
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
#
#import sys
#sys.path.append('../lib/debug')
# 
#

from simuPOP import *
pop=population(size=10, loci=[3,5])


################## PART 1 ###############################
#
# details about shared variables.
#
# Shared variables
# - have type string, double, array (of double) or dictionary
# - can be set using setVar function
# - can be retrieved using getVar function
# - usually set by statistics calculators. They
#   usually do not generate any output but set certain
#   shared variables.
# - all file names can have shared variables in the form
#   of %name $(name), index take form of [] or array and
#   {} for dictionary
# - output can display shared variables, including
#   shared array. (separator from output will be used
#   to separate array elements.
# - DataSource objects can collect shared variables and
#   provide data to whoever needs them.
#
#########################################################

#
# area=0 is the public, default area.
#

setStringVar("str1","name",area=0)
listVars(0)
setDoubleVar("double1",111)
setDoubleVar("double1_1",12, area=1)
setArrayVar("array1",[1,2,3,4])
setStringVar("str2","a string")
# a dictionary is a list of paired values (int:double)
setDictionaryVar("dict1",{1:1,2:3,4:5.6})
listVars(1)


# modify the value
setStringVar("str1","modified")
listVars(0)

# can set to area 1 (for replicate 1)
setArrayVar("array1_1", [1,2,3],1)

# 
print getStringVar("str1")
# you can not get a double var using getStringVar
#print getStringVar("double1")
print getDoubleVar("double1")
print getArrayVar("array1")
# double1_1 does not exist in area 0. 
#print getDoubleVar("double1_1")
listVars(0)
# since v4 is an array, we need to specify index
print getDoubleVar("double1_1",1)

# we can get the values of these variables by
print getStringVar("str1")
print getDoubleVar("double1")
print getArrayVar("array1")
print getDictionaryVar("dict1")

#
# now, in operator, how can we use these variables?
#
# Rules:
# %name or %(name)
#
# %grp, %rep, %gen is explained as group, replicate, generation
# others are shared variables.
# %name  string or double, or whole array/dictionary
# %name[index] array
# %name{index} dictionary
#
output("%gen%(str1)a%%\n").apply(pop)
output("%(gen)a%(str1)a%%\n").apply(pop)
# v4 is not in area 0
#output("%(double1_1)\n").apply(pop)
#
output("%double1_1\n").apply(pop,rep=1)
#
output("%(array1[0])\n").apply(pop,rep=1)
# can not specify index for scalar variable
#output("%(array1{2})a\n").apply(pop)
# can not use {} for array
#output("%(array1[2])a\n").apply(pop)
# can not use [] for dictionary
#output("%dict1[3]\n").apply(pop)
output("%(dict1{2})a\n").apply(pop)
#
# or all the values
output("%(dict1)a\n").apply(pop)


################## PART 3 ###############################
#
# save dedug output to another file.
#
# There are several outputs:
#  - python output --- can not be redirected
#  - simuPOP output that goes to standard output
#    including debug informations and output from dumper()
#    --- can be directed
#
##########################################################

# outputset simuPOP log info (not python sys.stdout) to a file
# the following goes to default output
dumper().apply(pop)

# redirected to a file. 
setLogOutput("session.log")
dumper().apply(pop)

setLogOutput()
dumper().apply(pop)

print open("session.log").read()
 
import os
os.remove("session.log")
