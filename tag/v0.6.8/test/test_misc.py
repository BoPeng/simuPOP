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

################## PART 1 ###############################
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
#dumper().apply(pop)

# redirected to a file. 
#setLogOutput("session.log")
#dumper().apply(pop)

#setLogOutput()
#dumper().apply(pop)

#print open("session.log").read()
 
#import os
#os.remove("session.log")
