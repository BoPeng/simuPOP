#!/usr/bin/env python
#
# Purpose:
#   testing inline expression for simupoop.
#
# Features:
#   data type:
#     int        34
#     double     3.4e5, .05
#     string     "string", 'atring'
#     array      [3,4,5], [4,5,6.4]  
#     dictionary {2:4.5, 10:8, 1020:9.45}
#
#   ops:
#     +, - , *, /, ==, !=, >, <, >=, <=, ^, &&, ||
#     [], {}   vector/dictionary index.
#     prefix --/++ 
#
#     All in the same meaning as those in C/C++
#     (please use () to group expressions since
#      precedence is not implemented strictly)
#
#     operations on different types are stricted
#     in the case of no ambiguity.
#
#   functions: (use listFunctions() function to list )
#     sin, cos, len, any, all, ... 
#
#   if/loop
#     if/else if/else
#     for(;;){}
#     foreach( %i in vector){ }
#     while(){}
#     return
#     break
#     continue
#
#   variable:
#     %a, %a[], %a{}
#     no type definition
#     all variables are glocal. Can be displayed by listVars(rep)
#
#   function call
#     name(arg1, arg2...)
#
#   function declaration
#     function name(arg1,arg2,...){ ...  }
#
#     function definition is permanent.
#
#   
#
# Major differences with C syntax
#   1. everything is vector or dictionary
#   2. all variables are glocal. No concept of scope.
#   3. no type (automatic type conversion)
#   4. no suffix ++/--, bit ops ( |, &, >> etc )
#   5. function can only pass values.
#
#
#
# Author:
#   Bo Peng (bpeng@rice.edu)
#
# Important:  set path to simuPOP.py
#
#  if you do not set PYTHONPATH, you can uncomment
#  the following two lines and specify path to
#  simuPOP.py here.
#
#import sys
#sys.path.append('../lib/debug')
# 
# load module.
from simuPOP import *

a="""
b
c
d
d
"""

#
simu = simulator(
    population(size=10, ploidy=2, loci=[2, 3 ]), 
    randomMating(), rep=5)
# this will not change generation number
simu.apply([ initByFreq([.2, .3, .5]) ])

#
# let us first set up some variables, but no output
#
simu.setGen(0)
simu.evolve([ alleleCounter(alleles=[[0]]) ], end=10)
#
# let us do some calculation
# anyway, you should be using the following
# within an operator, but for now, let us see
# what expression can do
#
# basic calculator like expressions but
# calculator and dictionary can be used.
#
print Expression("5+6").valueAsString()
#
# this is the full form, we can use a shortcut function
#
print evaluate("5+(6*6+2.3)")
#
# involve a vector?
# vector length should be the same
#print evaluate("[1,1]+2.*(9+[3,4,5])")
print evaluate("[[1,1],2]+2.*(9+[3,4*4,5])")
# we can use shared variable
print evaluate("%alleleNum0+1")
# string can be used.
print evaluate("'string'")
# no such variable, will create one.
print evaluate("%num")
print evaluate("[1, [1,2,3],1]")
#
# index, subset of a dictionary is a dictionary
print evaluate("%alleleNum0[1]+1")
# negative index, index will be the key if dictionary  is involved.
print evaluate("%alleleNum0[-3]")
#
# for vector, things are a bit awkard. since there is no -0
print evaluate("%a=[1,2,3];%a[-1]")
print evaluate("%a=[1,2,3]; %a[0]")
# multiple index
print evaluate("%a=[1,2,3];%b=[1,2]; %a[%b]")
# note that the index should be a vector...
print evaluate("%a=[1,2,3];%a[[1,2]]")
print evaluate("%a=[1,2,3];%a[[-1,-2]]")
# can not mix pos/neg indices
#print evaluate("%a=[1,2,3];%a[[1,-1]]")
# index can be applied to any expression
print evaluate("%a=[1,2,3];(%a+1)[[-1,-2]]")
print evaluate("%a=[1,2,3];sin(%a)[1]")
#
# this one does look weird.
print evaluate("%a=[1,2,3][1]")
#
# there are other kind of assignments as well
#
print  evaluate("%b+=1")
print  evaluate("%b-=1")
print  evaluate("%b*=[2,4]")
print  evaluate("%b/=[2.5,4]")
print  evaluate("%b[1]+=1")
# vector index version is not supported as yet.
#print  evaluate("%b[[1,2]]+=1")
# function
print evaluate("sum(%alleleNum0)+1")
print evaluate("len(%alleleNum0)")
print evaluate("mean(%alleleNum0)")
print evaluate("sin(%alleleNum0)")

print evaluate("sin(10)")
print evaluate("sin(10.0)")
print evaluate("min([-10.0, %alleleNum0])")
print evaluate("max([-10.0, %alleleNum0])")
#
print evaluate("any([1,0,0])")
print evaluate("all([1,0,0])")
print evaluate("asInt([5.6, sin(6), 7.5][2])")
print evaluate("asArray(5)")
print evaluate("asArray({1:5,2:3})")
print evaluate("asDictionary([2,3])")
print evaluate("asString({1:5,2:4})")
# all cross products.
print evaluate("prod([1,2,3],[4,5])")
# exclude homozygous items (key1==key2)
print evaluate("prod({1:2,2:3,3:2},{1:2,2:4},1)")
#
#
# default to all
print evaluate("if([1,0,0]){ [4,5,6]; }")
print evaluate("if([1,1,1]){ [4,5,6]; }")
print evaluate("if(!all([1,0.1,0])){ [4,5,6]; }")

print evaluate("![1,0,0]")
print evaluate("-[1,0.5,0]")

# you can use prefix ++/-- to a variable.
#
# NOTE:
#  1. prefix ++/-- can ONLY work with variables as a whole
#       and inc/dec them
#
#  2. there is NO suffix ++/--
#     This is because the behavior of suffix ++/-- is
#     difficult to implement as that of C/C++. And
#     different behavior will need to confusion.
#
print evaluate("++%alleleNum0")
listVars(1)
print evaluate("--%alleleNum0")
listVars(1)
#print evaluate("%alleleNum0++")
#
# element accss/change
print evaluate("%alleleNum0[1]")
print evaluate("++%alleleNum0[1]")
#

# it is recommended that you access dictinary by {}
print evaluate("++%alleleNum0{2}")
listVars(1)
#
# you can not do these to a value.
#print evaluate("++[1,0,0]")
#print evaluate("--[1,0,0]")
print evaluate("[1,0,0] || [0,1,0]")
print evaluate("[1,0,0] && [0,1,0]")
print evaluate("[1,1,0] > [0,1,0]")
print evaluate("[1,1,0] >= [0,1,0]")
print evaluate("[1,1,0] < [0,1,0]")
print evaluate("[1,1,0] <= [0,1,0]")
print evaluate("[1,1,0] != [0,1,0]")
print evaluate("[1,1,0] == [0,1,0]")
# note that vectors only store floating values so
# the floating point comparison is not safe.
# We adopt the comparison system that consider
# a and a+eps the same. (eps=1e-9 in simuPOP.
# therefore
print evaluate("[1,1,0] == [0,1+1e-10,0]")
print evaluate("[1,1,0] == [0,1+1e-8,0]")
# however, this rule does not apply to internal formats
# (non-array types) which differ 1 and 1+e-15 but
# not 1 and 1+e-16
print evaluate("1 == (1+1e-8)")

print evaluate("[5,3,2]^4")
print evaluate("[5,3,2]^.5")
#turnOnDebug(DBG_UTILITY)
#turnOffDebug(DBG_UTILITY)
# assignment
print evaluate("%d=5*sin(8)")
# d is no as array type
#print evaluate("%d[1]=5*sin(8)")
print evaluate("%alleleNum0[1]=5*sin(8)")
#
# if
print evaluate("if(0){%d=5;}else{%d=-5;}")
print evaluate("if(1){%d=5;}else{%d=-5;}")
#
# for
#
# AGAIN: do not use %f++!
#
print evaluate("for(%f=1;%f<5;++%f){ print(%f,'\n'); }")
print evaluate("for(%f=1;%f<5;++%f){ print(%f,'\n'); }")
print evaluate("%f=1;while(%f<5){ ++%f; print(%f,'\n');}")
#
# foreach
print evaluate("%f={1:2,3:5.6}; foreach(%v in %f){ print(%v,' '); }")

#
# what if I have a long program?
#
# Comments:
#
# Both C and C++ style comments
#   /* */ and // are offered
#
# However, python will translate \n before it is
# passed to the parser so the parser will fail
# in the case of
#  // print("\n");
# 
# use // will caution!
#
program = """
/* C-style comment is allowed. */
%a=[1,2,3,4,5];
for(%i=1; %i<len(%a); ++%i)
  %a[%i]+= %a[%i-1];
%a;
"""

print evaluate(program)

# although it is unlikely that the program will be
# long enough to justify a separate file, you can
# have your library of statistics and call them
#
# Note: function() definition has not been implemented
# you have to do everything with global variables.
#
FILE = open("program.txt","w")
FILE.write(program)
FILE.close()

print evaluate(file="program.txt")

import os
os.remove("program.txt")

#
########### some bigger programs ##########
#
# return statement
#
program = """
%a=[];
%j = 0;
for(%i=1; %i<10; ++%i)
{
  %j += %i;
  %a = [ %a, %j];
  print(%i, ": ", %a,"\n");
  if( %j > 10 )
     return %a;
  else
     --%j;
 }
%j;
"""

print evaluate(program)


# if
# else if
# else
#
program = """
%a=[];
%j = 0;
for(%i=1; %i<10; ++%i)
{
  %j += %i;
  %a = [ %a, %j];
  print(%i, ": ", %a,"\n");
  if( %j > 10 )
     return %a;
  else  if( %j == 4 )
    ++%j;   
  else
    %j += 6;
}
%j;
"""

print evaluate(program)

# 
# continue and break
# else
#
program = """
%a=[];
for(%i=1; %i<10; ++%i)
{
  if( len(%a) > 5 )
    break;
  else if( %i == 2) // skip 2
    continue;
  
  %a = [ %a, %i];
  print(%i, ": ", %a, "\n"); 
}
print("Final result: ", %a, "\n");
"""

print evaluate(program)

#
#
# function declaration
#
# Note: function declaration is permanent.
#
# i.e., once a function is defined, it can be used
# (or replced by another function with the same name.)
#
# to see a list of functions, use listFunctions() function
program = """
function test(b)
"this is a test for function call"
{
return b*5;
}
%a=[1,2,3];
test(%a);
"""

print evaluate(program)
listFunctions()

program = """
function test(a)
" re-definition of test function"
{
return a+5;
}
%a=[1,2,3];
test(%a);
"""
print evaluate(program)
listFunctions()

#
# simupop will provide a function library so that
# you can calculate a lot statistics by
# something like
#evaluate(file="stat.def").evaluate()

#
#
# we can also call python from within simuPOP
print pyEvalExpression("1+1")
print pyEvalExpression("%gen+1")

# if we do not want to replace var
print pyEvalExpression("%gen+1", replaceVars=0)

# this is wrong. since b=1+1 is a staement, not an expression
# print pyEvalExpression("b=1+1")

#
# should treat as statement
# note that statement does not have return value
pyEvalStatement("b=1+1")
print b


# What is the point of evaluating python
# code in simuPOP?
#
# use as ops!
pop=population(1)
simu=simulator(pop,randomMating())
simu.evolve( [ pyExpression("%gen")], end=10)

