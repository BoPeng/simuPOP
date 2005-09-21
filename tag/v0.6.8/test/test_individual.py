#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of individual 
#  of simuPOP
#
# Author:
#  Bo Peng (bpeng@rice.edu)
#
# June, 2004.
# 
# load module
from simuPOP import *
# for its testing functions and more
from simuUtil import *
from exceptions import *

# load population for testing.
pop = loadPopulation("pop1.bin", format="bin")
ind = pop.individual(0)

# test table:
#
# expr, statements, expected result, expected exception
tests = [
  ['', 'pop.individual(-1)', '', NotImplementedError],
  ['', 'pop.individual(0)', '', ''],
  ['', 'pop.individual(20, 1)', '', ''],
  ['', 'pop.individual(20, 0)', '', IndexError],  
  ['', 'pop.individual(90, 1)', '', IndexError],
  ['', 'ind=pop.individual(0)', '', ''],
  ['ind.ploidy()', '', 2, ''],
  ['ind.ploidyName()', '', 'diploid', ''],
  ['ind.numChrom()', '', 2, ''],
  ['ind.numLoci(5)', '', '', IndexError],
  ['ind.numLoci(0)', '', 5, ''],
  ['ind.numLoci(1)', '', 7, ''],
  ['ind.numLoci(2)', '', '', IndexError],
  ['ind.locusDist(10)', '', 6, ''],
  ['ind.absLocusIndex(1,5)', '', 10, ''],
  ['ind.chromLocusPair(10)', '', (1,5), ''],
  ['ind.locusDist( ind.absLocusIndex(1,2) )', '', 3, ''],
  ['len(ind.arrLociDist())', '', 12, ''],
  ['ad[0]','ad=ind.arrLociDist()\nad[0]=1.5', 1.5, ''],
  ['ind.chromBegin(1)', '', 5, ''],
  ['ind.chromEnd(0)', '', 5, ''],
  ['ind.info()', '', 0, ''],
  ['ind.info()', 'ind.setInfo(1)', 1, ''],
  ['ind.affected()', '', 0, ''],
  ['ind.affectedChar()', '', 'U', ''],
  ['ind.affectedChar()', 'ind.setAffected(1)', 'A', ''],
  ['ind.sex()', '', 2, ''],
  ['ind.sexChar()', '', 'F', ''],
  ['ind.sex()', '', Female, ''],
  ['ind.sex()', 'ind.setSex(Male)', Male, ''],
  ['ind.maxAllele()', '', 4, ''],
  ['ind.alleleName(0)', '', '_', ''],
  ['ind.alleleName(9)', '', '_', IndexError],
  ['len(ind.arrAlleles())', '', ind.genoSize(), ''],
  ['len(ind.arrAlleles(0))', '', ind.totNumLoci(), ''],
  ['len(ind.arrAlleles(1,1))', '', ind.numLoci(1), ''],
  ['ind.allele(2, 1, 1)','b=ind.arrAlleles(1,1)\nb[2]=2', 2, ''],
  ['ind.allele(4)', '', 4, ''],
  ['ind.allele(2, 0)', '', 2, ''],
  ['ind.allele(2, 1)', '', 1, '']
]

# perform test
for i in range(0, len(tests)):
  print "Testing ", i
  testExpr( expr=tests[i][0], stmts=tests[i][1],\
    res = tests[i][2], excpt = tests[i][3], d=vars() )
