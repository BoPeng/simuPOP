#!/usr/bin/env python
#
# Purpose:
#  testing of interfaces of population 
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
a = LoadPopulation("pop1.bin")

#pop=population(size=100, ploidy=2, loci=[5, 7], subPop=[20, 80],
#  lociDist=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
#  maxAllele=4, alleleNames=['_','A','C','T','G'])
#InitByFreq(pop, [.2, .8])
#pop.savePopulation("pop1.bin")


# test table:
#
# expr, statements, expected result, expected exception
#
# commented statements will cause memory leak but this is
# expected (very difficult to remove)
#

tests = [
  ['', """population(size=100, ploidy=2, loci=[5, 7], subPop=[20, 80],
  lociDist=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
  maxAllele=4, alleleNames=['_','A','C','T','G'])""", '', ''],
  ['', 'population(size=1000, ploidy=4, loci=[5, 7]*20, subPop=[20, 80]*20)', '', ValueError],
  ['', 'population(subPop=[20,20], ploidy=0)', '', ValueError],
#  ['', 'population(subPop=[2], lociDist=[1,2])', '', TypeError],
   ['', 'population(subPop=[2], loci=[2], lociDist=[[1,2]])', '', ''],
#  ['', 'population(subPop=[2], loci=2, lociDist=[[2,1]])', '', TypeError],
#  ['', 'population(subPop=[2], loci=[2], lociDist=[[2,1]])', '', ValueError],
  ['a.popSize()', '', 100, ''],
  ['a.ploidy()', '', 2, ''],
  ['a.ploidyName()', '', 'diploid', ''],
  ['a.numChrom()', '', 2, ''],
  ['a.numLoci(5)', '', '', IndexError],
  ['a.numLoci(0)', '', 5, ''],
  ['a.numLoci(1)', '', 7, ''],
  ['a.numLoci(2)', '', '', IndexError],
  ['a.locusDist(10)', '', 12, ''],
  ['len(a.arrLociDist())', '', 12, ''],
  ['ad[0]','ad=a.arrLociDist()\nad[0]=1.5', 1.5, ''],
  ['a.chromBegin(1)', '', 5, ''],
  ['a.chromEnd(0)', '', 5, ''],
  ['a.absLocusIndex(1,5)', '', 10, ''],
  ['a.chromLocusPair(10)', '', (1,5), ''],
  ['a.locusDist( a.absLocusIndex(1,2) )', '', 6, ''],
  ['a.numSubPop()', '', 2, ''],
  ['a.totNumLoci()', '', 12, ''],
  ['a.genoSize()', '', a.totNumLoci()*a.ploidy(), ''],
  ['a.absIndIndex(1,1)', '', 21, ''],
  ['a.subPopIndPair(4)', '', (0,4), ''],
  ['a.subPopSize(1)', '', 80, ''],
  ['a.subPopSize(2)', '', '', IndexError],
  ['a.subPopSizes()', '', (20,80), ''],
  ['a.subPopBegin(1)', '', 20, ''],
  ['a.subPopEnd(0)', '', 20, ''],
  ['', 'a.setIndInfoWithSubPopID()', '', ''],
  ['a.individual(20).info()', '', 1, ''],
  ['', 'a.savePopulation("a.txt")', '', ''],
  # judge format from extension
  ['', 'a.savePopulation("a.xml")', '', ''],
  ['', 'a.savePopulation("a.bin", format="bin")', '', ''],
  ['', 'b=LoadPopulation("a.txt")', '', ''],
  ['a.equalTo(b)', '', True, ''],  
  ['', 'b=LoadPopulation("a.xml", format="xml")', '', ''],
  ['a.equalTo(b)', '', True, ''],
  # judge format from extension
  ['', 'b=LoadPopulation("a.bin")', '', ''],
  ['a.equalTo(b)', '', True, ''],
  ['len(a.arrGenotype())', '', a.genoSize()*a.popSize(),''],
  ['len(a.arrGenotype(0))', '', a.genoSize()*a.subPopSize(0),''],
  ['', 'a.individual(0)', '', ''],
  ['', 'a.setIndInfo(range(0,100))', '', ''],
  ['a.popSize()','a.setIndInfo([-1]*50 + [0]*50)\na.setSubPopByIndInfo()', 50, ''],
  ['b.popSize()', 'b=a.newPopByIndInfo()', 50, ''],
  ['a.grp()', '', -1, ''],
  ['a.rep()', '', -1, ''],
  ['', 'a.vars()', '', ''],
  ['a.vars()["b"]', 'a.execute("b=1")', 1, ''],
  ['a.evaluate("b+4")', '', '5', ''],
  ['', 'ind=a.individual(0)','','']
]


# perform test
for i in range(0, min(500,len(tests))):
  print "Testing ", i
  testExpr( expr=tests[i][0], stmts=tests[i][1],\
    res = tests[i][2], excpt = tests[i][3], d=vars() )


# clean up  
import os
os.remove("a.txt")
os.remove("a.bin")
os.remove("a.xml")
# a.bin will be used by other test scripts



# test ancestry history features
pop = population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
InitByFreq(pop, [.2,.8])
Dump(pop, ancestralPops=True)
pop1 = population(subPop=[2,3], loci=[2,3], ancestralDepth=2)
InitByFreq(pop1, [.8,.2])
Dump(pop1)

print pop.ancestralDepth()
pop.pushAndDiscard(pop1)
print pop.ancestralDepth()
Dump(pop1)
Dump(pop, ancestralPops=True)

tmp = pop.clone()
pop.pushAndDiscard(tmp)
print pop.ancestralDepth()
Dump(pop, ancestralPops=True)

# test SavePopulations and LoadPopulations
pop = population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
InitByFreq(pop, [.2,.8])
pop1 = population(subPop=[2,3], loci=[2,3], ancestralDepth=2)
InitByFreq(pop1, [.8,.2])
SavePopulations([pop,pop1], 'pop01.bin')
# pop and pop1 is still alive :-)
Dump(pop)
Dump(pop1)
a = LoadPopulations('pop01.bin')
print len(a)
Dump(a[0])
Dump(a[1])

# test split and merge subpopulations
pop = population(subPop=[5,6,7])
InitByFreq(pop,[.2,.8])
Dump(pop)
SplitSubPop(pop, 1, [2,4],subPopID=[4,1])
Dump(pop)
pop = population(subPop=[5,6,7])
InitByFreq(pop,[.2,.8])
Dump(pop)
SplitSubPop(pop, 2, proportions=[.5,.5])
Dump(pop)
MergeSubPops(pop)
Dump(pop)
SplitSubPop(pop, 0, proportions=[.2,.3,.5])
Dump(pop)
MergeSubPops(pop,[0,2])
Dump(pop)
SplitSubPop(pop, 0, proportions=[.2,.3,.5])
Dump(pop)
MergeSubPops(pop,[2,0])
Dump(pop)
SplitSubPop(pop, 3, proportions=[.5,.5], subPopID=[-1,0])
Dump(pop)

#
#
## # testing serialization of shared vars
## pop = population(10)
## InitByFreq(pop, [.2,.4,.4])
## Stat(pop, alleleFreq=[0])
## d = pop.dvars()
## d.a = 1
## d.b = 1.0
## d.c = [1,2,3.5, "a"]
## d.d = {'1':3,'4':6}
## s =  pop.varsAsString()
## print s
## pop.varsFromString(s)
## listVars(pop.dvars())


#
# test remove loci
pop = population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
pop1 = population(subPop=[2,3], loci=[2,3], ancestralDepth=2)
InitByFreq(pop, [.2,.8])
InitByFreq(pop1, [.5,.5])
pop.pushAndDiscard(pop1)
Dump(pop, ancestralPops=1)
pop.removeLoci(remove=[0,1])
#pop.removeLoci(keep=[3])
Dump(pop, ancestralPops=1)


#
# load population for testing.
#a = LoadPopulation("pop1.bin")
# save and load
#SaveFstat(a,'a.dat')
#p = LoadFstat('a.dat')
#Dump(p)

turnOnDebug(DBG_UTILITY)
# copy of variables
pop =  population(subPop=[3,5], loci=[2,3], ancestralDepth=2)
InitByFreq(pop, [.2,.8])
Stat(pop, LD=[[1,2]])
pop.savePopulation('a.bin')
pop1 = LoadPopulation('a.bin')

pop1 = pop.clone()
pop1.dvars()
pop1.savePopulation('a.bin')
pop.pushAndDiscard(pop1)
pop.dvars().a=1
Dump(pop, ancestralPops=1)
pop.savePopulation('a.bin')

pop2 = pop.clone()
Dump(pop2, ancestralPops=1)
simu = simulator(pop, randomMating())
SavePopulation(simu.getPopulation(0), 'a.bin')
simu.population(0).dvars()
p = simu.getPopulation(0)
p.savePopulation('a.bin')
p.dvars().clear()
