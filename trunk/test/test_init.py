# test initialization

from simuUtil import *

pop = population(subPop=[5,10], loci=[2,4,6])

InitByFreq(pop, [.2, .3, .4, .1], atLoci=[2,4,6])
Dump(pop)
InitByFreq(pop, [.2, .3, .4, .1])
Dump(pop)
InitByFreq(pop, [.2, .3, .4, .1], identicalInds=1)
Dump(pop)
InitByValue(pop, value=[0]*pop.genoSize())
Dump(pop)
InitByFreq(pop, [.2, .8], identicalInds=1, subPop=[0])
Dump(pop)
InitByFreq(pop, [.2, .8], identicalInds=1, indRange=[6,8])
Dump(pop)
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]], identicalInds=1)
Dump(pop)
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]])
Dump(pop)
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]], indRange=[[0,3],[5,13]])
Dump(pop)
InitByValue(pop, value=[0]*pop.genoSize())
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
  indRange=[[0,3],[8,14]], atLoci=[2,3,5,6,9])
Dump(pop)
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
  indRange=[[0,3],[8,14]], atLoci=[2,3,5,6,9],
  identicalInds=1)
Dump(pop)

# init by value
InitByValue(pop, [1]*5 + [2]*7 + [3]*5 +[4]*7)
Dump(pop,alleleOnly=True)
# or one copy of chromosomes
InitByValue(pop, [6]*5 + [7]*7)
Dump(pop,alleleOnly=True)

InitByValue(pop, value=[0]*pop.genoSize())
InitByValue(pop, [1]*5 + [2]*7 + [3]*5 +[4]*7,
  indRange=[[2],[5]])
Dump(pop)
InitByValue(pop, value=[0]*pop.genoSize())
InitByValue(pop, value=[1,2,3,4,5,6], atLoci=[2,4,5],
  indRange=[3,5])
Dump(pop)
#
# by proportion
InitByValue(pop, value=[0]*pop.genoSize())
InitByValue(pop, value=[ [1]*5 +[2]*7, [3]*7+[4]*5 ] ,
  proportions=[.3,.7])
Dump(pop)
InitByValue(pop, value=[0]*pop.genoSize())
# whole ind
InitByValue(pop, value=[[1,2,3,4,5,6],[6,5,4,3,2,1]], atLoci=[2,4,5],
  proportions=[.3,.7],  indRange=[[3,6],[7,10]] )         
Dump(pop)
#
#end

# pyInit
help(pyInit.__init__)
def initAllele(ind, p, sp):
  return sp + ind + p

PyInit(pop, func=initAllele)
Dump(pop, dispWidth=2)
#end


# test atPloidy parameter
pop = population(subPop=[5,10], loci=[2,4,6])

InitByFreq(pop, [.2, .3, .4, .1], atLoci=[2,4,6], atPloidy=0)
Dump(pop)
InitByFreq(pop, [.2, .3, .4, .1], atPloidy=1)
Dump(pop)
InitByFreq(pop, [.2, .3, .4, .1], identicalInds=1, atPloidy=0)
Dump(pop)
InitByValue(pop, value=[0]*pop.genoSize())
Dump(pop)
InitByFreq(pop, [.2, .8], identicalInds=1, subPop=[0], atPloidy=1)
Dump(pop)
InitByFreq(pop, [.2, .8], identicalInds=1, indRange=[6,8], atPloidy=0)
Dump(pop)
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]], identicalInds=1, atPloidy=1)
Dump(pop)
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]], atPloidy=0)
Dump(pop)
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
           indRange=[[0,4],[8,14]], atPloidy=1)
Dump(pop)
InitByValue(pop, value=[0]*pop.genoSize())
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
  indRange=[[2,3],[8,10]], atLoci=[2,3,5,6,9], atPloidy=0)
Dump(pop)
InitByFreq(pop, alleleFreq=[[.2, .8],[.8,.2]],
  indRange=[[0,2],[8,10]], atLoci=[2,3,5,6,9],
  identicalInds=1, atPloidy=1)
Dump(pop)

# init by value


InitByValue(pop, value=[1]*5 + [2]*7 , atPloidy=1)
Dump(pop,alleleOnly=True)
# or one copy of chromosomes
InitByValue(pop, [6]*5 + [7]*7, atPloidy=0)
Dump(pop,alleleOnly=True)

InitByValue(pop, value=[0]*pop.genoSize())
InitByValue(pop, [1]*5 + [2]*7 ,
  indRange=[[2,3],[5,6]], atPloidy=1)
Dump(pop)
InitByValue(pop, value=[0]*pop.genoSize())
InitByValue(pop, value=[1,2,3,4,5,6], atLoci=[2,4,5],
  indRange=[3,5])
Dump(pop)
#
# by proportion
InitByValue(pop, value=[0]*pop.genoSize())
InitByValue(pop, value=[ [1]*5 +[2]*7, [3]*7+[4]*5 ] ,
  proportions=[.3,.7])
Dump(pop)
InitByValue(pop, value=[0]*pop.genoSize())
# whole ind
InitByValue(pop, value=[[1,2,3,4,5,6],[6,5,4,3,2,1]], atLoci=[2,4,5],
  proportions=[.3,.7],  indRange=[[3,6],[7,10]] )         
Dump(pop)
#
#end

# pyInit
help(pyInit.__init__)
def initAllele(ind, p, sp):
  return sp + ind + p

PyInit(pop, func=initAllele)
Dump(pop, dispWidth=2)
#end
