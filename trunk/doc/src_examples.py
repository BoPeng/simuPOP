#
# This script will genreate code/result pieces
# that will be inserted into simuPOP user's guide
# and reference manual
# 
# #file_begin file
#
# #file_end file
#
# will be used in this file so the running result can be
# separated into files specified by 'file'
#
# #PS after commands that will be executed.
#
# create directory log if not exist
import os, sys

if not os.path.isdir('log'):
  try:
    os.mkdir('log')
  except:
    print "Failed to make output directory log"
    sys.exit(1)

#
from simuPOP import *


#file log/src_genoStruTrait.log
# create a population, most parameters have default values
pop = population(size=5, ploidy=2, loci=[5,10],
    lociPos=[range(0,5),range(0,20,2)],
    alleleNames=['A','C','T','G'],
    subPop=[2,3], maxAllele=3)
print pop.popSize()
print pop.ploidy()
print pop.ploidyName()
print pop.numChrom()
print pop.locusPos(2)
print pop.alleleName(1)
# get the fourth individual of the population
ind = pop.individual(3)
# access genotypic structure info
print ind.ploidy()
print ind.numChrom()
print ind.numLoci(0)
print ind.genoSize()
# and from simulator level
simu = simulator(pop, randomMating(), rep=3)
print simu.numChrom()
#end

#file log/src_population.log
# use of population function
# a Wright-Fisher population
WF = population(size=100, ploidy=1, loci=[1])

# a diploid population of size 10
# there are two chromosomes with 5 and 7 loci respectively
pop = population(size=10, ploidy=2, loci=[5, 7], subPop=[2, 8])

# a population with SNP markers (with names A,C,T,G)
# range() are python functions
pop = population(size=5, ploidy=2, loci=[5,10],
    lociPos=[range(0,5),range(0,20,2)],
    alleleNames=['A','C','T','G'],
    subPop=[2,3], maxAllele=3)

#
# population structure functions
print pop.popSize()
print pop.numSubPop()
print pop.subPopSize(0)
print pop.subPopSizes()
print pop.subPopBegin(1)
print pop.subPopEnd(1)
print pop.subPopIndPair(3)
print pop.absIndIndex(1,1)

#
# functions of setting population structure
pop.setIndSubPopID([1,2,2,3,1])
pop.setSubPopByIndID()
pop.removeLoci(keep=range(2,7))
Dump(pop)

#
# save and load population
# save it in various formats, default format is "txt"
pop = population(1000, loci=[2, 5, 10])
pop.savePopulation("pop.txt")
pop.savePopulation("pop.txt", compress=False)
pop.savePopulation("pop.xml", format="xml")
pop.savePopulation("pop.bin", format="bin")

# load it in another population
pop1 = LoadPopulation("pop.xml", format="xml")
#end


#file log/src_individual.log
pop = population(500, loci=[2, 5, 10])
# get an individual
ind = pop.individual(9)
# oops, wrong index
ind = pop.individual(3)
# you can access genotypic structure info
print ind.ploidy()
print ind.numChrom()
# ...
# as well as genotype
print ind.allele(1) 
ind.setAllele(1,5)
print ind.allele(1)
# you can also use an overloaded function
# with a second parameter being the ploidy index
print ind.allele(1,1) # second locus at the second copy of chromosome
# other information
print ind.affected()
print ind.affectedChar()
ind.setAffected(1)
print ind.affectedChar()
print ind.sexChar()
#end


#file log/src_initByFreq.log
simu = simulator( 
    population(subPop=[2,3], loci=[5,7], maxAllele=1),
    randomMating(), rep=1)
simu.step([
    initByFreq(alleleFreq=[ [.2,.8],[.8,.2]]),
    dumper(alleleOnly=True)
  ])
#end

#file log/src_initByValue.log
simu = simulator(
    population(subPop=[2,3], loci=[5,7], maxAllele=9),
    randomMating(), rep=1)
simu.step([
    initByValue([1]*5 + [2]*7 + [3]*5 +[4]*7),
    dumper(alleleOnly=True)])
#end

#file log/src_pyInit.log
def initAllele(ind, p, sp):
  return sp + ind + p

simu = simulator( 
    population(subPop=[2,3], loci=[5,7]),
    randomMating(), rep=1)
simu.step([
    pyInit(func=initAllele),
    dumper(alleleOnly=True, dispWidth=2)])
#end


#file log/src_kamMutator.log
simu = simulator(population(size=5, loci=[3,5]), noMating())
simu.step([
    kamMutator( rate=[.2,.6,.5], loci=[0,2,6], maxAllele=9),
    dumper(alleleOnly=True)])
#end

#file log/src_smmMutator.log
simu = simulator(population(size=3, loci=[3,5]), noMating())
simu.step([
    initByFreq( [.2,.3,.5]),
    smmMutator(rate=1,  incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])
#end

#file log/src_gsmMutator.log
simu = simulator(population(size=3, loci=[3,5]), noMating())
simu.step([
    initByFreq( [.2,.3,.5]),
    gsmMutator(rate=1, p=.8, incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])

import random
def rndInt():
  return random.randrange(3,6)

simu.step([
    initByFreq( [.2,.3,.5]),
    gsmMutator(rate=1, func=rndInt, incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])

#end

#file log/src_pyMutator.log
def mut(x):
  return 8

simu = simulator(population(size=3, loci=[3,5]), noMating())
simu.step([
  pyMutator(rate=.5, loci=[3,4,5], func=mut),
  dumper(alleleOnly=True)])
#end


#file log/src_recombinator.log
simu = simulator(
    population(4, loci=[4,5,6], maxAllele=9,
    infoFields=['father_idx', 'mother_idx']),
    randomMating())
simu.step([
  parentsTagger(),
  ],
  preOps = [initByFreq([.2,.2,.4,.2]), dumper(alleleOnly=True) ],
  postOps = [ dumper(alleleOnly=True)]
)
simu.step([
  parentsTagger(),
  recombinator(rate=[1,1,1], afterLoci=[2,6,10])
  ],
  postOps = [ dumper(alleleOnly=True)]
)
#end


#file log/src_pyOperator.log
def dynaMutator(pop, param):
  ''' this mutator mutate common loci with low mutation rate
  and rare loci with high mutation rate, as an attempt to
  bring allele frequency of these loci at an equal level.'''
  # unpack parameter
  (cutoff, mu1, mu2) = param;
  Stat(pop, alleleFreq=range( pop.totNumLoci() ) )
  for i in range( pop.totNumLoci() ):
    # 1-freq of wild type = total disease allele frequency
    if 1-pop.dvars().alleleFreq[i][1] < cutoff:
      KamMutate(pop, maxAllele=2, rate=mu1, loci=[i])
    else:
      KamMutate(pop, maxAllele=2, rate=mu2, loci=[i])
  return True

pop = population(size=10000, ploidy=2, loci=[2, 3])

simu = simulator(pop, randomMating())

simu.evolve(
  preOps = [ 
    initByFreq( [.6, .4], loci=[0,2,4]),
    initByFreq( [.8, .2], loci=[1,3]) ],
  ops = [ 
    pyOperator( func=dynaMutator, param=(.5, .1, 0) ),
    stat(alleleFreq=range(5)),
    pyEval(r'"%f\t%f\n"%(alleleFreq[0][1],alleleFreq[1][1])', step=10)
    ],
  end = 30
)        
#end



#file log/src_mating.log
# arbitrary demographic model
def lin_inc(gen, oldsize=[]):
    return [10+gen]*5

simu = simulator(
    population(subPop=[5]*5, loci=[1]),
    randomMating(newSubPopSizeFunc=lin_inc)
)
simu.evolve(
    ops = [
        stat(popSize=True),
        pyEval(r'"%d %d\n"%(gen, subPop[0]["popSize"])'),
    ],
    end=5
)

#
# control the number of offspring per mating event
# famSizes is only defined when DBG_MATING is defined
TurnOnDebug(DBG_MATING)
simu = simulator(population(50, loci=[1]),
    randomMating(numOffspring=2, 
        maxNumOffspring=5,
        mode=MATE_UniformDistribution))
simu.step(ops=[])
print simu.population(0).dvars().famSizes
TurnOffDebug(DBG_MATING)
#end


#file log/src_mapSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
    randomMating())
s1 = .1
s2 = .2
simu.evolve([
    stat( alleleFreq=[0], genoFreq=[0]),
    mapSelector(locus=0, fitness={'0-0':(1-s1), '0-1':1, '1-1':(1-s2)}),
    pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    preOps=[  initByFreq(alleleFreq=[.2,.8])],
    end=300)
#end


#file log/src_maSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
    randomMating())
s1 = .1
s2 = .2
simu.evolve(
    preOps=[initByFreq(alleleFreq=[.2,.8])],
    ops = [
        stat( alleleFreq=[0], genoFreq=[0]),
        maSelector(locus=0, fitness=[1-s1, 1, 1-s2]),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    end=300)
#end

#file log/src_pySelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[3], infoFields=['fitness'] ),
    randomMating())

s1 = .2
s2 = .3
# the second parameter gen can be used for varying selection pressure
def sel(arr, gen=0):
  if arr[0] == 1 and arr[1] == 1:
    return 1 - s1
  elif arr[0] == 1 and arr[1] == 2:
    return 1
  elif arr[0] == 2 and arr[1] == 1:
    return 1
  else:
    return 1 - s2

# test func
print sel([1,1])

simu.evolve([
    stat( alleleFreq=[0], genoFreq=[0]),
    pySelector(loci=[0,1],func=sel),
    pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=25)
    ],
    preOps=[  initByFreq(alleleFreq=[.2,.8])],
    end=100)
#end

#file log/src_pySubset.log
simu = simulator(population(subPop=[2,3], loci=[3,4], infoFields=['fitness']),
    randomMating())
simu.step([
    initByFreq([.3,.5,.2]),
    pySubset( [1,-1,-1,1,-1] ),
    dumper(alleleOnly=True, stage=PrePostMating)
   ])
#end


#file log/src_mlPenetrance.log
pop = population(1000, loci=[3])
InitByFreq(pop, [0.3, 0.7])
pen = []
for loc in (0, 1, 2):
    pen.append( maPenetrance(locus=loc, wildtype=[1],
        penetrance=[0, 0.3, 0.6] ) )

# the multi-loci penetrance
MlPenetrance(pop, mode=PEN_Multiplicative, peneOps=pen)
Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#end


#file log/src_pyPenetrance.log
pop = population(1000, loci=[3])
InitByFreq(pop, [0.3, 0.7])
def peneFunc(geno):
    p = 1
    for l in range(len(geno)/2):
        p *= (geno[l*2]+geno[l*2+1])*0.3
    return p

PyPenetrance(pop, func=peneFunc, loci=(0, 1, 2))
Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#
# You can also define a function, that returns a penetrance
# function using given parameters
def peneFunc(table):
    def func(geno):
      return table[geno[0]][geno[1]]
    return func

# then, given a table, you can do
PyPenetrance(pop, loci=(0, 1, 2),
    func=peneFunc( ((0,0.5),(0.3,0.8)) ) )
#end



#turnOnDebug(DBG_ALL)
#turnOnDebug(DBG_SIMULATOR)

#file log/src_ifElse.log
from simuRPy import *
from simuUtil import *
numRep=4
popSize=100
endGen=50

simu = simulator(population(size=popSize, loci=[1]),
  randomMating(), rep=numRep)
simu.evolve(
  preOps = [ initByValue([1,1])],
  ops = [
    # penetrance, additve penetrance
    maPenetrance(locus=0, wildtype=[1], penetrance=[0,0.5,1]),
    # count number of affected
    stat(numOfAffected=True),
    # introduce disease if no one is affected
    ifElse(cond='numOfAffected==0',
      ifOp=kamMutator(rate=0.01, maxAllele=2)),
    # expose affected status
    pyExec('pop.exposeAffectedness()', exposePop=True),
    # plot affected status
    varPlotter(expr='affected',plotType="image", byRep=1, update=endGen, 
      varDim=popSize, win=endGen, numRep=numRep,
      title='affected status', saveAs="ifElse")
  ],
  end=endGen,
  dryrun=False
)
#end

#file log/src_none.log
# this may be set from command line option
savePop = False
# then, saveOp is defined accordingly
if savePop:
    saveOp = savePopulation(output='a.txt')
else:
    saveOp = noneOp()
simu = simulator(population(10), randomMating())
simu.evolve( [saveOp])
#end


#file log/src_rng.log
print ListAllRNG()
print rng().name()
SetRNG("taus2", seed=10)
print rng().name()
#end

#file log/src_rngrand.log
r=rng()
#help(RNG)
for n in range(1,10):
  print r.randBinomial(10, .7),
#end




