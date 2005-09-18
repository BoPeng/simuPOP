# test for simulator

from simuUtil import *

simu = simulator(population(size=10),
  randomMating(), rep=2)

simu.evolve( dryrun=1,
  preOps = [ initByFreq([ .2, .3, .5]) ],
  ops = [
    pyEval("gen"),
    tab(),
    endl(rep=REP_LAST)
  ]
)

# genotypic structure is available also in simulator level.
print simu.ploidy()


# set new size
def newSize(gen,sz):
  os = range(0,len(sz))
  for i in range(0, len(os)):
    os[i] = sz[i] * 1.2
  return os

simu = simulator(population(subPop=[10,20,30]),
  randomMating(newSubPopSizeFunc=newSize))

simu.evolve(
  ops = [
    stat( popSize = 1),
    pyEval('subPopSize'),
    endl()
    ],
  end = 10
  )

# ancestral population

simu = simulator(population(subPop=[10,20,30], ancestralDepth=3),
  randomMating(newSubPopSizeFunc=newSize))

simu.evolve(
  ops = [
    stat( popSize = 1),
    pyEval('subPopSize'),
    endl()
    ],
  end = 10
  )

Dump( simu.population(0), ancestralPops=True, infoOnly=True)


# test numOffsprings
simu = simulator(population(10, loci=[2]), randomMating(numOffsprings=2))
simu.step(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()])

# numOffspringsFunc
def nos(gen):
  return gen+1

simu = simulator(population(10, loci=[2]), randomMating(numOffspringsFunc=nos))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# MATE_NumOffspringsEachFamily
import random
def nos(gen):
  return random.randrange(1,4)

simu = simulator(population(10, loci=[2]), randomMating(numOffspringsFunc=nos))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffspringsFunc=nos, mode=MATE_NumOffspringsEachFamily))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# MATE_GeometricDistribution
simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffsprings=.3, mode=MATE_GeometricDistribution))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# MATE_BinomialDistribution
simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffsprings=.3, maxNumOffsprings=5, mode=MATE_BinomialDistribution))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# MATE_PoissonDistribution
simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffsprings=.3, mode=MATE_PoissonDistribution))
simu.evolve(preOps=[initByFreq([.4,.6])], ops=[parentsTagger(), dumper()], end=4)

# func with distribution
def p(gen):
  return min(gen/10., .9)

turnOnDebug(DBG_MATING)
simu = simulator(
  population(10, loci=[2]),
  randomMating(numOffspringsFunc=p,
    maxNumOffsprings=5,
    mode=MATE_BinomialDistribution))
simu.evolve(preOps=[initByFreq([.4,.6])],
  ops=[parentsTagger()], end=10)

