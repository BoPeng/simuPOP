# test variable collector

from simuUtil import *

simu = simulator(population(100, loci=[2]),
  randomMating(), rep=3)

simu.evolve(
  preOps = [initByFreq([.2,.3,.5])],
  ops = [
    stat(alleleFreq=[0,1],
         genoFreq=[0,1],
         LD=[[0,1]]),
    collector(expr='genoFreq[1][1]', name="genoFreq1")
    ],
  end = 10
)
