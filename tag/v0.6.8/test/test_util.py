# test utilities
#

from simuPOP import *
from simuUtil import *

simu = simulator(population(subPop=[100,200], loci=[4,3]),
  randomMating(), rep=2)

count = stat(alleles=range(0,5),
                genotypes=[(1,1,2),(2,1,3)],
                haplotypes=[(1,2,1,2),(1,2,2,3)])

ps = popStat([StatNumOfMale, StatNumOfFemale, StatPopSize])

simu.apply([
  initByFreq([.3,.5,.2]),
  ps, count, hetero(2) ]
)
