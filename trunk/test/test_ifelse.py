# test ifElse operator

from simuUtil import *

simu = simulator(population(10, loci=[2]), randomMating())

simu.evolve(
  preOps = [ initByFreq([.2,.8]) ],
  ops = [
    # count number of allels at this locus
    stat(numOfAlleles=[0]),
    # print it
    pyEval('numOfAlleles'),
    endl(),
    # terminate if fixation happend
    terminateIf('numOfAlleles[0]==1')
  ],
  end=100,
  dryrun=False
)

# now if we want to inject a mutation whenever fixation happens
simu.setGen(0)
simu.evolve(
  preOps = [ initByFreq([.2,.8]) ],
  ops = [
    # count number of allels at this locus
    stat(numOfAlleles=[0]),
    # print it
    pyEval('numOfAlleles'),
    endl(),
    ifElse('numOfAlleles[0]==1', kamMutator(rate=.01, maxAllele=2, atLoci=[0]))
  ],
  end=100,
  dryrun=False
)
