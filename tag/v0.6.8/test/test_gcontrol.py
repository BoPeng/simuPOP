# testing for Genomic Control

from simuUtil import *
from simuRPy import *
pop = population(subPop=[3000,3000], ploidy=2,
   loci = [5,20,6], maxAllele=2)
simu = simulator(pop, randomMating(), rep=2)

# SNP markers
init = initByFreq([.3,.7])

count = stat(
  genoFreq = [ x for x in range(0,pop.totNumLoci())],
  )

# turnOnDebug(DBG_STATOR)
simu.apply([init, count,
  GC(case=0, control=1, loci=range(0, pop.totNumLoci())  )])
GControl(case=0,control=1,rep=0, loci=range(0,pop.totNumLoci()))

saveInFstatFormat(simu.population(0), "a1.dat")

# the results totally agree with what fstat got
simu.apply([ saveFstat( outputExpr=r"'ran%d.dat' % rep")])

stat(1).listVars(2)

