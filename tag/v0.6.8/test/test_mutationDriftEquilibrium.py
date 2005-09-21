# This script test the mutation, drift equilibrium
# for a biallelic case
#

# infinite allele model
import simuOpt
simuOpt.setOptions(longAllele=True)
from simuPOP import *
from simuUtil import *
from simuRPy import *


N, numRep, mu = 5000, 4, 1e-4
simu = simulator(
  population(size=N, ploidy=2, loci=[1]),
  randomMating(),
  rep=numRep)
simu.evolve(
  preOps = [ initByFreq([.1]*10) ],
  ops = [
    #kamMutator(rate=mu, maxAllele=1000),
    smmMutator(rate=mu, maxAllele=1000),
    stat(homoFreq=[0]),
    varPlotter('homoFreq[0]', numRep=numRep, byRep=0, ylim=[0,1],
      update=50, win=2000),
    pyEval("gen, homoFreq[0]", rep=REP_LAST, step=50),
    endl(rep=REP_LAST, step=50),
    pause(stopOnKeyStroke=True),
    ]
  )
