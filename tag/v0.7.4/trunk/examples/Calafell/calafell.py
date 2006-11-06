#
# aimulation from F. Calcfell, EL Grigorenko, Kidd's paper
#
# Haplotype Evolution and Linkage Disequlibrium: A simulation
# study

from simuPOP  import *
from simuUtil import *
from simuRPy  import *
import math, random

# loci:
#   B1    biallelic RFLP    B1 B2
#   (CA)n STRP              CA11 - CA17
#   A1    biallelic RFLP:   A1 A2

#
# recombination
#   B1 - CA:  0.006%
#   CA - A1:  0.020%

recombine = recombinator(
  rates = [0.00006, 0.0002],
  afterLoci=[0,1] )
    
#
# Mutation:
#
#   CA  symmetric stepwise mutation model
#     mutation exceeds boundaries are discaded. (e.g. 11->10, 17->18))
#     rate of mutate one:  1e-4
#     mutate two steps: 2e-5
#     mutate three steps: 4e-6
#     ... at 1/5 of previous rate
#
#   A1 e-8
#   B1 e-8

# with r = 1e-4, call this func
def step():
  return int(-math.log(random.uniform(0,1))/math.log(5))+1

# test this step machanism:
# with 5 times less likely to go a step further
#s = [0] * 10000
#for i in range(0,10000):
#  s[i] = step()
## table:
#for i in range(0,10):
#  print s.count(i)

# mutators
mutate02 = kamMutator( atLoci=[0,2], rate=1e-8, maxAllele=2)
mutate1 =  gsmMutator( atLoci=[1], maxAllele=7, rate=1e-4, func = step)

#
# initial generation
#   b1-ca13-a1   40%
#   b2-ca14-a2   60%
#
init = initByValue(
  values = [ [ 1, 3, 1], [2, 4, 2] ],  # two genotypes
  proportions = [.4, .6 ]              # by proportion
)

#
# measurements:
#   model 1,2,3 at 1000, 2500 and 5000 gen
#   model 4, 5 at 5000 gen
#
# Meaaure: D' (-1 ~ 1 )
#   between A and B
#   between A and CA11-CA17
#   between B and CA11-CA17

LD_Freq = stat(
  LD= [[0,2]] + [[0,1,1,x] for x in range(0,8)]
    + [[1,2,x,1] for x in range(0,8)],
  alleleFreq=[0,1,2],
  step=10)   # calculate statistics every 10 steps

# demographics scenarios:
#
# 1. 2N=2000, 5000 gen
# 2. 2N=4000, 5000 gen
# 3. 2N=10000,5000 gen
# 4. 2N=2848 for 2786 gen, 2N=10,000 aferwards
# 5. 2N=3600 for 3400 gen, reaching 2N=10,000 at 5000 gen

# the following functions defines how population size
# will change with generation. The return values
# should be an array of subpopulation sizes. In this
# case, we should return an array of size one.
def scenario_1(gen,sz):
  return [2000]

def scenario_2(gen,sz):
  return [5000]

def scenario_3(gen,sz):
  return [10000]

def scenario_4(gen,sz):
  if gen < 2786:
    return 2848
  elif gen < 3000: # rapid linear growth
    return [(10000-2848)/(3000-2786)*(gen-2786)+2848]
  else:
    return [10000]

def scenario_5(gen, sz):
  if gen < 3400:
    return [3600]
  else:
    return [(10000-3600)/(5000-3400)*(gen-3400)+3600]


####################### 

def simulation(nRep, scenario, endGen, visualizers=[]):
  # one chromosome with three loci
  pop = population( subPop=scenario(0), ploidy=2, loci=[3])
    
  # simulator with random mating, scenario is a function accepting
  # generation as its argument and return subpopsize...
  simu = simulator(pop, randomMating(newSubPopSizeFunc=scenario),
    rep=nRep)
  # evolve
  simu.evolve(
    preOps = [ init ],
    ops = [
      LD_Freq,   # check LD' value and allele frequency
      recombine,
      mutate02,  # mutate loci 0 and 2
      mutate1,   # mutate locus 1
      # plot LD' between B1 and A1 loci
      varPlotter("LD_prime['0-2|1-1']", numRep=nRep,
        step=10, update=10, saveAs="LDprime"),
      # plot allele frequencies of B1 A1 
      varPlotter("[alleleFreq[0][1], alleleFreq[2][1]]",
        numRep=nRep, byRep=True, step=10, update=10, varDim=2,
        ylim=[0,1],saveAs="alleleFreq02"),
      # plot allele frequencies at CA locus
      varPlotter("alleleFreq[1][1:]",
        numRep=nRep, byVal=True, step=10, update=10, varDim=7,
        ylim=[0,1], saveAs="AlleleFreq1"),
      # keep track of the following values at different
      # generation
      collector( expr='alleleFreq[1]', name = 'CA_af', 
        at=[500,1000,1500,2000,2500,3000,3500,4000,4500,5000]),
      # report progress
      pyEval(r'"%d\n" % gen',rep=REP_LAST)
      ] + visualizers ,   
    end = endGen,
    dryrun=False
  )
  return simu
 
# run!
#simu = simulation(100, scenario_3, 5000)

#simu.saveSimulator("sec3_99_5000.dat", format="bin")

import math

def analyze(simu, nrep, gen):
  d = []
  mean = [0]*7
  std = [0]*7
  for i in range(0, nrep):
    d.append(simu.dvars(i))
    # calulate mean frequency of CA 11 - 17
  for a in range(0,7):
    fq = [0]*nrep
    sum1 = 0.
    sum2 = 0.
    for i in range(0, 100):
      fq[i] = d[i].CA_af[gen][a+1]
      sum1 += fq[i]
      sum2 += fq[i]*fq[i]
    
    mean[a] = sum1/100.
    std[a] = math.sqrt((sum2 - nrep*mean[a]*mean[a])/(nrep-1))
  return {'mean':mean, 'std':std}

#listVars(analyze(simu, 100, 1000))
#listVars(analyze(simu, 100, 2500))
#listVars(analyze(simu, 100, 5000))

simu = simulation(4, scenario_3, 1000)
