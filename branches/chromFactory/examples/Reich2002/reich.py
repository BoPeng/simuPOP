#!/usr/bin/env python
# simulation for Reich's paper:
#   On the allelic spectrum of human disease
#
import simuOpt
# optimized should be false during testing
simuOpt.setOptions(optimized=True)
from simuPOP import *
from simuRPy import *    # use R for plotting

initSize =  10000        # initial population size
finalSize = 1000000      # final population size
burnin = 500             # evolve with constant population size
endGen = 1000            # last generation
mu = 3.2e-5              # mutation rate
C_f0 = 0.2               # initial allelic frequency of *c*ommon disease
R_f0 = 0.001             # initial allelic frequency of *r*are disease
max_allele = 255         # allele range 1-255 (1 for wildtype)
C_s = 0.0001             # selection on common disease
R_s = 0.9                # selection on rare disease
psName = 'lin_exp'       # filename of saved figures 

C_f = [1-C_f0] + [x*C_f0 for x in [0.9, 0.02, 0.02, 0.02, 0.02, 0.02]]
R_f = [1-R_f0] + [x*R_f0 for x in [0.9, 0.02, 0.02, 0.02, 0.02, 0.02]]

## five different population change scenarios
#
# constant small population size
def c_small(gen, oldSize=[]):
  return [initSize]

# constant large population size
def c_large(gen, oldSize=[]):
  return [finalSize]

# linear growth after burn-in
def lin_exp(gen, oldSize=[]):
  if gen < burnin:
    return [initSize]
  elif gen % 10 != 0:
    return oldSize
  else:
    incSize = (finalSize-initSize)/(endGen-burnin)
    return [oldSize[0]+10*incSize]

# instantaneous population growth
def ins_exp(gen, oldSize=[]):
  if gen < burnin:
    return [initSize]
  else:
    return [finalSize]

# expoenential growth after burn-in
def exp_exp(gen, oldSize=[]):
  if gen < burnin:
    return [initSize]
  elif gen%10 != 0:
    return oldSize
  else:
    incRate = (math.log(finalSize)-math.log(initSize))/((endGen-burnin)/10)
    return  [int(initSize*math.exp((gen-burnin)/10.*incRate))]

def simulate(incSenario):
  simu = simulator(                                  # create a simulator
    population(subPop=incSenario(0), loci=[1,1]),    # inital population
    randomMating(newSubPopSizeFunc=incSenario)       # random mating
  )
  simu.evolve(              # start evolution
    preOps=[                # operators that will be applied before evolution
      # initialize locus 0 (for common disease)
      initByFreq(atLoci=[0], alleleFreq=C_f),
      # initialize locus 1 (for rare disease)
      initByFreq(atLoci=[1], alleleFreq=R_f),
    ],
    ops=[                   # operators that will be applied at each gen
      # report population size, for monitoring purpose only
      # count allele frequencies at both loci
      stat(popSize=True, alleleFreq=[0,1]),
      # report generation and popsize
      pyEval(r"'%d\t%d\n' % (gen, popSize)", step=5),
      # mutate: k-alleles mutation model
      kamMutator(rate=mu, maxAllele=max_allele),
      # selection on common and rare disease,
      mlSelector([        # multiple loci - multiplicative model
        maSelector(locus=0, fitness=[1,1,1-C_s], wildtype=[1]),
        maSelector(locus=1, fitness=[1,1,1-R_s], wildtype=[1])
      ], mode=SEL_Multiplicative),   
      # visualization of allelic frequencies/spectra
      # use a freqPlotter defined in simuRPy.py
      pyEval(
        # run once when this operator is created
        preStmts='''p=freqPlotter(max=max_allele, y1max=C_f0*1.5, y2max=R_f0*1.5,
          x1lab="common disease alleles", x2lab="rare disease alleles", 
          save=10, name=psName)''',
        # plot allele frequencies every 10 generations
        stmts='p.plot(gen, popSize, alleleFreq[0], alleleFreq[1])',
          step=10),
      # monitor time
      ticToc(step=5),
      # pause at any user key input (for presentation purpose)
      pause(stopOnKeyStroke=1)
    ],
    end=endGen
  )

#simulate(c_small)
#simulate(c_large)
simulate(lin_exp)
#simulate(ins_exp)
#simulate(exp_exp)

