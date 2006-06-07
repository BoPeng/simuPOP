#!/usr/bin/env python

"""

This is a small example of how to use pyOperator to do
custimized mating and condition checking (terminator) 
etc. This is *not* the efficient way to do that since
everything can be handled using a post-mating operator.

Evolution scenario:

A disease spread like this: when one of the chromosomes
of a diploid individual has this allele, it will spread
to the other copy.

  x x x - - - - -   ==>  x x x x x - - -
  x x x x x - - -   ==>  x x x x x - - -

Purpose of simulation:

  Given N, random mating, how many generations to 
  reach all X status?

Note that many simulations will fail because all x
chromosomes disappeared because of genetic drift.
  
"""

import simuOpt, types

options = [
  {'arg': 'h',
   'longarg': 'help',
   'default': False, 
   'description': 'Print this usage message.',
   'allowedTypes': [types.NoneType, type(True)],
   'jump': -1          # if -h is specified, ignore any other parameters.
  },
  {'arg': 'N:',
   'longarg': 'N=',
   'default': 10000,
   'configName': 'Population size',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Population size (10000):  ',
   'description': '''Population size. ''',
   'validate':  simuOpt.valueGT(0)
  },
  {'arg': 'r:',
   'longarg': 'rep=',
   'default': 100,
   'configName': 'Replicates',
   'allowedTypes': [types.IntType, types.LongType],
   'prompt': 'Replicats: (100) ',
   'description': '''Replicates''',
   'validate':  simuOpt.valueGT(0)
  }
]

from simuPOP import *

# this will be used by a duringMating pyOpertor,
# however, for this particular problem, you can 
# use a postMating operator and handle all invidiauls
# one by one. The performance will be *much* better.
#
def offGen(pop, off, dad, mom):
  ''' how to pass allele? 
    We can of course do everything by ourself, but if we do not
    set formOffGenotype of this pyOperator, we can let randomMating()
    generate offspring genotype as usual, we just need to change
    1 2 to 2 2.
  '''
  geno = off.arrGenotype()
  if geno[0] + geno[1] == 3: # in the case of 1 2 or 2 1
    geno[0] = 2
    geno[1] = 2
  return True

def allTwos(pop):
  ''' see if the genotype of the population is all 2'''
  geno = sum(pop.arrGenotype())
  # terminate if all 1 (set flag to fail
  if geno == pop.popSize()*2:   # all 1
    pop.dvars().succ = False
		# return false to terminate simulation
    return False
  # terminate if all 2 (set flag to success)
  elif geno == pop.popSize()*4: # all 2
    pop.dvars().succ = True
		# return false to terminate simulation
    return False
  else:
		# return true to continue simulation
    return True
     
def simu(N):
  ' run the simulation! '
  pop = population(N, loci=[1], ploidy=2)
  InitByValue(pop, value=[1])
  # you can also use a PointerMutator ...
  pop.individual(0).setAllele(2,0)
  pop.individual(0).setAllele(2,1)
  simu = simulator(pop, randomMating())
  simu.evolve(
    ops=[
      pyOperator(stage=DuringMating, func=offGen),
      pyOperator(func=allTwos),
    ]
  )
  return (simu.dvars(0).succ, simu.gen())

if __name__ == '__main__':
  allParam = simuOpt.getParam(options,
    '''This program mimic the evolution of the infection process
  where X- chromosome will be turn to XX. We are concerned  about 
  the speed at which all population becomes XX. ''', __doc__)
  
  if len(allParam) == 0:
    sys.exit(1)
  
  # -h or --help
  if allParam[0]:  
    print simuOpt.usage(options, __doc__)
    sys.exit(0)
  
  # get the parameters
  N = allParam[1]
  rep = allParam[2]

  succCount = 0
  genRecord = []
  for r in range(rep):
    print "Replicate ", r,
    (succ, gen) = simu(N) 
    if succ:
      succCount += 1
      genRecord.append(gen)
      print " success at generation ", gen
    else:
      print " failed"

  print "Population size: ", N
  print "Replicates: ", rep
  print "Successful counts: ", succCount
  if succCount > 0:
    print "Mean generation: ", sum(genRecord)*1.0/succCount
    print "Minimal generation: ", min(genRecord)
    print "Maxmimal generation: ", max(genRecord)
  
