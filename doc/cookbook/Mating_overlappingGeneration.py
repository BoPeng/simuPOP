#!/usr/bin/env python
# 
# Test the evolution of LD, with non-random mating schemes.
#

from simuPOP import *
#from rpy import *
import random

import random

size = 1000
maxAge = 4
pop = population(size, loci=[2], infoFields=['age'])
InitByFreq(pop, [0.2, 0.8])
pop.setIndInfo([random.randint(0, maxAge) for x in range(size)], 'age')
pop.setVirtualSplitter(infoSplitter('age', range(maxAge + 1)), 0)

def szFunc(gen, sz=[]):
    return [x+100 for x in sz]

simu = simulator(pop, heteroMating(
    [cloneMating(virtualSubPop=x, weight=-1) for x in (1, 2, 3, 4)] +
    [randomMating(virtualSubPop=3)],
    newSubPopSizeFunc=szFunc)
)
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        infoExec('age += 1', stage=PreMating),
        stat(popSize=True),
        pyEval(r"'Size of age groups: %s\n' % (','.join(['%d' % x for x in virtualPopSize[0]]))",
            stage=PostMating)
        ],
    end = 20
)
   
