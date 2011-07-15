#!/usr/bin/env python

# This example demonstrates how to use parameter postHook to draw
# a sample and output mutants at each generation after the burn-in period.
#


from simuOpt import setOptions
setOptions(gui=False)
from srv import simuRareVariants, saveMarkerInfoToFile, saveMutantsToFile

import simuPOP as sim
rng = sim.getRNG()

from simuPOP.sampling import drawRandomSample

import sys
if len(sys.argv) == 1:
    name = 'sample'
else:
    name = sys.argv[1]

import logging
logging.basicConfig(filename='%s.log' % name, level=logging.DEBUG)
logger = logging.getLogger(name)

# 170bp  exon
# 5420bp intron  ----5590
# 170bp  exon ...
# 5420bp intron  ----5590x2
# ...      
#                ----5590x8
# 170 bp exon    ----5590x*8 + 170 = 44890
#            (total length 5420*8 + 170*9) = 44890

N=(8000, 8000, 8000, 900000, 1500000)
G=(160000, 40, 370, 40)

def myDist(loc):
    '''Define a function that return fitness parameter of newly arising
    mutants.
    '''
    seg = loc % 5590
    # mutations in intron regions are all neutral
    if seg >= 170 or rng.randUniform() < 0.317:
        return 0
    # purifying selection
    else:
        return rng.randGamma(0.206, 0.146*2)


def drawSample(pop):
    '''
    Draw sample at every generation (after G[0]).
    '''
    if pop.dvars().gen < G[0]:
        return True
    sample = drawRandomSample(pop, 100)
    saveMutantsToFile(sample, '%s_%d.mut' % (name, pop.dvars().gen))
    return True
    
pop = simuRareVariants(regions=['chr1:1..44890'],
    N=N,
    G=G,
    mu=1.2e-8, steps=[1000,10,10,10],
    mutationModel="finite_sites",
    recRate=0.945e-8,
    selModel='multiplicative',
    selDist=myDist,
    selCoef=None,
    postHook=drawSample,
    popFile=('!"%s_%%d.pop" %% gen' % name, (-450, -40, -1)),
    statFile='%s.sts' % name,
    logger=logger
)
saveMarkerInfoToFile(pop, '%s.map' % name)
