#!/usr/bin/env python
# 
# Test the evolution of LD, with non-random mating schemes.
#

from simuPOP import *
#from rpy import *
import random

scheme = 3
#
# mating schemes one: (Amos et al, 1993, Science)
#
# 1. There are several pods, each has around 100 whales
# 2. Male and female in the same pod does not mate
# 3. Male from one pod will mate with female from other
#   pods
#
# Note: This is difficult to achieve using regular
#   can not use subpopulation and migration.
# 

def podParentsChooser(pop, sp):
    '''Choose parents of parents from different pods'''
    males = [x for x in range(pop.popSize()) if pop.individual(x).sex() == Male]
    females = [x for x in range(pop.popSize()) if pop.individual(x).sex() == Female]
    pods = [x.info('pod') for x in pop.individuals()]
    while True:
        # randomly choose a male
        male = males[random.randint(0, len(males)-1)]
        pod = pods[male]
        # randomly choose a female from different pod
        while True:
            female = females[random.randint(0, len(females)-1)]
            if pods[female] != pod:
                break
        yield (male, female)

numPods = 5
podSize = 100
pop = population(numPods * podSize, loci=[2], infoFields=['pod'])
InitByFreq(pop, [0.2, 0.8])
pop.setIndInfo([x for z in range(numPods) for x in [z]*podSize], 'pod')
pop.setVirtualSplitter(infoSplitter('pod', range(numPods)), 0)

simu = simulator(pop, pyMating(
    pyParentsChooser(podParentsChooser),
    mendelianOffspringGenerator(numOffspring=1)))
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        inheritTagger(mode=TAG_Maternal, infoFields=['pod']),
        stat(popSize=True, stage=PrePostMating),
        pyEval(r"'Size of pods: %s\n' % (','.join(['%d' % x for x in virtualPopSize[0]]))")
        ],
    end = 10
)
            

