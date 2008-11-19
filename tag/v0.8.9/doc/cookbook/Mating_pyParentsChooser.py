#!/usr/bin/env python

'''
File: Mating_pyParentsChooser.py
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates the use of a Python parents chooser
  to implement a non-random mating scheme.

  This mating scheme simulates the mating behavior of long-finned
  pilot whales (e.g. Amos et al, 1993, Science)

  1. There are several pods, each has around 100 whales
  2. Male and female in the same pod does not mate
  3. Male from one pod will mate with female from other
     pods

$Date$
$Revision$
$HeadURL$
'''

from simuPOP import *
import random

def podParentsChooser(pop, sp):
    '''Choose parents of parents from different pods'''
    males = [x for x in range(pop.popSize()) if pop.individual(x).sex() == Male]
    females = [x for x in range(pop.popSize()) if pop.individual(x).sex() == Female]
    pods = pop.indInfo('pod', True)
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


def simuPodMating(numPods, podSize):
    '''
    numPods     number of pods
    podSize     size of each pod
    '''
    pop = population(numPods * podSize, loci=[2], infoFields=['pod'])
    InitByFreq(pop, [0.2, 0.8])
    pop.setIndInfo([x for z in range(numPods) for x in [z]*podSize], 'pod')
    pop.setVirtualSplitter(infoSplitter('pod', range(numPods)))

    simu = simulator(pop, pyMating(
        pyParentsChooser(podParentsChooser),
        mendelianOffspringGenerator(numOffspring=1)))
    simu.evolve(
        preOps = [initByFreq([0.5, 0.5])],
        ops = [
            # offspring stays with their natal pod
            inheritTagger(mode=TAG_Maternal, infoFields=['pod']),
            # calculate size of pods
            stat(popSize=True, stage=PrePostMating),
            # print size of each pod
            pyEval(r"'Size of pods: %s\n' % (','.join(['%d' % x for x in virtualPopSize[0]]))")
            ],
        end = 10
    )

if __name__ == '__main__':
    simuPodMating(5, 100)
