#!/usr/bin/env python

'''
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

Date: 2008-12-14
'''

from simuPOP import *
import random

def podParentsChooser(pop, sp):
    '''Choose parents of parents from different pods'''
    males = [x for x in pop.individuals() if x.sex() == MALE]
    females = [x for x in pop.individuals() if x.sex() == FEMALE]
    while True:
        # randomly choose a male
        male = males[random.randint(0, len(males)-1)]
        pod = male.pod
        # randomly choose a female from different pod
        while True:
            female = females[random.randint(0, len(females)-1)]
            if female.pod != pod:
                break
        yield (male, female)


def simuPodMating(numPods, podSize):
    '''
    numPods     number of pods
    podSize     size of each pod
    '''
    pop = Population(numPods * podSize, loci=[2], infoFields=['pod'])
    initGenotype(pop, freq=[0.2, 0.8])
    pop.setIndInfo([x for z in range(numPods) for x in [z]*podSize], 'pod')
    pop.setVirtualSplitter(InfoSplitter('pod', range(numPods)))
    pop.evolve(
        initOps = [
            InitSex(),
            InitGenotype(freq=[0.5, 0.5])
        ],
            # calculate size of pods
        preOps = Stat(popSize=True, subPops=[(0,x) for x in range(numPods)]),
            # offspring stays with their natal pod
        matingScheme = HomoMating(
            PyParentsChooser(podParentsChooser),
            OffspringGenerator(numOffspring=1, ops=[
                MendelianGenoTransmitter(),
                InheritTagger(mode=MATERNAL, infoFields='pod')])),
        postOps = [
            # print size of each pod
            PyEval(r"'Size of pods: %s\n' % (','.join(['%d' % x for x in subPopSize]))")
            ],
        gen = 10
    )

if __name__ == '__main__':
    simuPodMating(5, 100)
