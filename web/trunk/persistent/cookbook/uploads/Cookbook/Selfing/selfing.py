#!/usr/bin/env python

'''
File: Mating_selfing.py
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates the use of selfing mating schemes.

Date: 2008-12-14
'''

from simuPOP import *

def simuSelfing(perc, N, n_rep, gen):
    '''
    perc    percentage of individuals under selfing mating schemes
    N       population size
    n_rep   Number of replicates per simulation
    gen     generations to run
    '''
    pop = Population(N, loci=[2])
    pop.setVirtualSplitter(ProportionSplitter([perc, 1-perc]))

    simu = Simulator(pop, rep=n_rep)

    simu.evolve(
        initOps= [
            InitSex(),
            InitGenotype(genotype=[0, 1, 1, 0]),
            PyExec('ld_hist=[]')  # record ld
        ],
        matingScheme = HeteroMating([
            SelfMating(subPops=[(0, 0)]),
            RandomMating(subPops=[(0, 1)], ops = Recombinator(rates=0.01))
        ]),
        postOps=[
            Stat(LD=[0,1]),
            PyExec('ld_hist.append(LD[0][1])')
        ],
        gen = gen
    )
    print simu.dvars(0).ld_hist
    return 0


if __name__ == '__main__':
    simuSelfing(.4, 1000, 10, 100)
