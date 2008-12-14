#!/usr/bin/env python

'''
File: Mating_selfing.py
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates the use of selfing mating schemes.

$Date$
$Revision$
$HeadURL$
'''

from simuPOP import *

def simuSelfing(perc, N, n_rep, gen):
    '''
    perc    percentage of individuals under selfing mating schemes
    N       population size
    n_rep   Number of replicates per simulation
    gen     generations to run
    '''
    pop = population(N, loci=[2])
    pop.setVirtualSplitter(proportionSplitter([perc, 1-perc]))

    simu = simulator(pop,
        heteroMating([
            selfMating(subPop=(0, 0)),
            randomMating(subPop=(0, 1))
        ]),
        rep=n_rep
    )

    simu.evolve(preOps=[
            initByValue([0, 1, 1, 0]),
            pyExec('ld_hist=[]')  # record ld
        ],
        ops=[
            recombinator(rate=0.01),
            stat(LD=[0,1]),
            pyExec('ld_hist.append(LD[0][1])')
        ],
        gen = gen
    )
    print simu.dvars(0).ld_hist
    return 0


if __name__ == '__main__':
    simuSelfing(.4, 1000, 10, 100)
