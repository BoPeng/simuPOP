#!/usr/bin/env python
#
# Author: Bo Peng (bpeng@mdanderson.org)
#

'''
This is an example of an evolutionary scenario where a simulation is restarted
if certain condition is met.
'''

from simuPOP import *

def simuRestartIfFailed(N, initFreq, freqRange, gen, genCheck):
    '''
    Evolve a population, restart the simulation if the allele frequency
    is not between freqRange at 'genCheck'.

    N:         population size
    initFreq:  initial frequency of allele 1
    freqRange: range of allele frequency
    gen:       total generations to evolve
    genCheck:  when to check allele frequency
    '''
    pop = Population(size=N, loci=[1])
    initSex(pop)
    initGenotype(pop, freq=[1 - initFreq, initFreq])
    # put freqRange as a population variable so that it can be
    # used in an expression
    pop.dvars().fr = freqRange
    while True:
        simu = Simulator(pop, steal=False)
        evolved = simu.evolve(
            matingScheme = RandomMating(),
            postOps = [
                Stat(alleleFreq=[0], at=genCheck),
                TerminateIf('alleleFreq[0][1] > fr[1] or alleleFreq[0][1] < fr[0]',
                    at=genCheck)
            ],
            gen = gen
        )
        if evolved[0] != gen:
            print 'Restart simulation due to allele frequency %.3f.' % \
                simu.dvars(0).alleleFreq[0][1]
        else:
            print 'Allele frequency %.3f is within specified range.' % \
                simu.dvars(0).alleleFreq[0][1]
            break

if __name__ == '__main__':
    simuRestartIfFailed(1000, 0.3, [0.35, 0.40], 100, 50)

