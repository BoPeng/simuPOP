#!/usr/bin/env python
#
# Author: Bo Peng (bpeng@mdanderson.org)
#

'''
This is an example of an evolutionary scenario where the genetic drift of a
specified allele follows a previous simulation.
'''

from simuPOP import *

def recordTrajectory(N, initFreq, gen):
    '''
    Evolve a population gen generations and record its allele
    frequency trajectory.
    '''
    simu = simulator(population(N, loci=[1]), randomMating())
    simu.evolve(
        initOps = [
            initSex(),
            initByFreq([1 - initFreq, initFreq]),
            # initialize an array in the population's local namespace
            pyExec('traj=[]')
        ],
        postOps = [
            stat(alleleFreq=[0]),
            pyExec('traj.append(alleleFreq[0][1])'),
        ],
        gen = gen
    )
    return simu.dvars(0).traj


def simuFollowTrajectory(N, locus, initFreq, traj):
    '''
    Evolve a population, initialize a specified locus with given
    allele frequency and force it to follow a particular allele
    frequency trajectory.
    '''
    # define a trajectory function
    def func(gen):
        return traj[gen]
    #
    simu = simulator(
        population(size=N, loci=[100]), # a larger simulation
        controlledRandomMating(loci=locus, alleles=1, freqFunc=func),
        rep = 5
    )
    simu.evolve(
        initOps = [
            initSex(),
            initByFreq([1 - initFreq, initFreq])
        ],
        postOps = [stat(alleleFreq=[locus], at=-1)],
        gen = len(traj)
    )
    print [simu.dvars(rep).alleleFreq[locus][1] for rep in range(5)]


if __name__ == '__main__':
    traj = recordTrajectory(5000, 0.3, 100)
    print 'Ending allele frequency is:', traj[-1]
    simuFollowTrajectory(5000, 50, 0.3, traj)

