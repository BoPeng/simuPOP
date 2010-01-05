#!/usr/bin/env python

# A spatial migration example.
#
# Suppose subpopulations are located at (x_i,y_i), migration rate between
# two subpopulations is exp(-r*d_ij) where d_ij is the Euclidean distances
# between two subpopulations and r= 1/mean dispersal.
#

from simuPOP import *
import math

def spatialMigrRates(xy, r):
    '''
    Return a migration matrix where migration rates between two
    subpopulations vary according to Euclidean distance between them.

    xy
        A list of (x,y) location for each subpopulation.

    r
        migrate rate between two subpopulations is exp(-r*d_ij) where
        d_ij is the Euclidean distance between subpopulations i and j.
    '''
    nSubPop = len(xy)
    rate = []
    for i in range(nSubPop):
        rate.append([])
        for j in range(nSubPop):
            if i == j:
                rate[-1].append(0)
                continue
            d_ij = math.sqrt((xy[i][0] - xy[j][0])**2 + (xy[i][1] - xy[j][1])**2)
            rate[-1].append(math.exp(-1 * r * d_ij))
    return rate

def printAlleleFreq(pop):
    'Print allele frequencies of all subpopulations'
    stat(pop, alleleFreq=[0], vars=['alleleFreq_sp'])
    print 'Allele frequencies at generation', pop.dvars().gen
    for i in range(10):
        for j in range(10):
            print '%.2f' % pop.dvars(10*i + j).alleleFreq[0][1],
        print
    return True


def simuSpatial():
    '''
    A example.
    '''
    xy = []
    for i in range(10):
        for j in range(10):
            xy.append((i, j))
    r = spatialMigrRates(xy, 3)
    pop = Population(size=[100]*100, loci=[1],
        infoFields='migrate_to')
    pop.evolve(
        # only subpopulation 55 has genotype 1, 1
        initOps = [
            InitSex(),
            InitGenotype(genotype=[1, 1], subPops=55),
        ],
        preOps = Migrator(rate=r),
        matingScheme = RandomSelection(),
        postOps = PyOperator(printAlleleFreq, at=3),
        gen = 10
    )
                

if __name__ == '__main__':
    simuSpatial()
