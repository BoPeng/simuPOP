#!/usr/bin/env python
from simuPOP import *
from simuPOP.utils import MigrIslandRates
from simuPOP.sampling import RandomSample

def calcFst(pop):
    'Calculate Fst and Gst for the whole population and a random sample'
    Stat(pop, structure=range(5), vars=['F_st', 'G_st'])
    sample = RandomSample(pop, size=[500]*pop.numSubPop())[0]
    Stat(sample, structure=range(5), vars=['F_st', 'G_st'])
    print 'Gen: %3d Gst: %.6f (all), %.6f (sample) Fst: %.6f (all) %.6f (sample)' \
        % (pop.dvars().gen,
           pop.dvars().G_st, sample.dvars().G_st,
           pop.dvars().F_st, sample.dvars().F_st)
    return True

simu = simulator(population([10000]*2, loci=[1]*5, infoFields='migrate_to'),
    randomMating())
simu.evolve(
    initOps = [
        initSex(),
        initByFreq([0.5, 0.5], loci=[0, 2]),
        initByFreq([0.2, 0.4, 0.4], loci=[1, 3, 4]),
    ],
    postOps = [
        # migrator(rate=MigrIslandRates(0.01, 3)),
        pyOperator(func=calcFst, step=20),
    ],
    gen = 500
)
