#!/usr/bin/env python
#
# python scripts for tutorial.lyx
#
from simuOpt import setOptions
setOptions(quiet=True)
#file log/tutorial_example1.log
from simuPOP import *
simu = simulator(
    population(size=1000, ploidy=2, loci=[2]),
    randomMating(),
    rep = 3)
simu.evolve(
    preOps = [initByValue([1,2,2,1])],  
    ops = [
        recombinator(rate=0.1),
        stat(LD=[0,1]),
        pyEval(r"'%3d   ' % gen", rep=0, step=10),
        pyEval(r"'%f    ' % LD[0][1]", step=10),
        pyEval(r"'\n'", rep=REP_LAST, step=10)
    ],
    end=100
)
#end
#file log/tutorial_example2.log
from simuPOP import *
from simuRPy import *
simu = simulator(
    population(size=1000, ploidy=2, loci=[2]),
    randomMating(),
    rep = 3)
simu.evolve(
    preOps = [initByValue([1,2,2,1])],  
    ops = [
        recombinator(rate=0.1),
        stat(LD=[0,1]),
        varPlotter('LD[0][1]', numRep=3, step=10,
                   ylim=[0,.25], xlab='generation',
                   ylab='D', title='LD Decay'),
    ],
    end=100
)
r.dev_print(file='log/LDdecay.eps')
#end
