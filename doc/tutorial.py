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
        varPlotter('LD[0][1]', numRep=3, step=10, saveAs='ld',
            ylim=[0,.25], lty=range(1, 4), col=range(2, 5),
            xlab='generation', ylab='D', title='LD Decay'),
    ],
    end=100
)
#end
#PS epstopdf ld10.eps
#PS epstopdf ld20.eps
#PS epstopdf ld30.eps
#PS epstopdf ld40.eps
#PS epstopdf ld50.eps
#PS epstopdf ld60.eps
#PS epstopdf ld70.eps
#PS epstopdf ld80.eps
#PS epstopdf ld90.eps
#PS epstopdf ld100.eps
