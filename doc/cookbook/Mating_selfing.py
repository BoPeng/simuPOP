#!/usr/bin/env python
# 
# Test the evolution of LD, with non-random mating schemes.
#

from simuPOP import *
#from rpy import *
import random

N=10000
n_rep = 10
gen = 100
percentages = [0, 0.2, 0.4, 0.6]

def avg_LD(perc, N, n_rep, gen):
    pop = population(N, loci=[2])
    pop.setVirtualSplitter(proportionSplitter([perc, 1-perc]), 0)
    simu = simulator(pop,
        heteroMating([
            selfMating(subPop=0, virtualSubPop=0),
            randomMating(subPop=0, virtualSubPop=1)
        ]),
        rep=n_rep
    )
    file = 'LD_%.2f.dat' % perc

    simu.evolve(preOps=[
            initByValue([0, 1, 1, 0]),
        ],
        ops=[
            recombinator(rate=0.01),
            stat(LD=[0,1]),
            pyEval('"%.3f" % LD[0][1]', output='>>' + file),
            pyOutput('\t', output='>>' + file),
            pyOutput('\n', rep=REP_LAST, output='>>' + file)
        ],
        end=gen
    )

    ld = open(file)
    avg_ld = []
    for line in ld.readlines():
        avg_ld.append(sum(map(float, line.split()))/n_rep)
    return avg_ld
    
ld_perc = []
if os.path.isfile('ld.res'):
    ld = open('ld.res')
    for line in ld.readlines():
        ld_perc.append(map(float, line.split()))
else:
    ld = open('ld.res', 'w')
    for perc in percentages:
        print 'Percentage of selfing is', perc
        ld_perc.append(avg_LD(perc, N, n_rep, gen))
        print >> ld, ' '.join(['%.4f' % x for x in ld_perc[-1]])
    ld.close()

r.postscript('ld.eps')
r.plot(0, 0, xlim=[0, gen], ylim=[0, 0.25])
for idx,ld in enumerate(ld_perc):
    r.lines(range(0, gen+1), ld)
    r.text(gen/2, ld[gen/2], '%.2f' % percentages[idx])
r.dev_off()
