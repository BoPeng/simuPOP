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
#file log/tutorial_population.log
pop = population(size=10, loci=[2, 3])
Dump(pop)
#end
#file log/tutorial_genostru.log
pop = population(subPop=[200, 300], loci=[3, 2], 
    maxAllele=3, ploidy=3, 
    lociPos=[[1, 3, 5], [2.5, 4]], 
    alleleNames=['A', 'C', 'T', 'G'])
pop.numLoci(0)
pop.totNumLoci()
pop.locusPos(4)
pop.subPopSize(1)
pop.popSize()
pop.ploidyName()
pop.individual(1).allele(1, 2)
#end
#file log/tutorial_pop_manipulate.log
# make a copy of pop
pop1 = pop.clone()
# remove loci 2, 3, 4
pop.removeLoci(keep=[0, 1])
# pop2 will have 3 chromosomes, with loci 2, 3, 2
pop2 = MergePopulationsByLoci(pops=[pop, pop1])
# randomly assign alleles using given allele frequencies
InitByFreq(pop2, [0.8, .2])
# calculate population allele frequency
Stat(pop2, alleleFreq=range(pop2.totNumLoci()))
# print allele frequency
print pop2.dvars().alleleFreq
# assign affection status using a penetrance model
MapPenetrance(pop2, locus=1, 
    penetrance={'0-0': 0.05, '0-1': 0.2, '1-1': 0.8})
# draw case control sample
(sample,) = CaseControlSample(pop2, cases=5, controls=5)
# have a look at the sample
Dump(sample, infoOnly=True)
Dump(sample, alleleOnly=True)
#end
