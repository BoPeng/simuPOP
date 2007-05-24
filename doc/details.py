#!/usr/bin/env python
#
# python scripts for details.lyx
#
#file log/details_load.log
from simuOpt import setOptions
setOptions(alleleType='long', optimized=False)
from simuPOP import *
#end
# part 1:
# 
#   Genotypic structure
#   
#file log/details_chrom_stru.log
pop = population(size=10, loci=[2, 4, 5])
print pop.numLoci()
# index starts at zero!
print pop.numLoci(1)
print pop.ploidy()
print pop.ploidyName()
print pop.chromBegin(1)
print pop.locusPos(3)
print pop.locusName(4)
#end
#file log/details_loci_pos.log
pop = population(size=10, loci=[2, 4], maxAllele=3,
    lociPos=[[1.5, 2.5], [1, 2, 5, 10]],
    lociNames=['loc%x' % x for x in range(6)],
    alleleNames=['A', 'T', 'C', 'G'])
print pop.locusPos(3)
print pop.locusName(4)
print pop.alleleName(1)
#end
#file log/details_dump.log
Dump(pop)
#end
#file log/details_init_dump.log
InitByFreq(pop, alleleFreq=[0.3, 0.7])
Dump(pop)
#end
#file log/details_pop_stru.log
pop = population(subPop=[2, 5, 6], loci=[2])
print pop.popSize()
print pop.subPopSizes()
print pop.subPopSize(1)
Dump(pop, infoOnly=True)
#end
#file log/details_subpop_mating.log
pop = population(subPop=[5, 6], loci=[2])
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initByFreq(alleleFreq=[0.2, 0.8], subPop=[0]),
        initByFreq([0, 0, 0, 0.5, 0.5], subPop=[1])
        ],
    ops = [
        dumper(alleleOnly=True, indRange=[[0, 3], [5, 7]]),
        recombinator(rate=0.1) ],
    end = 1
)
#end
