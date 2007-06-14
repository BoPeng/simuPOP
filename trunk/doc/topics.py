#!/usr/bin/env python
#
# python scripts for topics.lyx
#
from simuPOP import *
#file log/topics_split_merge.py
pop = population(subPop=[100, 200], loci=[1])
pop.splitSubPop(0, [20, 80])
# Note that subpop 1 is intact
print pop.subPopSizes()
pop.splitSubPopByProportion(1, [0.4, 0.6])
print pop.subPopSizes()
# merge
pop.mergeSubPops([1,2])
# Note that subpopulation 2 is not removed
print pop.subPopSizes()
pop.removeEmptySubPops()
print pop.subPopSizes()
pop1 = pop.clone()
pop.mergePopulation(pop1)
print pop.subPopSizes()
#end
#file log/topics_newSubPopSizeFunc.py
def demo(gen, oldsize):
    return [x+10 for x in oldsize]

simu = simulator(
    population(subPop=[100, 200]),
    randomMating(newSubPopSizeFunc=demo)
)
simu.evolve(
    ops = [
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    end=5
)
#end
#file log/topics_split_and_grow.py
def demo(gen, oldsize):
    if gen < 4:
        return [100]
    else:
        return [50+gen]*3

simu = simulator(
    population(size=100),
    randomMating(newSubPopSizeFunc=demo)
)
simu.evolve(
    ops = [
        splitSubPop(which=0, proportions=[0.2, 0.4, 0.4], at=[4]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    end=10
)
#end
#file log/topics_kam_mutation.log
simu = simulator(
    population(size=1000, loci=[1]),
    randomMating()
)
simu.step(ops = [
    kamMutator(rate=0.5, atLoci=[0], maxAllele = 5),
    stat(alleleFreq=[0]),
    pyEval(r'", ".join(["%.2f" % x for x in alleleFreq[0]]) + "\n"')
    ]
)
#end
#file log/topics_smm_mutation.log
simu = simulator(
    population(size=1000, loci=[1]),
    randomMating()
)
simu.evolve(
    preOps = [initByValue([5])],        
    ops = [
        smmMutator(rate=0.1, atLoci=[0], incProb=0.4),
        stat(alleleFreq=[0]),
        pyEval(r'", ".join(["%.2f" % x for x in alleleFreq[0]]) + "\n"')
    ],
    end=5
)
#end
#file log/topics_quantrait.log
pop = population(100, loci=[1, 1], infoFields=['qtrait'])
InitByFreq(pop, [0.4, 0.6])
def qtrait(geno):
    return sum(geno)

PyQuanTrait(pop, loci=[0, 1], func=qtrait)
for i in range(5):
    ind = pop.individual(i)
    print '%d %d %d %d: %.2f' % (ind.allele(0, 0),
        ind.allele(1, 0), ind.allele(0, 1),
        ind.allele(1, 1), ind.info('qtrait'))

#end
