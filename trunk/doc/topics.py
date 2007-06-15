#!/usr/bin/env python
#
# python scripts for topics.lyx
#
from simuPOP import *
#file log/topics_split_merge.log
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
#file log/topics_newSubPopSizeFunc.log
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
#file log/topics_split_and_grow.log
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
#file log/topics_varying_recombination.log
simu = simulator(
    population(size=1000, loci=[1000]),
    randomMating()
)
rates = [0.00001]*400 + [0.0001]*200 + [0.000001]*399
expr = r'"%.3f %.3f %.3f\n" % (LD[100][101], LD[500][501], LD[900][901])'
simu.evolve(
    preOps = [initByValue([1]*1000+[2]*1000)],
    ops = [
        recombinator(rate=rates, afterLoci=range(999)),
        stat(LD=[[100,101], [500, 501], [900, 901]]),
        pyEval(expr, step=100)
    ],
    end = 1000
)
#end
#file log/topics_recombination_intensity.log
dist = [0.01]*400+[0.1]*200+[0.01]*400
simu = simulator(
    population(size=1000, loci=[1000],
        lociPos=[sum(dist[:x]) for x in range(1000)]),
    randomMating()
)
expr = r'"%.3f %.3f %.3f\n" % (LD[100][101], LD[500][501], LD[900][901])'
simu.evolve(
    preOps = [initByValue([1]*1000+[2]*1000)],
    ops = [
        recombinator(intensity=0.001),
        stat(LD=[[100,101], [500, 501], [900, 901]]),
        pyEval(expr, step=100)
    ],
    end = 1000
)
#end
#file log/topics_migrator.log
simu = simulator(
    population(subPop=[1000]*5),
    randomMating()
)
simu.evolve(
    ops = [
        migrator(rate=[
            [0, 0.2, 0.1],
            [0, 0, 0.1],
            [0.2, 0.2, 0]],
            fromSubPop=[1,2,3], toSubPop=[1,2,3]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    end = 3
)
#end
#file log/topics_penetrance.log
def myPene(geno):
    return (0.01, 0.1, 0.3)[sum(geno)]

simu = simulator(
    population(size=20000, loci=[1]),
    randomMating()
)
expr = r'"%s (%.3f)\n" % (numOfAffected, '
        '1.*numOfAffected/popSize)'
simu.evolve(
    preOps = [initByFreq([0.9, 0.1])],
    ops = [
        pyPenetrance(locus=0, func=myPene),
        stat(numOfAffected=True, popSize=True),
        pyEval(expr, step=10),
    ],
    end=20
)      
#end
#file log/topics_R.log
from rpy import *
def association(pop, loci):
    Stat(pop, alleleFreq=loci)
    a1 = pop.dvars().alleleNum[loci[0]][0]
    a2 = 2*pop.popSize() - a1
    b1 = pop.dvars().alleleNum[loci[1]][0]
    b2 = 2*pop.popSize() - b1
    print '%4d %4d %4d %4d: %.3f' % (a1, a2, b1, b2, 
        r.chisq_test(with_mode(NO_CONVERSION, r.matrix)(
        (a1, a2, b1, b2), ncol=2, byrow=True))['p.value'])
    return True

simu = simulator(
    population(10000, loci=[2]),
    randomMating()
)
simu.evolve(
    preOps=[initByValue([0,1,1,0])],
    ops = [
        recombinator(rate=0.0001),
        pyOperator(func=association, param=[0,1], step=20)
    ],
    end=100
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
#file log/topics_pause.log
simu = simulator(
    population(size=50000, loci=[100]*10),
    randomMating()
)
simu.evolve(
    ops = [
        pause(stopOnKeyStroke=True),
        ticToc(step=10),
    ],
    end = 50
)
#end
#file log/topics_wxPython.log
pop = population(1000, loci=[3,5])
InitByFreq(pop, [.2, .8])
Stat(pop, alleleFreq=range(8), LD=[1,2])
from simuUtil import ListVars
ListVars(pop.vars())
import sys
sys.path.append('../scripts')
from simuViewPop import *
viewPop(pop)
#end
