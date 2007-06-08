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
#file log/details_population.log
pop = population(size=10, loci=[2, 3])
Dump(pop)
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
#file log/details_infofields.log
pop = population(size=20, loci=[1], infoFields=['birthday'])
InitByFreq(pop, alleleFreq=[0.2, .8])
#end
#file log/details_genostru.log
pop = population(subPop=[200, 300], loci=[3, 2], 
    maxAllele=3, ploidy=4, 
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
#file log/details_pop_manipulate.log
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
# save sample in Merlin QTDT format
from simuUtil import SaveQTDT
SaveQTDT(sample, output='sample', affectionCode=['U', 'A'], 
    fields=['affection'])
# have a look at the sample in Merlin-QTDT Format
print open('sample.map').read()
print open('sample.dat').read()
print open('sample.ped').read()
#end
#file log/details_pop_variable.log
pop = population(subPop=[5, 10], loci=[5])
InitByFreq(pop, [.6, .3, .1])
Stat(pop, alleleFreq=[1], genoFreq=[2])
print pop.dvars().alleleFreq[1][0]
from simuUtil import ListVars
ListVars(pop.dvars(), useWxPython=False)
#end
#file log/details_individual.log
pop = population(subPop=[5, 8], loci=[5], 
    infoFields=['penetrance'])
InitByFreq(pop, [.6, .3, .1])
MaPenetrance(pop, locus=2, penetrance=[0.05, 0.2, 0.5],
    wildtype=[0], infoFields=['penetrance'])
# iterate through all inviduals in subPop 1
for ind in pop.individuals(1):
    print 'Aff: %d Fit: %.3f Geno: %d %d' % \
        (ind.affected(), ind.info('penetrance'), \
        ind.allele(2, 0), ind.allele(2, 1))

#end
#file log/details_info.log
pop = population(100, loci=[5, 8],
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, randomMating(numOffspring=2))
simu.evolve(ops=[parentsTagger()], end=5)
ind = simu.population(0).individual(0)
ind1 = simu.population(0).individual(1)
print ind.info('father_idx'), ind.info('mother_idx')
print ind1.info('father_idx'), ind1.info('mother_idx')

#end
#file log/details_op_stage.log
simu = simulator(
    population(subPop=[20, 80], loci=[3]),
    randomMating())
simu.evolve(
    preOps = [initByFreq([0.2, 0.8])],
    ops = [
        kamMutator(maxAllele=10, rate=0.00005, atLoci=[0,2]),
        recombinator(rate=0.001),
        dumper(stage=PrePostMating),
        stat(alleleFreq=[1]),
    ],
    dryrun=True
)
#end
#file log/details_op_gen.log
simu = simulator(
    population(10000, loci=[3]),
    randomMating())
eval1 = r"'Gen: %3d  Freq: %f\n' % (gen, alleleFreq[1][0])"
eval2 = r"'Last Gen: %3d  Freq: %s\n' % (gen, alleleFreq[1])"
simu.evolve(
    preOps = [initByFreq([0.3, 0.7])],
    ops = [
        recombinator(rate=0.01, begin=10, end=30),
        stat(alleleFreq=[1], step=10),
        pyEval(eval1, step=10),
        pyEval(eval2, at=[-1])
    ],
    end = 50
)
#end
#file log/details_op_rep.log
simu = simulator(
    population(100, loci=[3]),
    randomMating(), 
    rep=5, grp=[1,1,2,2,2])
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        stat(alleleFreq=[1]),
        recombinator(rate=0.01, grp=1),
        recombinator(rate=0.01, grp=2),
        pyEval(r"'%.2f ' % alleleFreq[1][0]", grp=1),
        pyEval(r"'\n'", rep=REP_LAST),
    ],
    end=5
)
#end
#file log/details_op_output.log
simu = simulator(
    population(100, loci=[3]),
    randomMating(), 
    rep=5, grp=[1,1,2,2,2])
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        stat(alleleFreq=[1]),
        pyEval(r"'%.2f ' % alleleFreq[1][0]", 
            output='>>out'),
        pyEval(r"'\n'", rep=REP_LAST, output='>>out'),
        pyEval(r"'%.2f ' % alleleFreq[1][0]", 
            outputExpr="'>>out%d' % grp"),
    ],
    end=2
)
print open('out').read()
print open('out1').read()
print open('out2').read()
#end
#file log/details_mating_demo.log
def lin_inc(gen, oldsize=[]):
    return [10+gen]*5

simu = simulator(
    population(subPop=lin_inc(1), loci=[1]),
    randomMating(newSubPopSizeFunc=lin_inc)
)
simu.evolve(
    ops = [
        stat(popSize=True),
        pyEval(r'"%d %d\n"%(gen, subPop[0]["popSize"])'),
    ],
    end=5
)
#end
#file log/details_mating_offspring.log
simu = simulator(
    population(size=10000, loci=[1]),
    randomMating(),
)
simu.evolve(
    preOps = [initByFreq([0.1, 0.9])],
    ops = [ ], end=100
)
simu.setMatingScheme(randomMating(numOffspring=2))
simu.addInfoFields(['father_idx', 'mother_idx'])
simu.setAncestralDepth(1)
simu.step(ops=[parentsTagger()])
pop = simu.getPopulation(0)
MaPenetrance(pop, locus=0, penetrance=[0.05, 0.1, 0.5])
AffectedSibpairSample(pop, size=100)
#end

