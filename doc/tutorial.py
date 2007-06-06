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
# #file log/tutorial_example2.log
# from simuPOP import *
# from simuRPy import *
# simu = simulator(
#     population(size=1000, ploidy=2, loci=[2]),
#     randomMating(),
#     rep = 3)
# simu.evolve(
#     preOps = [initByValue([1,2,2,1])],  
#     ops = [
#         recombinator(rate=0.1),
#         stat(LD=[0,1]),
#         varPlotter('LD[0][1]', numRep=3, step=10, 
#             saveAs='ld', ylim=[0,.25], 
#             lty=range(1, 4), col=range(2, 5),
#             xlab='generation', ylab='D', 
#             title='LD Decay'),
#     ],
#     end=100
# )
# #end
# #PS epstopdf ld10.eps
# #PS epstopdf ld20.eps
# #PS epstopdf ld30.eps
# #PS epstopdf ld40.eps
# #PS epstopdf ld50.eps
# #PS epstopdf ld60.eps
# #PS epstopdf ld70.eps
# #PS epstopdf ld80.eps
# #PS epstopdf ld90.eps
# #PS epstopdf ld100.eps
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
# save sample in Merlin QTDT format
from simuUtil import SaveQTDT
SaveQTDT(sample, output='sample', affectionCode=['U', 'A'], 
    fields=['affection'])
# have a look at the sample in Merlin-QTDT Format
print open('sample.map').read()
print open('sample.dat').read()
print open('sample.ped').read()
#end
#file log/tutorial_pop_variable.log
pop = population(subPop=[5, 10], loci=[5])
InitByFreq(pop, [.6, .3, .1])
Stat(pop, alleleFreq=[1], genoFreq=[2])
print pop.dvars().alleleFreq[1][0]
from simuUtil import ListVars
ListVars(pop.dvars(), useWxPython=False)
#end
#file log/tutorial_individual.log
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
#file log/tutorial_info.log
pop = population(100, loci=[5, 8],
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, randomMating(numOffspring=2))
simu.evolve(ops=[parentsTagger()], end=5)
ind = simu.population(0).individual(0)
ind1 = simu.population(0).individual(1)
print ind.info('father_idx'), ind.info('mother_idx')
print ind1.info('father_idx'), ind1.info('mother_idx')

#end
#file log/tutorial_op_stage.log
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
#file log/tutorial_op_gen.log
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
#file log/tutorial_op_rep.log
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
#file log/tutorial_op_output.log
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
#file log/tutorial_mating_demo.log
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
#file log/tutorial_mating_offspring.log
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
#file log/tutorial_hapmap.log
genes = [
    ("p53exon4", "rs1042522"),
    ("p53_6", "rs1625895"),
    ("xpdex23", "rs1799793"),
    ("xpdex10", "rs13181"),
    ("xpa", "rs1800975"),
    ("xpg1104", "rs17655"),
    ("xpf662", "rs2020955"),
    ("ercc61097", "rs2228526"),
    ("ercc61230", "rs4253211"),
    ("xpc_939", "rs2228001"),
    ("ccnh", "rs2266690"),
    ("rad23", "rs1805329"),
    ("ercc1", "rs3212986"),
    ("xpc_499", "rs2228000"), 
]

names = [x[1] for x in genes]
pops = []
for i in range(1, 23):
    print "Loading hapmap chromosome %d..." % i
    pop = LoadPopulation('hapmap_%d.bin' % i)
    markers = []
    for name in names:
        try:
            idx = pop.locusByName(name)
            markers.append(idx)
        except:
            pass
    if len(markers) > 0:
        markers.sort()
        pop.removeLoci(keep=markers)
        pops.append(pop)

all = MergePopulationsByLoci(pops)
#end
