#!/usr/bin/env python

#begin_file log/importSimuPOP.py
import simuOpt
simuOpt.setOptions(optimized=False, alleleType='long', quiet=False)
import simuPOP as sim
#end_file

#begin_file log/standard.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
#expect_error
pop = sim.population(10, loci=2)
pop.locusPos(10)
pop.individual(20).setAllele(1, 0)
#end_file

#begin_file log/simpleExample.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=1000, loci=2)
simu = sim.simulator(pop, sim.randomMating(ops=sim.recombinator(rates=0.01)), rep=3)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([1, 2, 2, 1])
    ],
    postOps = [
        sim.stat(LD=[0, 1]),
        sim.pyEval(r"'%.2f\t' % LD[0][1]", step=10),
        sim.pyOutput('\n', reps=-1, step=10)
    ],
    gen=100
)
#end_file

#begin_file log/help.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
help(sim.population.addInfoFields)
#end_file

#begin_file log/absIndex.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[10, 20], loci=[5, 7])
print pop.chromLocusPair(7)
print pop.absLocusIndex(1, 1)
print pop.absIndIndex(2, 1)
print pop.subPopIndPair(25)
#end_file

#begin_file log/iterator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=2, loci=[5, 6])
sim.InitByFreq(pop, [0.2, 0.3, 0.5])
for ind in pop.individuals():
    for loc in range(pop.chromBegin(1), pop.chromEnd(1)):
        print ind.allele(loc),
    print 
#end_file

#begin_file log/carray.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=2, loci=[3, 4])
sim.InitByFreq(pop, [0.3, 0.5, 0.2])
ind = pop.individual(0)
arr = ind.genotype()    # a carray to the underlying genotype
geno = list(arr)        # a list of alleles
print arr
print geno
arr.count(1)           # count
arr.index(2)           # index 
ind.setAllele(5, 3)    # change underlying genotype using setAllele
print arr              # arr is change
print geno             # but not geno
arr[2:5] = 4           # can use regular Python slice operation
print ind.genotype()
#end_file

#begin_file log/defdict.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population([100]*2, loci=1)
sim.InitByFreq(pop, [0, 0.2, 0.8], subPops=0)
sim.InitByFreq(pop, [0.2, 0.8], subPops=1)
sim.Stat(pop, alleleFreq=0, vars=['alleleFreq_sp'])
for sp in range(2):
    print 'Subpop %d (with %d alleles): ' % (sp, len(pop.dvars(sp).alleleFreq[0])),
    for a in range(3):
        print '%.2f ' % pop.dvars(sp).alleleFreq[0][a],
    print
#end_file

#begin_file log/genoStru.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2, 3], ploidy=2, loci=[5, 10],
    lociPos=range(0, 5) + range(0, 20, 2), chromNames=['Chr1', 'Chr2'],
    alleleNames=['A', 'C', 'T', 'G'])
# access genotypic information from the sim.population
pop.ploidy()
pop.ploidyName()
pop.numChrom()
pop.locusPos(2)
pop.alleleName(1)
# access from an individual
ind = pop.individual(2)
ind.numLoci(1)
ind.chromName(0)
ind.locusName(1)
# utility functions
ind.chromBegin(1)
ind.chromByName('Chr2')
# loci pos can be unordered within each chromosome
pop = sim.population(loci=[2, 3], lociPos=[3, 1, 1, 3, 2],
    lociNames=['loc%d' % x for x in range(5)])
pop.lociPos()
pop.lociNames()
#end_file

#begin_file log/haplodiploid.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2,5], ploidy=sim.Haplodiploid, loci=[3, 5])
sim.InitByFreq(pop, [0.3, 0.7])
sim.Dump(pop)
#end_file

#begin_file log/chromType.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=6, ploidy=2, loci=[3, 3, 6, 4, 4, 4],
    chromTypes=[sim.Autosome]*2 + [sim.ChromosomeX, sim.ChromosomeY] + [sim.Customized]*2)
sim.InitByFreq(pop, [0.3, 0.7])
sim.Dump(pop, structure=False) # does not display genotypic structure information
#end_file

#begin_file log/infoField.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(10, loci=[20], ancGen=1,
    infoFields=['father_idx', 'mother_idx'])
simu = sim.simulator(pop, sim.randomMating(ops=sim.recombinator(rates=0.01)))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([0]*20+[1]*20)
    ],
    duringOps = sim.parentsTagger(),
    gen = 1
)
pop = simu.extract(0)
pop.indInfo('mother_idx')  # mother of all offspring
ind = pop.individual(0)
mom = pop.ancestor(int(ind.mother_idx), 1)
print ind.genotype(0)
print mom.genotype(0)
print mom.genotype(1)
#end_file

#begin_file log/individual.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population([5, 4], loci=[2, 5], infoFields='x')
# get an individual
ind = pop.individual(3)
ind.ploidy()            # access to genotypic structure
ind.numChrom()
ind.affected()
ind.setAffected(True)   # access affection sim.status,
ind.sex()               # sex,
ind.setInfo(4, 'x')     # and information fields
ind.x = 5               # the same as ind.setInfo(4, 'x')
ind.info('x')           # get information field x
ind.x                   # the same as ind.info('x')
#end_file

#begin_file log/individual_genotype.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population([2, 1], loci=[2, 5])
for ind in pop.individuals(1):
    for marker in range(pop.totNumLoci()):
        ind.setAllele(marker % 2, marker, 0)
        ind.setAllele(marker % 2, marker, 1)
        print '%d %d ' % (ind.allele(marker, 0), ind.allele(marker, 1))

ind = pop.individual(1)
geno = ind.genotype(1)      # the second homologous copy
geno
geno[2] = 3
ind.genotype(1)
geno[2:4] = [3, 4]          # direct modification of the underlying genotype
ind.genotype(1)
# set genotype (genotype, ploidy, chrom)
ind.setGenotype([2, 1], 1, 1)
geno
#end_file

#begin_file log/subPop.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[3, 4, 5], ploidy=1, loci=1, infoFields='x')
# individual 0, 1, 2, ... will have an allele 0, 1, 2, ...
pop.setGenotype(range(pop.popSize()))
#
pop.subPopSize(1)
# merge subpopulations
pop.mergeSubPops([1, 2])
# split subpopulations
pop.splitSubPop(1, [2, 7])
pop.subPopSizes()
# set information field to each individual's new subpopulation ID
pop.setIndInfo([0, 1, 2, -1, 0, 1, 2, -1, -1, 0, 1, 2], 'x')
# this manually triggers an migration, individuals with negative values
# at this information field are removed.
pop.setSubPopByIndInfo('x')
sim.Dump(pop, width=2, structure=False)
#end_file

#begin_file log/subPopName.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[3, 4, 5], subPopNames=['x', 'y', 'z'])
pop.removeSubPops([1])
pop.subPopNames()
pop.subPopByName('z')
pop.splitSubPop(1, [2, 3])
pop.subPopNames()
pop.setSubPopName('z-1', 1)
pop.subPopNames()
pop.subPopByName('z')
#end_file

#begin_file log/virtualSplitter.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[200, 400], loci=[30], infoFields='x')
# assign random information fields
sim.InitSex(pop)
sim.InitInfo(pop, lambda: random.randint(0, 3), infoFields='x')
# define a virtual splitter by sex
pop.setVirtualSplitter(sim.sexSplitter())
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 1])    # Size of VSP 0 in subpopulation 0
# define a virtual splitter by information field 'x'
pop.setVirtualSplitter(sim.infoSplitter(field='x', values=[0, 1, 2, 3]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # Size of VSP 0 in subpopulation 0
pop.subPopSize([1, 0])    # Size of VSP 0 in subpopulation 1
#end_file

#begin_file log/virtualSubPop.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(10, loci=[2, 3], infoFields='Sex')
sim.InitSex(pop)
pop.setVirtualSplitter(sim.sexSplitter())
# initialize male and females with different genotypes. Set sim.initSex
# to False because this operator will by default also initialize sex.
sim.InitByValue(pop, [[0]*5, [1]*5], subPops=([0, 0], [0, 1]))
# set Sex information field to 0 for all males, and 1 for all females
pop.setIndInfo([sim.Male], 'Sex', [0, 0])
pop.setIndInfo([sim.Female], 'Sex', [0, 1])
# Print individual genotypes, followed by values at information field Sex
sim.Dump(pop, structure=False)
#end_file


#begin_file log/advancedVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[2000, 4000], loci=[30], infoFields='x')
# assign random information fields
sim.InitSex(pop)
sim.InitInfo(pop, lambda: random.randint(0, 3), infoFields='x')
#
# 1, use a combined splitter
pop.setVirtualSplitter(sim.combinedSplitter(splitters = [
    sim.sexSplitter(),
    sim.infoSplitter(field='x', values=[0, 1, 2, 3])
]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.Male
pop.subPopSize([1, 4])    # individuals in sp 1 with value 2 at field x
#
# use a product splitter that defines additional VSPs by sex and info
pop.setVirtualSplitter(sim.productSplitter(splitters = [
    sim.sexSplitter(names=['M', 'F']),  # give a new set of names
    sim.infoSplitter(field='x', values=[0, 1, 2, 3])
]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.Male with value 1 in sp 0
pop.subPopSize([1, 5])    # sim.Female with value 1 in sp 1
#
# use a combined splitter to join VSPs defined by a
# product splitter
pop.setVirtualSplitter(sim.combinedSplitter([
    sim.productSplitter([
        sim.sexSplitter(),
        sim.infoSplitter(field='x', values=[0, 1, 2, 3])])],
    vspMap = [[0,1,2], [4,5,6], [7]],
    names = ['sim.Male x<=3', 'sim.Female x<=3', 'sim.Female x=4']))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.Male with value 0, 1, 2 at field x
pop.subPopSize([1, 1])    # sim.Female with value 0, 1 or 2 at field x
#end_file


#begin_file log/accessIndividual.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
# create a sim.population with two generations. The current generation has values
# 0-9 at information field x, the parental generation has values 10-19.
pop = sim.population(size=[5, 5], loci=[2, 3], infoFields='x', ancGen=1)
pop.setIndInfo(range(11, 20), 'x')
pop1 = pop.clone()
pop1.setIndInfo(range(10), 'x')
pop.push(pop1)
#
ind = pop.individual(5)       # using absolute index
ind.x
ind.x       # the same as ind.x
# use a for loop, and relative index
for idx in range(pop.subPopSize(1)):
    print pop.individual(idx, 1).x,

# It is usually easier to use an iterator
for ind in pop.individuals(1):
    print ind.x,

# Access individuals in VSPs
pop.setVirtualSplitter(sim.infoSplitter(cutoff=[3, 7], field='x'))
for ind in pop.individuals([1, 1]):
    print ind.x,

# Access individuals in ancetral generations
pop.ancestor(5, 1).x        # absolute index
pop.ancestor(0, 1, 1).x     # relative index
# Or make ancestral generation the current generation and use 'individual'
pop.useAncestralGen(1)
pop.individual(5).x         # absolute index
pop.individual(0, 1).x      # relative index
# 'ancestor' can still access the 'present' (generation 0) generation
pop.ancestor(5, 0).x
# access individual by ID
pop.addInfoFields('ind_id')
sim.TagID(pop)
[int(ind.ind_id) for ind in pop.individuals()]
# access individual by ID. Note that individual 12 is in the parental generation
pop.indByID(12).x
#end_file

#begin_file log/batchAccess.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[4, 6], loci=2, infoFields='x')
pop.setIndInfo([random.randint(0, 10) for x in range(10)], 'x')
pop.indInfo('x')
pop.setGenotype([0, 1, 2, 3], 0)
pop.genotype(0)
pop.setVirtualSplitter(sim.infoSplitter(cutoff=[3], field='x'))
pop.setGenotype([0])    # clear all values
pop.setGenotype([5, 6, 7], [1, 1])
pop.indInfo('x', 1)
pop.genotype(1)
#end_file

#begin_file log/popInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(10)
pop.setInfoFields(['a', 'b'])
pop.addInfoFields('c')
pop.addInfoFields(['d', 'e'])
pop.infoFields()
#
# information fields can be accessed in batch mode
pop.setIndInfo([1], 'c')
# as well as individually.
for ind in pop.individuals():
    ind.e = ind.c + 1

print pop.indInfo('e')
#end_file

#begin_file log/ancestralPop.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(500, loci=1), sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    # start recording ancestral generations at generation 18
    preOps = sim.setAncestralDepth(2, at=[-2]),
    postOps = [
        sim.stat(alleleFreq=0, begin=-3),
        sim.pyEval(r"'%.3f\n' % alleleFreq[0][0]", begin=-3)
    ],
    gen = 20
)
pop = simu.population(0)
# start from current generation
for i in range(pop.ancestralGens(), -1, -1):
  pop.useAncestralGen(i)
  sim.Stat(pop, alleleFreq=0)
  print '%d   %.3f' % (i, pop.dvars().alleleFreq[0][0])

# restore to the current generation  
pop.useAncestralGen(0)  
#end_file

#begin_file log/addRemoveLoci.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(10, loci=3, chromNames=['chr1'])
# 1 1 1, 
pop.setGenotype([1])
# 1 1 1, 0 0 0
pop.addChrom(lociPos=[0.5, 1, 2], lociNames=['rs1', 'rs2', 'rs3'],
    chromName='chr2')
pop1 = sim.population(10, loci=3, chromNames=['chr3'],
    lociNames=['rs4', 'rs5', 'rs6'])
# 2 2 2,
pop1.setGenotype([2])
# 1 1 1, 0 0 0, 2 2 2
pop.addChromFrom(pop1)
# 1 1 1, 0 0 0, 2 0 2 2 0
pop.addLoci(chrom=[2, 2], pos=[1.5, 3.5], lociNames=['rs7', 'rs8'])
# 1 1 1, 0 0 0, 2 0 2 0
pop.removeLoci([8])
sim.Dump(pop)
#end_file

#begin_file log/recode.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(5, loci=[5], alleleNames=['A', 'T', 'C', 'G'])
sim.InitByFreq(pop, [0.2, 0.3, 0.4, 0.1])
sim.Dump(pop, structure=False)
print pop.genotype()
pop.recodeAlleles([0, 3, 1, 2], alleleNames=['A', 'C', 'G', 'T'])
sim.Dump(pop, structure=False)
print pop.genotype()
#end_file

#begin_file log/extract.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[10, 10], loci=[5, 5], infoFields=['x', 'y'])
sim.InitByValue(pop, range(10))
pop.setIndInfo([-1]*4 + [0]*3 + [-1]*3 + [2]*4 + [-1]*3 + [1]*4, 'x')
pop1 = pop.extract(field='x', loci=[1, 2, 3, 6, 7], infoFields='x')
sim.Dump(pop1, structure=False)
#end_file

#begin_file log/popVars.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
from pprint import pprint
pop = sim.population(100, loci=2)
sim.InitByFreq(pop, [0.3, 0.7])
print pop.vars()    # No variable now
pop.dvars().myVar = 21
print pop.vars()
sim.Stat(pop, popSize=1, alleleFreq=0)
# pprint prints in a less messy format
pprint(pop.vars())
# print number of allele 1 at locus 0
print pop.vars()['alleleNum'][0][1]
# use the dvars() function to access dictionary keys as attributes
print pop.dvars().alleleNum[0][1]
print pop.dvars().alleleFreq[0]
#end_file

#begin_file log/expression.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=1),
    sim.randomMating(), 5)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    postOps = [
        sim.stat(alleleFreq=0),
        sim.terminateIf('len(alleleFreq[0]) == 1')
    ]
)
#end_file

#begin_file log/savePop.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=5, chromNames=['chrom1'])
pop.dvars().name = 'my sim.population'
pop.save('sample.pop')
pop1 = sim.LoadPopulation('sample.pop')
pop1.chromName(0)
pop1.dvars().name
#begin_ignore
import os
os.remove('sample.pop')
#end_ignore
#end_file

#begin_file log/applicableGen.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(1000, loci=[20]), sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.8, 0.2])
    ],
    preOps = [
        sim.pyEval(r"'At the beginning of gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1])
    ],
    postOps = [
        sim.stat(alleleFreq=0, begin=80, step=10),
        sim.pyEval(r"'At the end of gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        sim.pyEval(r"'At the end of gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1])
    ],
    finalOps = sim.savePopulation(output='sample.pop'),
    gen=100
)
#begin_ignore
import os
os.remove('sample.pop')
#end_ignore
#end_file

#begin_file log/dryrun.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=[20]), sim.randomMating())
simu.evolve(
    initOps = sim.initByFreq([0.2, 0.8]),
    preOps = [
        sim.pyEval(r"'Around gen %d: alleleFreq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1])
    ],
    postOps = [
        sim.stat(alleleFreq=0, begin=80, step=10),
        sim.pyEval(r"'After gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        sim.pyEval(r"'Around gen %d: alleleFreq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1])
    ],
    finalOps = [sim.savePopulation(output='sample.pop')],
    gen=100,
    dryrun = True
)
#end_file

#begin_file log/replicate.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=[20]), sim.randomMating(), 5)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.2, 0.8])
    ],
    postOps = [
        sim.stat(alleleFreq=0, step=10),
        sim.pyEval('gen', step=10, reps=0),
        sim.pyEval(r"'\t%.2f' % alleleFreq[0][0]", step=10, reps=(0, 2, -1)),
        sim.pyOutput('\n', step=10, reps=-1)
    ],
    gen=30,
)
#end_file

#begin_file log/output.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=1000, loci=2),
    sim.randomMating(ops=sim.recombinator(rates=0.01)), rep=3)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([1, 2, 2, 1])
    ],
    postOps = [
        sim.stat(LD=[0, 1]),
        sim.pyEval(r"'%.2f\t' % LD[0][1]", step=20, output='>>LD.txt'),
        sim.pyOutput('\n', reps=-1, step=20, output='>>LD.txt'),
        sim.pyEval(r"'%.2f\t' % R2[0][1]", output='R2.txt'),
        sim.pyEval(r"'%.2f\t' % LD[0][1]", step=20, output="!'>>LD_%d.txt' % rep"),
    ],
    gen=100
)
print open('LD.txt').read()
print open('R2.txt').read()    # Only the last write operation succeed.
print open('LD_2.txt').read()  # Each replicate writes to a different file.
#begin_ignore
import os
for file in ['LD.txt', 'LD_0.txt', 'LD_1.txt', 'LD_2.txt', 'R2.txt']:
    os.remove(file)

#end_ignore
#end_file

#begin_file log/outputFunc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import logging
# logging to a file simulation.log, with detailed debug information
logging.basicConfig(
    filename='simulation.log',
    level=logging.DEBUG,
    format='%(levelname)s: %(message)s',
    filemode='w'
)
# logging to standard output with less information
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
console.setFormatter(formatter)
logger = logging.getLogger('')
logger.addHandler(console)
#
simu = sim.simulator(sim.population(size=1000, loci=2),
    sim.randomMating(ops=sim.recombinator(rates=0.01)))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([1, 2, 2, 1])
    ],
    postOps = [
        sim.stat(LD=[0, 1]),
        sim.pyEval(r"'LD: %d, %.2f' % (gen, LD[0][1])", step=20,
            output=logger.info),   # send LD to console and a logfile
        sim.pyEval(r"'R2: %d, %.2f' % (gen, R2[0][1])", step=20,
            output=logger.debug),  # send R2 only to a logfile
    ],
    gen=100
)
print open('simulation.log').read()
#begin_ignore
logging.shutdown()
import os
os.remove('simulation.log')
#end_ignore
#end_file

#begin_file log/transmitter.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=10000, loci=2),
    sim.randomMating(ops=[
        sim.mendelianGenoTransmitter(end=29),
        sim.recombinator(rates=0.01, begin=30),
    ])
)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([1, 2, 2, 1])
    ],
    postOps = [
        sim.stat(LD=[0, 1]),
        sim.pyEval(r"'gen %d, LD: %.2f\n' % (gen, LD[0][1])", step=20)
    ],
    gen=100
)
#end_file

#begin_file log/hybrid.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def myPenetrance(geno, gen):
    'A three-locus heterogeneity penetrance model'
    if sum(geno) < 2:
        return 0
    else:
        return sum(geno)*0.1

simu = sim.simulator(sim.population(1000, loci=[20]*3), sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.8, 0.2])
    ],
    postOps = [
        sim.pyPenetrance(func=myPenetrance, loci=[10, 30, 50]),
        sim.stat(numOfAffected=True),
        sim.pyEval(r"'%d: %d\n' % (gen, numOfAffected)")
    ],
    gen = 5
)
#end_file

#begin_file log/pyOperator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def dynaMutator(pop, param):
    '''This mutator mutates commom loci with low mutation rate and rare
    loci with high mutation rate, as an attempt to raise allele frequency
    of rare loci to an higher level.'''
    # unpack parameter
    (cutoff, mu1, mu2) = param;
    sim.Stat(pop, alleleFreq=range(pop.totNumLoci()))
    for i in range(pop.totNumLoci()):
        # Get the frequency of allele 1 (disease allele)
        if pop.dvars().alleleFreq[i][1] < cutoff:
            sim.KamMutate(pop, k=2, rates=mu1, loci=[i])
        else:
            sim.KamMutate(pop, k=2, rates=mu2, loci=[i])
    return True

simu = sim.simulator(sim.population(size=10000, loci=[2, 3]),
    sim.randomMating())
simu.evolve(
    initOps = [ 
        sim.initSex(),
        sim.initByFreq([.99, .01], loci=[0, 2, 4]),
        sim.initByFreq([.8, .2], loci=[1, 3])
    ],
    preOps = sim.pyOperator(func=dynaMutator, param=(.2, 1e-2, 1e-5)),
    postOps = [
        sim.stat(alleleFreq=range(5), step=10),
        sim.pyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)                
#end_file

#begin_file log/pyDuringMatingOperator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def rejectInd(off):
    'reject an individual if it off.allele(0) == 1'
    return off.allele(0) == 0

simu = sim.simulator(sim.population(size=100, loci=1),
    sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    duringOps = sim.pyOperator(func=rejectInd, offspringOnly=True),
    gen = 1
)
# You should see no individual with allele 1 at locus 0, ploidy 0.
simu.population(0).genotype()[:20]
#end_file

#begin_file log/newOperator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
class dynaMutator(sim.pyOperator):
    '''This mutator mutates commom loci with low mutation rate and rare
    loci with high mutation rate, as an attempt to raise allele frequency
    of rare loci to an higher level.'''
    def __init__(self, cutoff, mu1, mu2, *args, **kwargs):
        self.cutoff = cutoff
        self.mu1 = mu1
        self.mu2 = mu2
        sim.pyOperator.__init__(self, func=self.mutate, *args, **kwargs)
    #
    def mutate(self, pop):
        sim.Stat(pop, alleleFreq=range(pop.totNumLoci()))
        for i in range(pop.totNumLoci()):
            # Get the frequency of allele 1 (disease allele)
            if pop.dvars().alleleFreq[i][1] < self.cutoff:
                sim.KamMutate(pop, k=2, rates=self.mu1, loci=[i])
            else:
                sim.KamMutate(pop, k=2, rates=self.mu2, loci=[i])
        return True

simu = sim.simulator(sim.population(size=10000, loci=[2, 3]),
    sim.randomMating())
simu.evolve(
    initOps = [ 
        sim.initSex(),
        sim.initByFreq([.99, .01], loci=[0, 2, 4]),
        sim.initByFreq([.8, .2], loci=[1, 3])
    ],
    preOps = dynaMutator(cutoff=.2, mu1=1e-2, mu2=1e-5),
    postOps = [
        sim.stat(alleleFreq=range(5), step=10),
        sim.pyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)          
#end_file

#begin_file log/funcInitByFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def InitByFreq(pop, *args, **kwargs):
    sim.initByFreq(*args, **kwargs).apply(pop)

pop = sim.population(1000, loci=[2,3])
sim.InitByFreq(pop, [.2, .3, .5])
#end_file

#begin_file log/migrSize.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=[500, 1000], infoFields='migrate_to'),
    sim.randomMating())

simu.evolve(
    initOps = sim.initSex(),
    preOps = sim.migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
    postOps = [
        sim.stat(popSize=True),
        sim.pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end_file

#begin_file log/migrFixedSize.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=[500, 1000], infoFields='migrate_to'),
    sim.randomMating(subPopSize=[500, 1000]))
simu.evolve(
    initOps = sim.initSex(),
    preOps = [
        sim.migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        sim.stat(popSize=True),
        sim.pyEval(r'"%s\n" % subPopSize')
    ],
    postOps = [
        sim.stat(popSize=True),
        sim.pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end_file

#begin_file log/demoFunc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def demo(gen, pop=None):
    return [500 + gen*10, 1000 + gen*10]

simu = sim.simulator(
    sim.population(size=[500, 1000], infoFields='migrate_to'),
    sim.randomMating(subPopSize=demo))
simu.evolve(
    initOps = sim.initSex(),
    preOps = sim.migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
    postOps = [
        sim.stat(popSize=True),
        sim.pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end_file


#begin_file log/advancedDemoFunc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def demo(gen, pop):
    if gen < 2:
        return 1000 + 100 * gen
    if gen == 2:
        # this happens right before mating at generation 2
        size = pop.popSize()
        pop.splitSubPop(0, [size/2, size - size/2]) 
    # for generation two and later
    return [x + 50 * gen for x in pop.subPopSizes()]

simu = sim.simulator(sim.population(1000),
    sim.randomSelection(subPopSize=demo))
simu.evolve(
    preOps = [
        sim.stat(popSize=True),
        sim.pyEval(r'"Gen %d:\t%s (before mating)\t" % (gen, subPopSize)')
    ],
    postOps = [
        sim.stat(popSize=True),
        sim.pyEval(r'"%s (after mating)\n" % subPopSize')
    ],
    gen = 5
)
#end_file

#begin_file log/numOff.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def checkNumOffspring(ms):
    '''Check the number of offspring for each family using
       information field father_idx
    '''
    simu = sim.simulator(
        sim.population(size=[30], infoFields=['father_idx', 'mother_idx']),
        matingScheme=ms)
    simu.evolve(
        initOps = sim.initSex(),
        duringOps = sim.parentsTagger(),
        gen=1)
    # get the parents of each offspring
    parents = [(x, y) for x, y in zip(simu.population(0).indInfo('mother_idx'),
        simu.population(0).indInfo('father_idx'))]
    # Individuals with identical parents are considered as siblings.
    famSize = []
    lastParent = (-1, -1)
    for parent in parents:
        if parent == lastParent:
            famSize[-1] += 1
        else:
            lastParent = parent
            famSize.append(1)
    return famSize

# Case 1: produce the given number of offspring
checkNumOffspring(sim.randomMating(numOffspring=2))
# Case 2: Use a Python function
import random
def func(gen):
    return random.randint(5, 8)

checkNumOffspring(sim.randomMating(numOffspring=func))
# Case 3: A geometric distribution
checkNumOffspring(sim.randomMating(numOffspring=(sim.GeometricDistribution, 0.3)))
# Case 4: A Possition distribution
checkNumOffspring(sim.randomMating(numOffspring=(sim.PoissonDistribution, 3)))
# Case 5: A Binomial distribution
checkNumOffspring(sim.randomMating(numOffspring=(sim.BinomialDistribution, 0.1, 10)))
# Case 6: A uniform distribution
checkNumOffspring(sim.randomMating(numOffspring=(sim.UniformDistribution, 2, 6)))
#end_file

#begin_file log/sexMode.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def checkSexMode(ms):
    '''Check the assignment of sex to offspring'''
    simu = sim.simulator(
        sim.population(size=[40]),
        matingScheme=ms)
    simu.evolve(initOps = sim.initSex(), gen=1)
    # return individual sex as a string
    return ''.join([ind.sexChar() for ind in simu.population(0).individuals()])

# Case 1: sim.NoSex (all male, sim.randomMating will not continue)
checkSexMode(sim.randomMating(sexMode=sim.NoSex))
# Case 2: sim.RandomSex (sim.Male/Female with probability 0.5)
checkSexMode(sim.randomMating(sexMode=sim.RandomSex))
# Case 3: sim.ProbOfMale (Specify probability of male)
checkSexMode(sim.randomMating(sexMode=(sim.ProbOfMale, 0.8)))
# Case 4: sim.NumOfMale (Specify number of male in each family)
checkSexMode(sim.randomMating(numOffspring=3, sexMode=(sim.NumOfMale, 1)))
# Case 5: sim.NumOfFemale (Specify number of female in each family)
checkSexMode(sim.randomMating(
    numOffspring=(sim.UniformDistribution, 4, 6),
    sexMode=(sim.NumOfFemale, 2))
)
#end_file

#begin_file log/monogamous.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(20, infoFields=['father_idx', 'mother_idx']),
    sim.monogamousMating(numOffspring=2, sexMode=(sim.NumOfMale, 1)))
simu.evolve(
    initOps = sim.initSex(sex=(sim.Male, sim.Female)),
    duringOps = sim.parentsTagger(),
    gen = 5
)
pop = simu.extract(0)
[ind.sex() for ind in pop.individuals()]
[int(ind.father_idx) for ind in pop.individuals()]
[int(ind.mother_idx) for ind in pop.individuals()]
# count the number of distinct parents
len(set(pop.indInfo('father_idx')))
len(set(pop.indInfo('mother_idx')))
#end_file

#begin_file log/polygamous.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, infoFields=['father_idx', 'mother_idx']),
    sim.polygamousMating(polySex=sim.Male, polyNum=2))
simu.evolve(
    initOps = sim.initSex(),
    duringOps = sim.parentsTagger(),
    gen = 5
)
pop = simu.extract(0)
[int(ind.father_idx) for ind in pop.individuals()][:20]
[int(ind.mother_idx) for ind in pop.individuals()][:20]
#end_file

#begin_file log/randomSelection.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, ploidy=1, loci=[5, 5], ancGen=1,
    infoFields='parent_idx'),
    sim.randomSelection())
simu.evolve(
    initOps = sim.initByFreq([0.3, 0.7]),
    duringOps = sim.parentsTagger(infoFields='parent_idx'),
    gen = 5
)
pop = simu.extract(0)
ind = pop.individual(0)
par = pop.ancestor(int(ind.parent_idx), 1)
print ind.sex(), ind.genotype()
print par.sex(), par.genotype()
#end_file

#begin_file log/alphaMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(1000, loci=5, 
    infoFields=['father_idx', 'mother_idx', 'fitness']),
    sim.alphaMating(alphaSex=sim.Male, alphaNum=2))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    preOps = sim.maSelector(loci=0, fitness=[0.8, 0.8, 1]),
    duringOps = sim.parentsTagger(),
    postOps = [
        sim.stat(alleleFreq=0),
        sim.pyEval(r'"%.2f\n" % alleleFreq[0][1]', step=5)
    ],
    gen = 20,
)
pop = simu.extract(0)
[int(ind.father_idx) for ind in pop.individuals()][:10]
[int(ind.mother_idx) for ind in pop.individuals()][:10]
#end_file

#begin_file log/haplodiploidMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(10, ploidy=sim.Haplodiploid, loci=[5, 5],
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(sim.sexSplitter())
simu = sim.simulator(pop, sim.haplodiploidMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([0]*10, subPops=[(0, 0)]),
        sim.initByValue([1]*10+[2]*10, subPops=[(0, 1)])
    ],
    preOps = sim.dumper(structure=False),
    duringOps = sim.parentsTagger(),
    postOps = sim.dumper(structure=False),
    gen = 1
)
#end_file

#begin_file log/selfMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(20, loci=8)
# every chromosomes are different. :-)
for idx, ind in enumerate(pop.individuals()):
    ind.setGenotype([idx*2], 0)
    ind.setGenotype([idx*2+1], 1)

simu = sim.simulator(pop, sim.selfMating(ops=sim.recombinator(rates=0.01)))
simu.evolve( gen = 1)
sim.Dump(simu.population(0), width=3, structure=False, max=10)
#end_file

#begin_file log/heteroMatingSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000, 1000], loci=2,
    infoFields=['father_idx', 'mother_idx'])
simu = sim.simulator(pop, sim.heteroMating(
    [sim.randomMating(numOffspring=2, subPops=0),
     sim.randomMating(numOffspring=4, subPops=1)
    ])
)
simu.evolve(
    initOps = sim.initSex(),
    duringOps = sim.parentsTagger(),
    gen=10
)
pop = simu.extract(0)
[int(ind.father_idx) for ind in pop.individuals(0)][:10]
[int(ind.father_idx) for ind in pop.individuals(1)][:10]
#end_file

#begin_file log/heteroMatingVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000], loci=2,
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(sim.proportionSplitter([0.2, 0.8]))
simu = sim.simulator(pop, sim.heteroMating(
    matingSchemes = [
        sim.selfMating(subPops=[(0, 0)]),
        sim.randomMating(subPops=[(0, 1)])
    ])
)
simu.evolve(
    initOps = sim.initSex(),
    duringOps = sim.parentsTagger(),
    gen = 10
)
pop = simu.extract(0)
[int(ind.father_idx) for ind in pop.individuals(0)][:15]
[int(ind.mother_idx) for ind in pop.individuals(0)][:15]
#end_file

#begin_file log/heteroMatingWeight.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000], loci=2,
    infoFields='mark')
pop.setVirtualSplitter(sim.rangeSplitter([[0, 500], [200, 1000]]))
def markOff(param):
    '''define a Python during mating operator that marks
       individual information field 'mark'
    '''
    def func(off, param):
        off.mark = param
        return True
    return sim.pyOperator(func=func, param=param, offspringOnly=True)

simu = sim.simulator(pop, sim.heteroMating(
    matingSchemes = [
        sim.randomMating(subPops=0, weight=-0.5, ops=[markOff(0), sim.mendelianGenoTransmitter()]),
        sim.randomMating(subPops=[(0, 0)], weight=2, ops=[markOff(1), sim.mendelianGenoTransmitter()]),
        sim.randomMating(subPops=[(0, 1)], weight=3, ops=[markOff(2), sim.mendelianGenoTransmitter()])
    ])
)
simu.evolve(
    initOps = sim.initSex(),
    gen = 10
)
marks = list(simu.extract(0).indInfo('mark'))
marks.count(0.)
marks.count(1.)
marks.count(2.)
#end_file

#begin_file log/randomMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def randomMating(numOffspring = 1., sexMode = sim.RandomSex,
        preOps = sim.mendelianGenoTransmitter(), subPopSize = [],
        subPops = AllAvail, weight = 0, selectionField = 'fitness'):
    'A basic diploid sexual random mating scheme.'
    return sim.homoMating(
        chooser = sim.randomParentsChooser(True, selectionField),
        generator = sim.offspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPops = subPops,
        weight = weight)
#end_file

#begin_file log/sequentialSelfing.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=5*3, infoFields='parent_idx'),
    sim.homoMating(sim.sequentialParentChooser(),
    sim.offspringGenerator(ops=sim.selfingGenoTransmitter())))
simu.evolve(
    initOps = [sim.initByFreq([0.2]*5)],
    preOps = sim.dumper(structure=False, max=5),
    duringOps = sim.parentsTagger(infoFields='parent_idx'),
    postOps = sim.dumper(structure=False, max=5),
    gen = 1
)
#end_file

#begin_file log/controlledOffGenerator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def traj(gen):
    return [0.5 + gen * 0.01]

simu = sim.simulator(sim.population(1000, loci=[10]*2),
    sim.homoMating(sim.randomParentChooser(),
        sim.controlledOffspringGenerator(loci=5,
            alleles=[0], freqFunc=traj,
            ops = sim.selfingGenoTransmitter()))
)
# evolve the sim.population while keeping allele frequency 0.5
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    postOps = [
        sim.stat(alleleFreq=[5, 15]),
        sim.pyEval(r'"%.2f\t%.2f\n" % (alleleFreq[5][0], alleleFreq[15][0])')
    ],
    gen = 5
)
#end_file

#begin_file log/mitochondrial.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(10, loci=[5]*5,
    # one autosome, two sex chromosomes, and two mitochondrial chromosomes
    chromTypes=[sim.Autosome, sim.ChromosomeX, sim.ChromosomeY] + [sim.Customized]*2,
    infoFields=['father_idx', 'mother_idx'])
simu = sim.simulator(pop, sim.randomMating(ops= [
    sim.recombinator(rates=0.1),
    sim.mitochondrialGenoTransmitter()]))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.4] + [0.2]*3)
    ],
    duringOps = sim.parentsTagger(),
    postOps = sim.dumper(structure=False),
    gen = 2
)
#end_file

#begin_file log/sexSpecificRec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
class sexSpecificRecombinator(sim.pyOperator):
    def __init__(self, intensity=0, rates=0, loci=[], convMode=sim.NoConversion,
            maleIntensity=0, maleRates=0, maleLoci=[], maleConvMode=sim.NoConversion,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.recombinator = sim.recombinator(rates, intensity, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = sim.recombinator(maleRates, maleIntensity,
            maleLoci, maleConvMode)
        #
        self.initialized = False
        #
        sim.pyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)
    #
    def transmitGenotype(self, pop, off, dad, mom):
        # Recombinators need to be initialized. Basically, they cache some
        # sim.population properties to speed up genotype transmission.
        if not self.initialized:
            self.recombinator.initialize(pop)
            self.maleRecombinator.initialize(pop)
            self.initialized = True
        # Form the first homologous copy of offspring.
        self.recombinator.transmitGenotype(mom, off, 0)
        self.maleRecombinator.transmitGenotype(dad, off, 1)
        return True

pop = sim.population(10, loci=[15]*2, infoFields=['father_idx', 'mother_idx'])
simu = sim.simulator(pop, sim.randomMating(
    ops = sexSpecificRecombinator(rates=0.1, maleRates=0)))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.4] + [0.2]*3)
    ],
    duringOps = sim.parentsTagger(),
    postOps = sim.dumper(structure=False),
    gen = 2
)
#end_file

#begin_file log/infoChooser.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=[10],
    infoFields=['father_idx', 'mother_idx', 'sibling'])
pop.setVirtualSplitter(sim.sexSplitter())
def locate_sibling(pop):
    '''The sim.population is arranged as MFMFMFMF... where MF are siblings, so the
    sibling of males are 1, 3, 5, .. and the slibling of females are 0, 2, 4, ...
    '''
    pop.setIndInfo([2*x+1 for x in range(pop.popSize()/2)], 'sibling', (0, 0))
    pop.setIndInfo([2*x for x in range(pop.popSize()/2)], 'sibling', (0, 1))

simu = sim.simulator(pop, sim.consanguineousMating(func=locate_sibling, infoFields='sibling',
    numOffspring=2, sexMode=(sim.NumOfMale, 1)))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.2, 0.8])
    ],
    duringOps = sim.parentsTagger(),
    postOps = sim.dumper(structure=False, max=6, at=[-1]),
    gen = 2
)
#end_file

#begin_file log/generator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def func():
    i = 1
    all = 0
    while i <= 5:
        all += 1./i
        i += 1
        yield all 

for i in func():
    print '%.3f' % i,
#end_file

#begin_file log/pyParentsChooser.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
from random import randint
def randomChooser(pop, sp):
    males = []
    females = []
    # identify males and females in each social rank
    for rank in range(3):
        males.append([x for x in pop.individuals(sp) \
            if x.sex() == sim.Male and x.rank == rank])
        females.append([x for x in pop.individuals(sp) \
            if x.sex() == sim.Female and x.rank == rank])
    #
    while True:
        # choose a rank randomly
        rank = int(pop.individual(randint(0, pop.subPopSize(sp) - 1), sp).rank)
        yield males[rank][randint(0, len(males[rank]) - 1)], \
            females[rank][randint(0, len(females[rank]) - 1)]

def setRank(pop, dad, mom, off):
    'The rank of offspring can increase or drop to zero randomly'
    off.rank = (dad.rank + randint(-1, 1)) % 3

pop = sim.population(size=[1000, 2000], loci=1, infoFields='rank')
simu = sim.simulator(pop, sim.homoMating(
    sim.pyParentsChooser(randomChooser),
    sim.offspringGenerator(ops=sim.mendelianGenoTransmitter())))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initInfo(lambda : randint(0, 2), infoFields='rank')
    ],
    gen = 5
)    
#end_file


#begin_file log/cppParentChooser.py
#expect_error because some system does not have swig or vc to compile
#begin_ignore
classFile = open('log/myParentsChooser.h', 'w')
classFile.write('''#include <stdlib.h>
#include <vector>
#include <utility>
using std::pair;
using std::vector;
class myParentsChooser
{
public:
        // A constructor takes all locations of male and female.
        myParentsChooser(const std::vector<int> & m, const std::vector<int> & f)
                : male_idx(m), female_idx(f)
        {
                srand(time(0));
        }

        pair<unsigned long, unsigned long> chooseParents()
        {
                unsigned long male = rand() % male_idx.size();
                unsigned long female = rand() % male_idx.size();
                return std::make_pair(male, female);
        }
private:
        vector<int> male_idx;
        vector<int> female_idx;
};
''')
classFile.close()
interFile = open('log/myParentsChooser.i', 'w')
interFile.write('''%module myParentsChooser
%{
#include "myParentsChooser.h"
%}
// std_vector.i for std::vector
%include "std_vector.i"
%template() std::vector<int>;
// stl.i for std::pair
%include "stl.i"
%template() std::pair<unsigned long, unsigned long>;
%include "myParentsChooser.h"
''')
interFile.close()
setupFile = open('log/setup.py', 'w')
setupFile.write('''from distutils.core import setup, Extension
import sys
# Under linux/gcc, lib stdc++ is needed for C++ based extension.
if sys.platform == 'linux2':
    libs = ['stdc++']
else:
    libs = []
setup(name = "myParentsChooser",
    description = "A sample parent chooser",
    py_modules = ['myParentsChooser'],  # will be generated by SWIG
    ext_modules = [
        Extension('_myParentsChooser',
            sources = ['myParentsChooser.i'],
            swig_opts = ['-O', '-templatereduce', '-shadow',
                '-python', '-c++', '-keyword', '-nodefaultctor'],
            include_dirs = ["."],
    )
  ]
)
''')
setupFile.close()
import os, sys
sys.path.append('log')
try:
    import myParentsChooser
except:
    os.chdir('log')
    os.system('python setup.py build_ext --swig-opts="-O -templatereduce -shadow -c++ -keyword -nodefaultctor" install --install-purelib="." --install-platlib="."')
    os.chdir('..')
#end_ignore

#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore

# The class myParentsChooser is defined in module myParentsChooser
try:
    from myParentsChooser import myParentsChooser
except ImportError:
    # if failed to import the C++ version, use a Python version
    import random
    class myParentsChooser:
        def __init__(self, maleIndexes, femaleIndexes):
            self.maleIndexes = maleIndexes
            self.femaleIndexes = femaleIndexes
        def chooseParents(self):
            return self.maleIndexes[random.randint(0, len(self.maleIndexes)-1)],\
                self.femaleIndexes[random.randint(0, len(self.femaleIndexes)-1)]

def parentsChooser(pop, sp):
    'How to call a C++ level parents chooser.'
    # create an object with needed information (such as x, y) ...
    pc = myParentsChooser(
        [x for x in range(pop.popSize()) if pop.individual(x).sex() == sim.Male],
        [x for x in range(pop.popSize()) if pop.individual(x).sex() == sim.Female])
    while True:
        # return indexes of parents repeatedly
        yield pc.chooseParents()

pop = sim.population(100, loci=1)
simu = sim.simulator(pop,
    sim.homoMating(sim.pyParentsChooser(parentsChooser),
    sim.offspringGenerator(ops=sim.mendelianGenoTransmitter()))
)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    gen = 100
)
#end_file

#begin_file log/simuGen.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(50, loci=[10], ploidy=1),
    sim.randomSelection(), rep=3)
simu.evolve(gen = 5)
simu.gen()
simu.evolve(
    initOps = [sim.initByFreq([0.5, 0.5])],
    postOps = [
        sim.stat(alleleFreq=5),
        sim.ifElse('alleleNum[5][0] == 0',
            sim.pyEval(r"'Allele 0 is lost in rep %d at gen %d\n' % (rep, gen)")),
        sim.ifElse('alleleNum[5][0] == 50',
            sim.pyEval(r"'Allele 0 is fixed in rep %d at gen %d\n' % (rep, gen)")),
        sim.terminateIf('len(alleleNum[5]) == 1'),
    ],
)
simu.gen()
#end_file

#begin_file log/describe.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
pop = sim.population(10000, loci=100, infoFields=['ind_id', 'age'])
pop.setVirtualSplitter(sim.infoSplitter(field='age', cutoff=[20, 50, 75]))

def outputStat(pop):
    'Calculate and output statistics, ignored'
    return True

simu = sim.simulator(pop,
    sim.heteroMating([
        sim.cloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
        sim.randomMating(ops=[
            sim.idTagger(),
            sim.recombinator(intensity=1e-4)
        ], subPops=[(0,1)]),
    ])
)
# describe this evolutionary process
print simu.describe(
    initOps = [
        sim.initSex(),
        sim.initInfo(lambda: random.randint(0, 75), infoFields='age'),
        sim.initByFreq([0.5, 0.5]),
        sim.idTagger(),
        sim.pyOutput('Prevalence of disease in each age group:\n'),
    ],
    preOps = sim.infoExec('age += 1'),
    duringOps = sim.pedigreeTagger(),
    postOps = [
        sim.maPenetrance(loci=0, penetrance=[0.01, 0.1, 0.3]),
        sim.pyOperator(func=outputStat)
    ],
    gen = 100
)     
#end_file

#begin_file log/twoStage.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
# First stage: use the standard random mating scheme, do not use any
# information field for efficiency considerations.
simu = sim.simulator(sim.population(500, loci=[10]), sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    gen = 50
)
# Second stage: track parents and produce more offspring per mating
# event. In preparation for sim.pedigree ascertainment.
for pop in simu.populations():
    pop.addInfoFields(['father_idx', 'mother_idx'])
    pop.setAncestralDepth(1)

simu.setMatingScheme(sim.randomMating(numOffspring=2))
simu.evolve(
    duringOps = sim.parentsTagger(),
    postOps = sim.maPenetrance(loci=0, penetrance=(0.2, 0.4, 0.5)),
    gen = 5
)
# Sample affected sibpairs
pop = simu.extract(0)
from simuPOP.sampling import AffectedSibpairSample
sample = AffectedSibpairSample(pop, size=5)[0]
[int(ind.father_idx) for ind in sample.individuals()]
#end_file

#begin_file log/changeStru.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
def mutator(pop, param):
    'Parameter has a length of region and a mutation rate at each basepair'
    region, rate = param
    # there are certainly more efficient algorithm, but this 
    # example just mutate each basepair one by one....
    for i in range(region):
        if random.random() < rate:
            try:
                idx = pop.addLoci(chrom=0, pos=i)[0]
            except:
                # position might duplicate
                continue
            # choose someone to mutate
            ind = pop.individual(random.randint(0, pop.popSize() - 1))
            ind.setAllele(1, idx)
    return True

# The sim.populations start with no loci at all.
simu = sim.simulator(sim.population(1000, loci=[]), sim.randomMating(), rep=3)
simu.evolve(
    initOps = sim.initSex(),
    postOps = sim.pyOperator(func=mutator, param=(10000, 2e-6)),
    gen = 200
)
for pop in simu.populations():
    print pop.totNumLoci(), pop.lociPos()
#end_file

#begin_file log/simuFunc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=[5, 10], infoFields='x'),
    sim.randomMating(), rep=5)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.4, 0.6])
    ],
    gen=10
)
# clone
cloned = simu.clone()
# save and load, using a different mating scheme
simu.save("sample.sim")
loaded = sim.LoadSimulator("sample.sim", sim.randomMating(numOffspring=2))
# 
simu.numRep()
loaded.numRep()
for pop1,pop2 in zip(cloned.populations(), loaded.populations()):
    assert pop1 == pop2

# continue to evolve
simu.evolve(gen=10)
simu.gen()
#begin_ignore
import os
os.remove('sample.sim')
#end_ignore
#end_file

#begin_file log/initSex.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000, 1000])
sim.InitSex(pop, maleFreq=0.3, subPops=0)
sim.InitSex(pop, sex=[sim.Male, sim.Female, sim.Female], subPops=1)
sim.Stat(pop, numOfMale=True, vars='numOfMale_sp')
print pop.dvars(0).numOfMale
print pop.dvars(1).numOfMale
#end_file

#begin_file log/initByFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2, 3], loci=[5, 7])
sim.InitByFreq(pop, alleleFreq=[[.2, .8], [.8, .2]])
sim.Dump(pop, structure=False)
#end_file

#begin_file log/initByFreqIdenticalInds.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2, 3], loci=[5, 7])
sim.InitByFreq(pop, alleleFreq=[.2, .8], identicalInds=True)
sim.Dump(pop, structure=False)
#end_file

#begin_file log/initByValue.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2, 3], loci=[5, 7])
sim.InitByValue(pop, [1]*5 + [2]*7 + [3]*5 +[4]*7)
sim.Dump(pop, structure=False)
#end_file

#begin_file log/initByValueProp.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[6, 8], loci=[5, 7])
pop.setVirtualSplitter(sim.sexSplitter())
# initialize sex and the first two loci
sim.InitSex(pop)
sim.InitByValue(pop, loci=range(5), value=range(10))
# initialize all males
sim.InitByValue(pop, loci=range(5, 12), value=[2]*7,
    subPops=[(0, 0), (1, 0)])
# initialize females by proportion
sim.InitByValue(pop, loci=range(5, 12), ploidy=1, value=[[3]*7, [4]*7],
    subPops=[(0, 1), (1, 1)], proportions=[0.4, 0.6])
sim.Dump(pop, structure=False)
#end_file

#begin_file log/initInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[5], loci=[2], infoFields=['sex', 'age'])
pop.setVirtualSplitter(sim.sexSplitter())
sim.InitSex(pop)
sim.InitInfo(pop, 0, subPops=[(0,0)], infoFields='sex')
sim.InitInfo(pop, 1, subPops=[(0,1)], infoFields='sex')
sim.InitInfo(pop, lambda: random.randint(20, 70), infoFields='age')
sim.Dump(pop, structure=False)
#end_file

#begin_file log/dumper.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[10, 10], loci=[20, 30], infoFields='gen',
    ancGen=-1)
pop.setVirtualSplitter(sim.sexSplitter())
pop1 = pop.clone()
sim.InitByFreq(pop, [0]*20 + [0.1]*10)
pop.setIndInfo(1, 'gen')
sim.InitByFreq(pop1, [0]*50 + [0.1]*10)
pop1.setIndInfo(2, 'gen')
pop.push(pop1)
sim.Dump(pop, width=3, loci=[5, 6, 30], subPops=([0, 0], [1, 1]),
    max=10, structure=False, ancGen=-1)
#end_file

#begin_file log/savePopulation.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=2),
    sim.randomMating(), rep=5)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.2, 0.8])
    ],
    postOps = sim.savePopulation(output="!'snapshot_%d_%d.pop' % (rep, gen)",
            step = 10),
    gen = 50
)
#begin_ignore
import os
for rep in range(5):
    for gen in range(0, 50, 10):
        os.remove('snapshot_%d_%d.pop' % (rep, gen))

#end_ignore
#end_file

#begin_file log/setAncDepth.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, infoFields=['father_idx', 'mother_idx']),
    sim.randomMating(), rep=5)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.3, 0.7])
    ],
    preOps = sim.setAncestralDepth(2, at=-2),
    duringOps = sim.parentsTagger(begin=-2),
    gen = 100
)
pop = simu.population(3)
print pop.ancestralGens()
print pop.ancestor(10, 1).father_idx
#end_file


#begin_file log/ifElseFixed.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=1000, loci=1),
    sim.randomMating())
verbose = True
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5]),
    ],
    postOps = sim.ifElse(verbose,
        ifOps = [
            sim.stat(alleleFreq=0),
            sim.pyEval(r"'Gen: %3d, allele freq: %.3f\n' % (gen, alleleFreq[0][1])",
                step=5)
        ],
        begin=10),
    gen = 30
)
#end_file


#begin_file log/ifElse.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=1000, loci=1),
    sim.randomMating(), rep=4)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5]),
        sim.pyExec('below40, above60 = 0, 0')
    ],
    postOps = [
        sim.stat(alleleFreq=0),
        sim.ifElse('alleleFreq[0][1] < 0.4',
            sim.pyExec('below40 += 1')),
        sim.ifElse('alleleFreq[0][1] > 0.6',
            sim.pyExec('above60 += 1')),
        sim.ifElse('len(alleleFreq[0]) == 1',
            sim.pyExec('stoppedAt = gen')),
        sim.terminateIf('len(alleleFreq[0]) == 1')
    ]
)
for pop in simu.populations():
    print 'Overall: %4d, below 40%%: %4d, above 60%%: %4d' % \
        (pop.dvars().stoppedAt, pop.dvars().below40, pop.dvars().above60)
#end_file

#begin_file log/terminateIf.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=100, loci=1),
    sim.randomMating(), rep=10)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5]),
    ],
    postOps = [
        sim.stat(alleleFreq=0),
        sim.terminateIf('len(alleleFreq[0]) == 1', stopAll=True)
    ]
)

#end_file

#begin_file log/debug.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=1), sim.randomMating(), rep=5)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.1, 0.9])
    ],
    postOps = [
        sim.stat(alleleFreq=0),
        sim.ifElse('alleleNum[0][0] == 0',
            ifOps = [
                sim.turnOnDebug("DBG_MUTATOR"),
                sim.pointMutator(loci=0, allele=0, inds=0),
            ],
            elseOps = sim.turnOffDebug("DBG_MUTATOR")),
    ],
    gen = 100
)
#begin_ignore
sim.TurnOffDebug("DBG_MUTATOR")
#end_ignore
#end_file

#begin_file log/pause.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100), sim.randomMating(), rep=10)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    postOps = [sim.pause(stopOnKeyStroke=str(x), reps=x) for x in range(10)],
    gen = 100
)
#end_file

#begin_file log/ticToc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(10000, loci=[100]*5), sim.randomMating(), rep=2)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.1, 0.9])
    ],
    postOps = [
        sim.stat(alleleFreq=0),
        sim.ticToc(step=50, reps=-1),
    ],
    gen = 101
)
#end_file

#begin_file log/pyExec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=1),
    sim.randomMating(), rep=2)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.2, 0.8]),
        sim.pyExec('traj=[]')
    ],
    postOps = [
        sim.stat(alleleFreq=0),
        sim.pyExec('traj.append(alleleFreq[0][1])'),
    ],
    gen=5
)
# print trajectory
print ', '.join(['%.3f' % x for x in simu.dvars(0).traj])
#end_file

#begin_file log/pyEval.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(1000, loci=1,
    infoFields=['mother_idx', 'father_idx']),
    sim.randomMating())
simu.evolve(
    initOps = sim.initSex(),
    duringOps = sim.parentsTagger(),
    postOps = [
        sim.stat(alleleFreq=0),
        sim.pyEval(r'"gen %d, #father %d, #mother %d\n"' \
            ' % (gen, numFather, numMother)',
            stmts="numFather = len(set(pop.indInfo('father_idx')))\n"
                "numMother = len(set(pop.indInfo('mother_idx')))",
            exposePop='pop')
    ],
    gen=3
)
#end_file

#begin_file log/infoEval.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(20, loci=1, infoFields='a')
pop.setVirtualSplitter(sim.infoSplitter('a', cutoff=[3]))
sim.InitByFreq(pop, [0.2, 0.8])
pop.setIndInfo([random.randint(2, 5) for x in range(20)], 'a')
sim.InfoEval(pop, 'a', subPops=[(0, 0)]);print
sim.InfoEval(pop, 'ind.allele(0, 0)', exposeInd='ind');print
# use sim.population variables
pop.dvars().b = 5
sim.InfoEval(pop, '"%d " % (a+b)', usePopVars=True);print
#end_file

#begin_file log/infoExec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=1, infoFields=['a', 'b', 'c'])
sim.InitByFreq(pop, [0.2, 0.8])
sim.InfoExec(pop, 'a=1')
print pop.indInfo('a')[:10]
sim.InfoExec(pop, 'b=ind.sex()', exposeInd='ind')
print pop.indInfo('b')[:10]
sim.InfoExec(pop, 'c=a+b')
print pop.indInfo('c')[:10]
pop.dvars().d = 5
sim.InfoExec(pop, 'c+=d', usePopVars=True)
print pop.indInfo('c')[:10]
#end_file

#begin_file log/migrateByProb.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=[1000]*3, infoFields='migrate_to'),
    sim.randomMating())
simu.evolve(
    initOps = sim.initSex(),
    preOps = sim.migrator(rate=[
            [0, 0.1, 0.1],
            [0, 0, 0.1],
            [0, 0.1, 0]
        ]), 
    postOps = [
        sim.stat(popSize=True),
        sim.pyEval('subPopSize'),
        sim.pyOutput('\n')
    ],
    gen = 5
)        
#end_file

#begin_file log/migrateByPropAndCount.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=[1000]*3, infoFields='migrate_to'),
    sim.randomMating())
simu.evolve(
    initOps = sim.initSex(),
    preOps = sim.migrator(rate=[[0.1], [0.2]],
            mode=sim.ByProportion,
            subPops=[1, 2],
            toSubPops=[3]),
    postOps = [
        sim.stat(popSize=True),
        sim.pyEval('subPopSize'),
        sim.pyOutput('\n')
    ],
    gen = 5
)        
#
simu.evolve(
    preOps = sim.migrator(rate=[[50, 50], [100, 50]],
            mode=sim.ByCounts,
            subPops=[3, 2],
            toSubPops=[2, 1]),
    postOps = [
        sim.stat(popSize=True),
        sim.pyEval('subPopSize'),
        sim.pyOutput('\n')
    ],
    gen = 5
)        
#end_file

#begin_file log/migrateVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000]*2, infoFields='migrate_to')
pop.setVirtualSplitter(sim.sexSplitter())
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    # 500 males and 500 females
    initOps = sim.initSex(sex=[sim.Male, sim.Female]),
    preOps = [
        sim.migrator(rate=[
            [0, 0.10],
            [0, 0.05],
            ],
            mode = sim.ByProportion,
            subPops=[(0, 0), (0, 1)]),
        sim.stat(popSize=True, numOfMale=True, vars='numOfMale_sp'),
        sim.pyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMale'], subPopSize[0], "
            "subPop[1]['numOfMale'], subPopSize[1])"),
    ],
    postOps = [
        sim.stat(popSize=True, numOfMale=True, vars='numOfMale_sp'),
        sim.pyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMale'], subPopSize[0], "
            "subPop[1]['numOfMale'], subPopSize[1])"),
    ],
    gen = 2
)   
#end_file

#begin_file log/manualMigration.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population([10]*2, infoFields='migrate_to')
pop.setIndInfo([0, 1, 2, 3]*5, 'migrate_to')
sim.Migrate(pop, mode=sim.ByIndInfo)
pop.subPopSizes()
#end_file

#begin_file log/splitBySize.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(1000), sim.randomSelection())
simu.evolve(
    preOps = [
        sim.splitSubPops(subPops=0, sizes=[300, 300, 400], at=2),
        sim.stat(popSize=True),
        sim.pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 4
)
#end_file

#begin_file log/splitByProp.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def demo(gen, pop=None):
    if gen < 2:
        return 1000 + 100 * gen
    else:
        return [x + 50 * gen for x in pop.subPopSizes()]

simu = sim.simulator(sim.population(1000),
    sim.randomSelection(subPopSize=demo))
simu.evolve(
    preOps = [
        sim.splitSubPops(subPops=0, proportions=[.5]*2, at=2),
        sim.stat(popSize=True),
        sim.pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 4
)
#end_file

#begin_file log/splitByInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population([1000]*3, subPopNames=['a', 'b', 'c'], infoFields='x')
pop.setIndInfo([random.randint(0, 3) for x in range(1000)], 'x')
print pop.subPopSizes()
print pop.subPopNames()
sim.SplitSubPops(pop, subPops=[0, 2], infoFields=['x'])
print pop.subPopSizes()
print pop.subPopNames()
#end_file

#begin_file log/mergeSubPops.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population([500]*2),
    sim.randomSelection())
simu.evolve(
    preOps = [
        sim.mergeSubPops(subPops=[0, 1], at=3),
        sim.stat(popSize=True),
        sim.pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 5
)
#end_file

#begin_file log/resizeSubPops.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population([500]*2),
    sim.randomSelection())
simu.evolve(
    preOps = [
        sim.resizeSubPops(proportions=(1.5, 2), at=3),
        sim.stat(popSize=True),
        sim.pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 5
)
#end_file

#begin_file log/recRate.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=[1000], loci=[100]),
    sim.randomMating(ops = [
        sim.recombinator(rates=0.01, reps=0),
        sim.recombinator(rates=[0.01]*10, loci=range(50, 60), reps=1),
    ]), rep=2)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([0]*100 + [1]*100)
    ],
    postOps = [
        sim.stat(LD=[[40, 55], [60, 70]]),
        sim.pyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[60][70])'),
        sim.pyOutput('\n', reps=-1)
    ],
    gen = 5
)
#end_file

#begin_file log/recIntensity.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=[1000], loci=3, lociPos=[0, 1, 1.1]),
    sim.randomMating(ops=sim.recombinator(intensity=0.01)))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([0]*3 + [1]*3)
    ],
    postOps = [
        sim.stat(LD=[[0, 1], [1, 2]]),
        sim.pyEval(r'"%.3f\t%.3f\n" % (LD_prime[0][1], LD_prime[1][2])', step=10)
    ],
    gen = 50
)
#end_file

#begin_file log/conversion.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=[1000], loci=[100]),
    sim.randomMating(ops=[
        sim.recombinator(rates=0.01, loci=50, reps=0),
        sim.recombinator(rates=0.01, loci=50, reps=1, convMode=(sim.NumMarkers, 1, 10)),
    ]), rep=2)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([0]*100 + [1]*100)
    ],
    postOps = [
        sim.stat(LD=[[40, 55], [40, 70]]),
        sim.pyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[40][70])'),
        sim.pyOutput('\n', reps=-1)
    ],
    gen = 5
)
#end_file

#begin_file log/trackRec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=[1000, 2000], infoFields='ind_id')
simu = sim.simulator(pop, sim.randomMating(ops = [
    sim.idTagger(),
    sim.recombinator(rates=0.001, output='>>rec.log', infoFields='ind_id')])
)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.idTagger(),
    ],
    gen = 5
)
rec = open('rec.log')
# print the first three lines of the log file
print ''.join(rec.readlines()[:4])
#begin_ignore
rec.close()
import os
os.remove('rec.log')
#end_ignore
#end_file

#begin_file log/matrixMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=[2000], loci=1),
    sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.2, 0.3, 0.5])
    ],
    preOps = sim.matrixMutator(rate = [
            [0, 1e-5, 1e-5],
            [1e-4, 0, 1e-4],
            [1e-3, 1e-3, 0]
        ]),
    postOps = [
        sim.stat(alleleFreq=0, step=100),
        sim.pyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(3)]) + '\n'",
            step=100),
    ],
    gen=1000
)
#end_file

#begin_file log/kamMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=1*3)
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = sim.initSex(),
    postOps = [
        sim.kamMutator(k=5, rates=[1e-2, 1e-3], loci=[0, 1]),
        sim.stat(alleleFreq=range(3), step=100),
        sim.pyEval(r"', '.join(['%.3f' % alleleFreq[x][0] for x in range(3)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/snpMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=[1, 1], infoFields='fitness')
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = sim.initSex(),
    preOps = [
        sim.snpMutator(u=0.001),
        sim.maSelector(loci=0, fitness=[1, 0.99, 0.98]),
    ],
    postOps = [
        sim.stat(alleleFreq=[0, 1], step=100),
        sim.pyEval(r"'%.3f\t%.3f\n' % (alleleFreq[0][1], alleleFreq[1][1])",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/acgtMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=1,
    alleleNames=['A', 'C', 'G', 'T'])
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([.1, .1, .1, .7])
    ],
    preOps = [
        sim.acgtMutator(rate=[1e-4, 0.5], model='K80'),
        sim.stat(alleleFreq=0, step=100),
        sim.pyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(4)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/smmMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=1000, loci=[1, 1]), sim.randomMating())
simu.evolve(
    # all start from allele 50
    initOps = [
        sim.initSex(),
        sim.initByFreq( [0]*50 + [1])
    ],
    preOps = [
        sim.smmMutator(rates=1e-3, loci=0),
        sim.smmMutator(rates=1e-3, incProb=0.6, loci=1,
            mutStep=(sim.GeometricDistribution, 0.2)),
    ],
    gen=100
)
# count the average number tandem repeats at both loci
cnt0 = cnt1 = 0
for ind in simu.population(0).individuals():
    cnt0 += ind.allele(0, 0) + ind.allele(0, 1)
    cnt1 += ind.allele(1, 0) + ind.allele(1, 1)

print 'Average number of repeats at two loci are %.2f and %.2f.' % \
    (cnt0/2000., cnt1/2000.)
#end_file

#begin_file log/pyMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
def incAllele(allele):
    return allele + random.randint(1, 5)

simu = sim.simulator(sim.population(size=1000, loci=[20]),
    sim.randomMating())
simu.evolve(
    initOps = sim.initSex(),
    postOps = sim.pyMutator(func=incAllele, rates=[1e-4, 1e-3],
            loci=[2, 10]),
    gen = 1000
)
# count the average number tandem repeats at both loci
def avgAllele(pop, loc):
    ret = 0
    for ind in pop.individuals():
        ret += ind.allele(loc, 0) + ind.allele(loc, 1)
    return ret / (pop.popSize() * 2.)

pop = simu.population(0)
print 'Average number of repeats at two loci are %.2f and %.2f.' % \
    (avgAllele(pop, 2), avgAllele(pop, 10))
#end_file

#begin_file log/mixedMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(5000, loci=[1, 1]),
    sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([50, 50])
    ],
    preOps = [
        # the first locus uses a pure stepwise mutation model
        sim.smmMutator(rates=0.001, loci=0),
        # the second locus uses a mixed model
        sim.mixedMutator(rates=0.001, loci=1, mutators=[        
            sim.kamMutator(rates=1, k=100),
            sim.smmMutator(rates=1)
        ], prob=[0.1, 0.9])],
    gen = 20
)
# what alleles are there?
geno0 = []
geno1 = []
for ind in simu.population(0).individuals():
    geno0.extend([ind.allele(0, 0), ind.allele(0, 1)])
    geno1.extend([ind.allele(1, 0), ind.allele(1, 1)])

print 'Locus 0 has alleles', ', '.join([str(x) for x in set(geno0)])
print 'Locus 1 has alleles', ', '.join([str(x) for x in set(geno1)])
#end_file

#begin_file log/contextMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(5000, loci=[3, 3]),
    sim.randomMating())
simu.evolve(
    # initialize locus by 0, 0, 0, 1, 0, 1
    initOps = [
        sim.initSex(),
        sim.initByValue([1, 1], loci=[3, 5])
    ],
    preOps = [
        sim.contextMutator(mutators=[
            sim.snpMutator(u=0.1),
            sim.snpMutator(u=1),
            ],
            contexts=[(0, 0), (1, 1)],
            loci=[1, 4],
            rates=0.01
        ),
        sim.stat(alleleFreq=[1, 4], step=5),
        sim.pyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    gen = 20
)
#end_file

#begin_file log/pyContextMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
simu = sim.simulator(sim.population(5000, loci=[3, 3]),
    sim.randomMating())
def contextMut(allele, context):
    if context == [0, 0]:
        if allele == 0 and random.random() < 0.1:
            return 1
    elif context == [1, 1]:
        if allele == 0:
            return 1
    # do not mutate
    return allele

simu.evolve(
    # initialize locus by 0, 0, 0, 1, 0, 1
    initOps = [
        sim.initSex(),
        sim.initByValue([1, 1], loci=[3, 5])
    ],
    preOps = [
        sim.pyMutator(func=contextMut, context=1,
            loci=[1, 4],  rates=0.01
        ),
        #sim.snpMutator(u=0.01, v= 0.01, loci=[1, 4]),
        sim.stat(alleleFreq=[1, 4], step=5),
        sim.pyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    gen = 20
)
#end_file

#begin_file log/pointMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=1, infoFields='fitness')
simu = sim.simulator(pop, sim.randomSelection())
simu.evolve(
    initOps = sim.pyOutput('Introducing alleles at generation'),
    preOps = sim.maSelector(loci=0, wildtype=0, fitness=[1, 1.05, 1.1]),
    postOps = [
        sim.stat(alleleFreq=0),
        sim.ifElse('alleleNum[0][1] == 0', ifOps=[
            sim.pyEval(r"' %d' % gen"),
            sim.pointMutator(inds=0, loci=0, allele=1),
        ]),
        sim.ifElse('alleleFreq[0][1] > 0.05', ifOps=[
            sim.pyEval(r"'.\nTerminate at generation %d at allele freq %.3f.\n'" +
                " % (gen, alleleFreq[0][1])"),
            sim.terminateIf('True'),
        ])
    ],
)
#end_file

#begin_file log/mutatorVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
def fragileX(geno, gen):
    '''A disease model where an individual has increased risk of 
    affected if the number of tandem repeats exceed 75.
    '''
    # Alleles A1, A2.
    maxRep = max(geno)
    if maxRep < 50:
        return 0
    else:
        # individuals with allele >= 70 will surely be affected
        return min(1, (maxRep - 50)*0.05)

def avgAllele(pop):
    'Get average allele by affection sim.status.'
    sim.Stat(pop, alleleFreq=(0,1), subPops=[(0,0), (0,1)],
        numOfAffected=True, vars=['alleleNum', 'alleleNum_sp'])
    avg = []
    for alleleNum in [\
            pop.dvars((0,0)).alleleNum[0],  # first locus, unaffected
            pop.dvars((0,1)).alleleNum[0],  # first locus, affected
            pop.dvars().alleleNum[1],       # second locus, overall
        ]:
        alleleSum = numAllele = 0
        for idx,cnt in enumerate(alleleNum):
            alleleSum += idx * cnt
            numAllele += cnt
        if numAllele == 0:
            avg.append(0)
        else:
            avg.append(alleleSum * 1.0 /numAllele)
    # unaffected, affected, loc2
    pop.dvars().avgAllele = avg
    return True

pop = sim.population(10000, loci=[1, 1])
pop.setVirtualSplitter(sim.affectionSplitter())
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([50, 50])
    ],
    postOps = [
        # determine affection sim.status for each offspring (duringMating)
        sim.pyPenetrance(func=fragileX, loci=0),
        # unaffected offspring, mutation rate is high to save some time
        sim.smmMutator(rates=1e-3, loci=1),
        # unaffected offspring, mutation rate is high to save some time
        sim.smmMutator(rates=1e-3, loci=0, subPops=[(0, 0)]),
        # affected offspring have high probability of mutating upward
        sim.smmMutator(rates=1e-2, loci=0, subPops=[(0, 1)],
           incProb=0.7, mutStep=3),
        # number of affected
        sim.pyOperator(func=avgAllele, step=20),
        sim.pyEval(r"'Gen: %3d #Aff: %d AvgRepeat: %.2f (unaff), %.2f (aff), %.2f (unrelated)\n'"
            + " % (gen, numOfAffected, avgAllele[0], avgAllele[1], avgAllele[2])",
            step=20),
    ],
    gen = 101
)
#end_file

#begin_file log/alleleMapping.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=1)
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0]*4 + [0.1, 0.2, 0.3, 0.4])
    ],
    postOps = [
        sim.kamMutator(k=4, rates=1e-4, mapIn=[0]*4 + range(4),
            mapOut=[4, 5, 6, 7]),
        sim.stat(alleleFreq=0, step=100),
        sim.pyEval(r"', '.join(['%.2f' % alleleFreq[0][x] for x in range(8)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/infiniteSites.py
import simuOpt
simuOpt.setOptions(alleleType='long')
#begin_ignore
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore

def infSitesMutate(pop, param):
    '''Apply an infinite mutation model'''
    (startPos, endPos, rate) = param
    # for each individual
    for ind in pop.individuals():
        # for each homologous copy of chromosomes
        for p in range(2):
            # using a geometric distribution to determine
            # the first mutation location
            loc = sim.GetRNG().randGeometric(rate)
            # if a mutation happens, record the mutated location
            if startPos + loc < endPos:
                try:
                    # find the first non-zero location
                    idx = ind.genotype(p).index(0)
                    # record mutation here
                    ind.setAllele(startPos + loc, idx, ploidy=p)
                except:
                    print 'Warning: more than %d mutations have accumulated' % pop.totNumLoci()
                    pass
    return True

pop = sim.population(size=[2000], loci=[100])
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = sim.initSex(),
    preOps = [
        # recombine in a 10Mb region at rate 1e-8
        sim.pyOperator(func=infSitesMutate, param=(1, 10000000, 1e-8)),
    ],
    gen = 100
)
# now, we get a sim.population. Let us have a look at the 'alleles'.
pop = simu.extract(0)
# print the first five mutation locations
print pop.individual(0).genotype()[:5]
# how many alleles are there (does not count 0)?
print len(set(pop.genotype())) - 1
# Allele count a simple count of alleles.
cnt = {}
for allele in pop.genotype():
    if allele == 0:
        continue
    if cnt.has_key(allele):
        cnt[allele] += 1
    else:
        cnt[allele] = 1

# highest allele frequency?
print max(cnt.values()) *0.5 / pop.popSize()
#end_file


#begin_file log/statSuffix.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population([5000]*3, loci=5), sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    postOps = [
        sim.stat(structure=range(5), subPops=(0, 1), suffix='_01', step=40),
        sim.stat(structure=range(5), subPops=(1, 2), suffix='_12', step=40),
        sim.stat(structure=range(5), subPops=(0, 2), suffix='_02', step=40),
        sim.stat(structure=range(5), step=40),
        sim.pyEval(r"'Fst=%.3f (pairwise: %.3f %.3f %.3f)\n' % (F_st, F_st_01, F_st_12, F_st_02)",
            step=40),
    ],
    gen = 200
)
#end_file

#begin_file log/statCount.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(10000, loci=1)
pop.setVirtualSplitter(sim.combinedSplitter(
    [sim.sexSplitter(), sim.affectionSplitter()]))
sim.InitSex(pop)
sim.InitByFreq(pop, [0.2, 0.8])
sim.MaPenetrance(pop, loci=0, penetrance=[0.1, 0.2, 0.5])
# Count sim.population size
sim.Stat(pop, popSize=True, subPops=[(0, 0), (0, 2)])
# popSize is the size of two VSPs, does not equal to total sim.population size.
# Because two VSPs overlap (all males and all unaffected), popSize can be
# greater than real sim.population size.
print pop.dvars().subPopSize, pop.dvars().popSize
# print popSize of each virtual subpopulation.
sim.Stat(pop, popSize=True, subPops=[(0, 0), (0, 2)], vars='popSize_sp')
# Note the two ways to access variable in (virtual) subpopulations.
print pop.dvars((0,0)).popSize, pop.dvars().subPop[(0,2)]['popSize']
# Count number of male (should be the same as the size of VSP (0,0).
sim.Stat(pop, numOfMale=True)
print pop.dvars().numOfMale
# Count the number of affected and unaffected male individual
sim.Stat(pop, numOfMale=True, subPops=[(0, 2), (0, 3)], vars='numOfMale_sp')
print pop.dvars((0,2)).numOfMale, pop.dvars((0,3)).numOfMale
# or number of affected male and females
sim.Stat(pop, numOfAffected=True, subPops=[(0, 0), (0, 1)], vars='numOfAffected_sp')
print pop.dvars((0,0)).numOfAffected, pop.dvars((0,1)).numOfAffected
# These can also be done using a sim.productSplitter...
pop.setVirtualSplitter(sim.productSplitter(
    [sim.sexSplitter(), sim.affectionSplitter()]))
sim.Stat(pop, popSize=True, subPops=[(0, x) for x in range(4)])
# counts for male unaffected, male affected, female unaffected and female affected
print pop.dvars().subPopSize
#end_file

#begin_file log/statAlleleFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(10000, loci=1)
pop.setVirtualSplitter(sim.affectionSplitter())
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(loci=0, alleleFreq=[0.8, 0.2])
    ],
    postOps = [
        sim.maPenetrance(penetrance=[0.1, 0.4, 0.6], loci=0),
        sim.stat(alleleFreq=0, subPops=[(0, 0), (0, 1)],
            vars=['alleleFreq', 'alleleFreq_sp']),
        sim.pyEval(r"'Gen: %d, freq: %.2f, freq (aff): %.2f, freq (unaff): %.2f\n' % " + \
            "(gen, alleleFreq[0][1], subPop[(0,1)]['alleleFreq'][0][1]," + \
            "subPop[(0,0)]['alleleFreq'][0][1])"),
    ],
    gen = 5
)
#end_file

#begin_file log/statGenoFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=[1, 1, 1], chromTypes=[sim.Autosome, sim.ChromosomeX, sim.ChromosomeY])
sim.InitByFreq(pop, [0.01, 0.05, 0.94])
sim.Stat(pop, genoFreq=[0, 1])
print 'Available genotypes on autosome:', pop.dvars().genoFreq[0].keys()
for i in range(3):
    for j in range(3):
        print '%d-%d: %.3f' % (i, j, pop.dvars().genoFreq[0][(i,j)])

print 'Genotype frequency on chromosome X:\n', \
    '\n'.join(['%s: %.3f' % (x,y) for x,y in pop.dvars().genoFreq[1].iteritems()])
#end_file

#begin_file log/statHeteroFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=1)
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    postOps = [
        sim.stat(heteroFreq=0, step=10),
        sim.pyEval(r"'Gen: %d, HeteroFreq: %.2f\n' % (gen, heteroFreq[0])", step=20)
    ],
    gen = 100
)
#end_file

#begin_file log/statHaploFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import ViewVars
pop = sim.population(100, loci=3)
sim.InitByFreq(pop, [0.2, 0.4, 0.4], loci=0)
sim.InitByFreq(pop, [0.2, 0.8], loci=2)
sim.Stat(pop, genoFreq=[0, 1, 2], haploFreq=[0, 1, 2],
    vars=['genoNum', 'haploFreq'])
ViewVars(pop.vars(), gui=False)
#end_file

#begin_file log/statInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population([500], infoFields='anc')
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sim.infoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
#
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        # anc is 0 or 1
        sim.initInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    duringOps = sim.inheritTagger(mode=sim.Mean, infoFields='anc'),
    postOps = [
        sim.stat(popSize=True, meanOfInfo='anc', varOfInfo='anc',
            subPops=[(0,x) for x in range(5)]),
        sim.pyEval(r"'Anc: %.2f (%.2f), #inds: %s\n' %" + \
            "(meanOfInfo['anc'], varOfInfo['anc'], " + \
            "', '.join(['%4d' % x for x in subPopSize]))")
    ],
    gen = 5,
)
#end_file

#begin_file log/statLD.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population([1000]*2, loci=3)
sim.InitByFreq(pop, [0.2, 0.8], subPops=0)
sim.InitByFreq(pop, [0.8, 0.2], subPops=1)
sim.Stat(pop, LD=[[0, 1, 0, 0], [1, 2]],
    vars=['LD', 'LD_prime', 'R2', 'LD_ChiSq', 'LD_ChiSq_p', 'CramerV',
        'LD_prime_sp', 'LD_ChiSq_p_sp'])
from pprint import pprint
pprint(pop.vars())
#end_file

#begin_file log/statAssociation.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import *
from simuPOP.sampling import CaseControlSample
def assoTest(pop):
    'Draw case-control sample and apply association tests'
    sample = CaseControlSample(pop, cases=500, controls=500)[0]
    sim.Stat(sample, association=(0, 2), vars=['Allele_ChiSq_p', 'Geno_ChiSq_p', 'Armitage_p'])
    print 'Allele test: %.2e, %.2e, Geno test: %.2e, %.2e, Trend test: %.2e, %.2e' \
        % (sample.dvars().Allele_ChiSq_p[0], sample.dvars().Allele_ChiSq_p[2],
        sample.dvars().Geno_ChiSq_p[0], sample.dvars().Geno_ChiSq_p[2],
        sample.dvars().Armitage_p[0], sample.dvars().Armitage_p[2])
    return True

simu = sim.simulator(sim.population(size=100000, loci=3),
    sim.randomMating(ops=sim.recombinator(loci=[0, 1], rates=[0.01, 0.005])))
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([[0]*3, [1]*3], proportions=[0.5, 0.5])
    ],
    postOps = [
        sim.maPenetrance(loci=1, penetrance=[0.1, 0.2, 0.4]),
        sim.pyOperator(func=assoTest, step=20),
    ],
    gen = 100
)
#end_file

#begin_file log/statStructure.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import MigrIslandRates
simu = sim.simulator(sim.population([5000]*3, loci=10, infoFields='migrate_to'),
    sim.randomMating(), rep=2)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    preOps = sim.migrator(rate=MigrIslandRates(0.01, 3), reps=1),
    postOps = [
        sim.stat(structure=range(10), step=40),
        sim.pyEval("'Fst=%.3f (rep=%d without migration) ' % (F_st, rep)", step=40, reps=0),
        sim.pyEval("'Fst=%.3f (rep=%d with migration) ' % (F_st, rep)", step=40, reps=1),
        sim.pyOutput('\n', reps=-1, step=40)
    ],
    gen = 200
)
#end_file

#begin_file log/statHWE.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population([1000], loci=1), sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([[0,0], [0, 1], [1,1]], proportions=[0.4, 0.4, 0.2])
    ],
    preOps = [
        sim.stat(HWE=0, genoFreq=0),
        sim.pyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
            'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])'),
    ],
    postOps = [
        sim.stat(HWE=0, genoFreq=0),
        sim.pyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
            'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])'),
    ],
    gen = 1
)
#end_file

#begin_file log/inheritTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000]*10, loci=1, infoFields='x')
# tag the first individual of each subpopulation.
for sp in range(pop.numSubPop()):
    pop.individual(0, sp).x = 1

simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = sim.initSex(),
    duringOps = sim.inheritTagger(mode=sim.Maximum, infoFields='x'),
    postOps = [
        sim.stat(sumOfInfo='x', vars=['sumOfInfo_sp']),
        sim.pyEval(r'", ".join(["%3d" % subPop[i]["sumOfInfo"]["x"] for i in range(10)])+"\n"'),
    ],
    gen = 5
)
#end_file

#begin_file log/summaryTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
sim.idTagger().reset(0)
#end_ignore
simu = sim.simulator(
    sim.population(1000, loci=1, infoFields=['fitness', 'avgFitness']),
    sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5]),
    ],
    preOps = sim.maSelector(loci=0, wildtype=0, fitness=[1, 0.99, 0.95]),
    duringOps = sim.summaryTagger(mode=sim.Mean, infoFields=['fitness', 'avgFitness']),
    postOps = [
        sim.stat(alleleFreq=0, meanOfInfo='avgFitness', step=10),
        sim.pyEval(r"'gen %d: allele freq: %.3f, average fitness of parents: %.3f\n' % "
            "(gen, alleleFreq[0][1], meanOfInfo['avgFitness'])", step=10)
    ],
    gen = 50,
)
#end_file


#begin_file log/idTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
sim.idTagger().reset(0)
#end_ignore
simu = sim.simulator(
    sim.population(10, infoFields='ind_id', ancGen=1),
    sim.randomSelection())
simu.evolve(
    initOps = sim.idTagger(),
    duringOps = sim.idTagger(),
    gen = 1
)
pop = simu.extract(0)
print [int(ind.ind_id) for ind in pop.individuals()]
pop.useAncestralGen(1)
print [int(ind.ind_id) for ind in pop.individuals()]
sim.TagID(pop) # re-assign ID
print [int(ind.ind_id) for ind in pop.individuals()]
#end_file

#begin_file log/pedigreeTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
sim.idTagger().reset(0)
#end_ignore
simu = sim.simulator(
    sim.population(100, infoFields=['ind_id', 'father_id', 'mother_id']),
    sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.idTagger()
    ],
    duringOps = [
        sim.idTagger(),
        sim.pedigreeTagger(output=">>sim.pedigree.txt")
    ],
    gen = 100
)
ped = open('sim.pedigree.txt')
print ''.join(ped.readlines()[100:105])
#begin_ignore
ped.close()
import os
os.remove('sim.pedigree.txt')
#end_ignore
#end_file

#begin_file log/pyTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
def randomMove(values):
    '''Pass parental information fields to offspring'''
    x1, y1, x2, y2 = values
    # shift right with high concentration of alleles... 
    x = random.normalvariate((x1+x2)/2., 0.1)
    y = random.normalvariate((y1+y2)/2., 0.1)
    return (x, y)

pop = sim.population(1000, loci=[1], infoFields=['x', 'y'])
pop.setVirtualSplitter(sim.genotypeSplitter(loci=0, alleles=[[0, 0], [0,1], [1, 1]]))
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5]),
        sim.initInfo(random.random, infoFields=['x', 'y'])
    ],
    duringOps = sim.pyTagger(func=randomMove, infoFields=['x', 'y']),
    postOps = [
        sim.stat(minOfInfo='x', maxOfInfo='x'),
        sim.pyEval(r"'Range of x: %.2f, %.2f\n' % (minOfInfo['x'], maxOfInfo['x'])")
    ],
    gen = 5
)

#end_file

#begin_file log/otherTagging.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=[1], infoFields=['aff', 'numAff'])
# define virtual subpopulations by affection sim.status
pop.setVirtualSplitter(sim.affectionSplitter())
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5]),
    ],
    preOps = [
        # get affection sim.status for parents
        sim.maPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4]),
        # set 'aff' of parents
        sim.infoExec('aff = ind.affected()', exposeInd='ind'),
    ],
        # get number of affected parents for each offspring and store in numAff
    duringOps = sim.summaryTagger(mode=sim.Summation, infoFields=['aff', 'numAff']),
    postOps = [
        # get affection sim.status for offspring
        sim.maPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4]),
        # calculate mean 'numAff' of offspring, for unaffected and affected subpopulations.
        sim.stat(meanOfInfo='numAff', subPops=[(0,0), (0,1)], vars=['meanOfInfo_sp']),
        # print mean number of affected parents for unaffected and affected offspring.
        sim.pyEval(r"'sim.Mean number of affected parents: %.2f (unaff), %.2f (aff)\n' % "
            "(subPop[(0,0)]['meanOfInfo']['numAff'], subPop[(0,1)]['meanOfInfo']['numAff'])")
    ],
    gen = 5
)

#end_file



#begin_file log/mapPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=2000, loci=2)
sim.InitByFreq(pop, [.2, .8])
sim.MapPenetrance(pop, loci=0,
    penetrance={(0,0):0, (0,1):.2, (1,1):.3})
sim.Stat(pop, genoFreq=0, numOfAffected=1, vars='genoNum')
# number of affected individuals
pop.dvars().numOfAffected
# which should be roughly (#01 + #10) * 0.2 + #11 * 0.3
(pop.dvars().genoNum[0][(0,1)] + pop.dvars().genoNum[0][(1,0)]) * 0.2 \
+ pop.dvars().genoNum[0][(1,1)] * 0.3
#end_file

#begin_file log/maPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(5000, loci=3)
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.9] + [0.02]*5)
    ],
    postOps = [
        sim.maPenetrance(loci=0, penetrance=(0.01, 0.2, 0.3)),
        sim.stat(numOfAffected=True, vars='propOfAffected'),
        sim.pyEval(r"'Gen: %d Prevalence: %.1f%%\n' % (gen, propOfAffected*100)"),
    ],
    gen = 5
)
#end_file

#begin_file log/mlPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(5000, loci=3)
sim.InitByFreq(pop, [0.2]*5)
# the multi-loci penetrance
sim.MlPenetrance(pop, mode=sim.Multiplicative,
    ops = [sim.maPenetrance(loci=loc,
        penetrance=[0, 0.3, 0.6]) for loc in range(3)])
# count the number of affected individuals.
sim.Stat(pop, numOfAffected=True)
pop.dvars().numOfAffected
#end_file

#begin_file log/pyPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=2000, loci=[1]*2, infoFields=['p', 'smoking'])
pop.setVirtualSplitter(sim.infoSplitter(field='smoking', values=[0,1]))
simu = sim.simulator(pop, sim.randomMating())
# the second parameter gen can be used for varying selection pressure
def penet(arr, smoking, gen):
    #     BB     Bb      bb
    # AA  0.01   0.01    0.01
    # Aa  0.01   0.03    0.03
    # aa  0.01   0.03    0.05
    #
    # arr is (A1 A2 B1 B2)
    if arr[0] + arr[1] == 1 and arr[2] + arr[3] != 0:
        v = 0.03   # case of AaBb
    elif arr[0] + arr[1] == 2 and arr[2] + arr[3] == 1:
        v = 0.03   # case of aaBb
    elif arr[0] + arr[1] ==2 and arr[2] + arr[3] == 2:
        v = 0.05   # case of aabb
    else:                
        v = 0.01   # other cases
    if smoking[0]:
        return v * 2
    else:
        return v

simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(alleleFreq=[.5, .5]),
        sim.pyOutput('Calculate prevalence in smoker and non-smokers'),
    ],
    postOps = [
        # set smoking status randomly
        sim.initInfo(lambda : random.randint(0,1), infoFields='smoking'),
        # assign affection status
        sim.pyPenetrance(loci=[0, 1], func=penet, paramFields='smoking'),
        sim.stat(numOfAffected=True, subPops=[(0,0),(0,1)], 
            vars='propOfAffected_sp', step=20),
        sim.pyEval(r"'Non-smoker: %.2f%%\tSmoker: %.2f%%\n' % "
            "(subPop[(0,0)]['propOfAffected']*100, subPop[(0,1)]['propOfAffected']*100)",
            step=20)
    ],
    gen = 50
)

#end_file



#begin_file log/pyQuanTrait.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=5000, loci=2, infoFields=['qtrait1', 'qtrait2', 'age'])
pop.setVirtualSplitter(sim.infoSplitter(field='age', cutoff=[40]))
simu = sim.simulator(pop, sim.randomMating())
def qtrait(geno, age, gen):
    'Return two traits that depends on genotype and age'
    return random.normalvariate(age[0] * sum(geno), 10), random.randint(0, 10*sum(geno))

simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.2, 0.8]),
    ],
    postOps = [
        # use random age for simplicity
        sim.initInfo(lambda:random.randint(20, 75), infoFields='age'),
        sim.pyQuanTrait(loci=(0,1), func=qtrait, paramFields='age',
            infoFields=['qtrait1', 'qtrait2']),
        sim.stat(meanOfInfo=['qtrait1'], subPops=[(0,0), (0,1)],
            vars='meanOfInfo_sp'),
        sim.pyEval(r"'Mean of trait1: %.3f (age < 40), %.3f (age >=40)\n' % "
            "(subPop[(0,0)]['meanOfInfo']['qtrait1'], subPop[(0,1)]['meanOfInfo']['qtrait1'])"),
    ],
    gen = 5
)

#end_file

#begin_file log/selectParents.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(4000, loci=1, infoFields='fitness')
simu = sim.simulator(pop, sim.randomMating(), rep=3)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    preOps = sim.mapSelector(loci=0, fitness={(0,0):1, (0,1):0.98, (1,1):0.97}),
    postOps = [
        sim.stat(alleleFreq=0, step=10),
        sim.pyEval("'Gen:%3d ' % gen", reps=0, step=10),
        sim.pyEval(r"'%.3f\t' % alleleFreq[0][1]", step=10),
        sim.pyOutput('\n', reps=-1, step=10)
    ],
    gen = 50
)
#end_file

#begin_file log/selectOffspring.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(10000, loci=1)
simu = sim.simulator(pop, sim.randomMating(), rep=3)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.5, 0.5])
    ],
    duringOps = sim.mapSelector(loci=0, fitness={(0,0):1, (0,1):0.98, (1,1):0.97}),
    postOps = [
        sim.stat(alleleFreq=0, step=10),
        sim.pyEval("'Gen:%3d ' % gen", reps=0, step=10),
        sim.pyEval(r"'%.3f\t' % alleleFreq[0][1]", step=10),
        sim.pyOutput('\n', reps=-1, step=10)
    ],
    gen = 50
)
#end_file

#begin_file log/mapSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=1000, loci=1, infoFields='fitness'),
    sim.randomMating())
s1 = .1
s2 = .2
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(alleleFreq=[.2, .8])
    ],
    preOps = sim.mapSelector(loci=0, fitness={(0,0):1-s1, (0,1):1, (1,1):1-s2}),
    postOps = [
        sim.stat(alleleFreq=0),
        sim.pyEval(r"'%.4f\n' % alleleFreq[0][0]", step=100)
    ],
    gen=301
)

#end_file


#begin_file log/maSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=1000, loci=1, infoFields='fitness'),
    sim.randomMating())
s1 = .1
s2 = .2
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(alleleFreq=[.2] * 5)
    ],
    preOps = sim.maSelector(loci=0, fitness=[1-s1, 1, 1-s2]),
    postOps = [
        sim.stat(alleleFreq=0),
        sim.pyEval(r"'%.4f\n' % alleleFreq[0][0]", step=100)
    ],
    gen = 301)
#end_file

#begin_file log/maSelectorHaploid.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=10000, ploidy=1, loci=[1,1], infoFields='fitness'),
    sim.randomSelection())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(alleleFreq=[.5, .5])
    ],
    # fitness values for AB, Ab, aB and ab
    preOps = sim.maSelector(loci=[0,1], fitness=[1, 1, 1, 0.95]),
    postOps = [
        sim.stat(haploFreq=[0, 1], step=25),
        sim.pyEval(r"'%.3f\t%.3f\t%.3f\t%.3f\n' % (haploFreq[(0,1)][(0,0)],"
                "haploFreq[(0,1)][(0,1)], haploFreq[(0,1)][(1,0)],"
                "haploFreq[(0,1)][(1,1)])", step=25)
    ],
    gen = 100
)
#end_file

#begin_file log/mlSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=10000, loci=2, infoFields='fitness'),
    sim.randomMating(), rep=3)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(alleleFreq=[.5, .5])
    ],
    preOps = [
        sim.mlSelector([
            sim.mapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}),
            sim.mapSelector(loci=1, fitness={(0,0):1, (0,1):0.9, (1,1):.8}),
            ], mode = sim.Additive, reps=0),
        sim.mapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}, reps=1),
        sim.mapSelector(loci=1, fitness={(0,0):1, (0,1):0.9, (1,1):.8}, reps=2)
    ],
    postOps = [
         sim.stat(alleleFreq=[0,1]),
         sim.pyEval(r"'REP %d:\t%.3f\t%.3f\t' % (rep, alleleFreq[0][1], alleleFreq[1][1])"),
         sim.pyOutput('\n', reps=-1),
    ],
    gen = 5
)
#end_file

#begin_file log/pySelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
import random
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=2000, loci=[1]*2, infoFields=['fitness', 'smoking']),
    sim.randomMating()
)
s1 = .02
s2 = .03
# the second parameter gen can be used for varying selection pressure
def sel(arr, smoking, gen):
    #     BB  Bb   bb
    # AA  1   1    1
    # Aa  1   1-s1 1-s2
    # aa  1   1    1-s2
    #
    # arr is (A1 A2 B1 B2)
    if arr[0] + arr[1] == 1 and arr[2] + arr[3] == 1:
        v = 1 - s1  # case of AaBb
    elif arr[2] + arr[3] == 2:
        v = 1 - s2  # case of ??bb
    else:                
        v = 1       # other cases
    if smoking[0]:
        return v * 0.9
    else:
        return v

simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(alleleFreq=[.5, .5])
    ],
    preOps = sim.pySelector(loci=[0, 1], func=sel, paramFields='smoking'),
    postOps = [
        # set smoking status randomly
        sim.initInfo(lambda : random.randint(0,1), infoFields='smoking'),
        sim.stat(alleleFreq=[0, 1], step=20),
        sim.pyEval(r"'%.4f\t%.4f\n' % (alleleFreq[0][1], alleleFreq[1][1])", step=20)
    ],
    gen = 50
)
#end_file


#begin_file log/peneSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=2000, loci=1, infoFields='fitness'),
    sim.randomMating()
)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(alleleFreq=[.5, .5])
    ],
    preOps = [
        sim.maPenetrance(loci=0, penetrance=[0.01, 0.1, 0.2]),
        sim.stat(numOfAffected=True, step=25, vars='propOfAffected'),
        sim.pyEval(r"'Percent of affected: %.3f\t' % propOfAffected", step=50),
        sim.infoExec('fitness = not ind.affected()', exposeInd='ind')
    ],
    postOps = [
        sim.stat(alleleFreq=0),
        sim.pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=50)
    ],
    gen=151
)
#end_file


#begin_file log/vspSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[5000, 5000], loci=1, infoFields='fitness')
pop.setVirtualSplitter(sim.sexSplitter())
simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq(alleleFreq=[.5, .5])
    ],
    preOps = [
        sim.maSelector(loci=0, fitness=[1, 1, 0.98], subPops=[(0,0), (1,1)]),
        sim.maSelector(loci=0, fitness=[1, 1, 1], subPops=[(0,1), (1,0)]),
    ],
    postOps = [
        sim.stat(alleleFreq=[0], subPops=[(0,0), (0,1), (1,0), (1,1)],
            vars='alleleFreq_sp', step=50),
        sim.pyEval(r"'%.4f\t%.4f\t%.4f\t%.4f\n' % "
            "tuple([subPop[x]['alleleFreq'][0][1] for x in ((0,0),(0,1),(1,0),(1,1))])",
            step=50)
    ],
    gen=151
)
#end_file


#begin_file log/forwardTrajectory.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import trajectory, ForwardTrajectory

traj = ForwardTrajectory(N=[2000, 4000], fitness=[1, 0.99, 0.98],
    beginGen=0, endGen=100, beginFreq=[0.2, 0.3],
    endFreq=[[0.1, 0.11], [0.2, 0.21]])
traj.plot('log/forwardTrajectory.png', plot_ylim=[0, 0.5], col_sp=['red', 'blue'],
    plot_main='Simulated trajectory (forward-time)')
simu = sim.simulator(
    sim.population(size=[2000, 4000], loci=10, infoFields='fitness'),
    sim.controlledRandomMating(
        ops=[sim.recombinator(rates=0.01)],
        loci=5, alleles=1, freqFunc=traj.func())
)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByFreq([0.8, 0.2], subPops=0),
        sim.initByFreq([0.7, 0.3], subPops=1),
        sim.pyOutput('Sp0: loc2\tloc5\tSp1: loc2\tloc5\n'),
    ],
    postOps = [
        sim.stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=20),
        sim.pyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
            "subPop[0]['alleleFreq'][5][1], subPop[1]['alleleFreq'][2][1],"
            "subPop[1]['alleleFreq'][5][1])", step=20)
    ],
    gen = 101
)
#end_file

#begin_file log/backTrajectory.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import trajectory, BackwardTrajectory
from math import exp
def Nt(gen, pop=None):
    'An exponential sim.population growth demographic model.'
    return int((10**4) * exp(.00115 * gen))

def fitness(gen, sp):
    'sim.Constant positive selection pressure.'
    return [1, 1.01, 1.02]

# simulate a trajectory backward in time, from generation 1000
traj = BackwardTrajectory(N=Nt, fitness=fitness, nLoci=2,
     endGen=1000, endFreq=[0.1, 0.2])
traj.plot('log/backTrajectory.png', plot_ylim=[0, 0.3], plot_xlim=[0, 1000],
    col_loc=['red', 'blue'], plot_main='Simulated trajectory (backward-time)')
print 'Trajectory simulated with length %s ' % len(traj.traj)
pop = sim.population(size=Nt(0), loci=[1]*2)
# save trajectory function in the sim.population's local namespace
# so that the sim.pyEval operator can access it.
pop.dvars().traj = traj.func()
simu = sim.simulator(pop, sim.controlledRandomMating(loci=[0, 1], alleles=[1, 1],
        subPopSize=Nt, freqFunc=traj.func()))
simu.evolve(
    initOps = [sim.initSex()],
    preOps = traj.mutators(loci=[0, 1]),
    postOps = [
        sim.stat(alleleFreq=[0, 1], begin=500, step=100),
        sim.pyEval(r"'%4d: %.3f (exp: %.3f), %.3f (exp: %.3f)\n' % (gen, alleleFreq[0][1],"
            "traj(gen)[0], alleleFreq[1][1], traj(gen)[1])",
            begin=500, step=100)
    ],
    gen=1001  # evolve 1001 generations to reach the end of generation 1000
)
#end_file

#begin_file log/varPlotter.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import varPlotter
pop = sim.population(size=1000, loci=2)
simu = sim.simulator(pop, sim.randomMating(ops=sim.recombinator(rates=0.01)), rep=3)
simu.evolve(
    initOps = [
        sim.initSex(),
        sim.initByValue([1, 2, 2, 1])
    ],
    postOps = [
        sim.stat(LD=[0, 1]),
        varPlotter('LD[0][1]', step=5, update=40, saveAs='log/rpy.png',
            legend=['Replicate %d' % x for x in range(3)],
            ylab='LD between marker 1 and 2',
            ylim=[0, 0.25], main='LD decay', lty_rep=[1, 2, 3],
        ),
    ],
    gen=100
)
#end_file

#begin_file log/rpyByRep.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import varPlotter
pop = sim.population(size=1000, loci=1*4)
simu = sim.simulator(pop, sim.randomMating(), rep=3)
simu.evolve(
    initOps = [sim.initSex()] +
        [sim.initByFreq([0.1*(x+1), 1-0.1*(x+1)], loci=x) for x in range(4)],
    postOps = [
        sim.stat(alleleFreq=range(4)),
        varPlotter('[alleleFreq[x][0] for x in range(4)]', byRep=True,
            update=10, saveAs='log/rpy_byRep.png',
            legend=['Locus %d' % x for x in range(4)],
            ylab='Allele frequency',
            ylim=[0, 1],
            main_rep=['Genetic drift, replicate %d' % x for x in range(3)],
        ),
    ],
    gen=100
)
#end_file

#begin_file log/rpyByDim.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import varPlotter
pop = sim.population(size=1000, loci=1*4)
simu = sim.simulator(pop, sim.randomMating(), rep=3)
def drawFrame(r, dim=None, **kwargs):
    '''Draw a frame around subplot dim. Parameter r is defined in the rpy
    module and is used for calling R functions. Parameter dim is the dimension
    index. Other parameters are ignored.
    '''
    r.axis(1)
    r.axis(2)
    r.grid()
    r.mtext({0:'A', 1:'B', 2:'C', 3:'D'}[dim], adj=1)

simu.evolve(
    initOps = [sim.initSex()]+
        [sim.initByFreq([0.1*(x+1), 1-0.1*(x+1)], loci=x) for x in range(4)],
    postOps = [
        sim.stat(alleleFreq=range(4)),
        varPlotter('[alleleFreq[x][0] for x in range(4)]', byDim=True,
            update=10, saveAs='log/rpy_byDim.png',
            legend=['Replicate %d' % x for x in range(3)],
            ylab='Allele frequency',
            ylim=[0, 1],
            main_dim=['Genetic drift, freq=%.1f' % ((x+1)*0.10) for x in range(4)],
            col_rep=['red', 'blue', 'black'],
            lty_rep=[1, 2, 3],
            # the default png dimension is 800x600
            dev_print_width=600, dev_print_height=500,
            # do not draw axes in r.plot, leaving the job to drawFrame
            plot_axes=False,
            # plot frame, grid etc after each r.plot call
            plotHook = drawFrame,
        ),
    ],
    gen=100
)
#end_file

#begin_file log/scatterPlotter.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import scatterPlotter
import random
pop = sim.population([500], infoFields=['x', 'y', 'anc'])
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sim.infoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
#
def passInfo(fields):
    'Parental fields will be passed as x1, y1, anc1, x2, y2, anc2'
    x1, y1, anc1, x2, y2, anc2 = fields
    anc = (anc1 + anc2)/2.
    x = (x1 + x2)/2 + random.normalvariate(anc - 0.5, 0.1)
    y = (y1 + y2)/2 + random.normalvariate(0, 0.1)
    return x, y, anc

simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        # random geographic location
        sim.initInfo(random.random, infoFields=['x', 'y']),
        # anc is 0 or 1
        sim.initInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    postOps = [
        scatterPlotter(['x', 'y'], 
            saveAs = 'log/scatterPlotter.png',
            subPops = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)],
            ylim = [0, 1.2],
            main = "!'Ancestry distribution of individuals at generation %d' % gen",
            legend = ['anc < 0.2', '0.2 <= anc < 0.4', '0.4 <= anc < 0.6',
                '0.6 <= anc < 0.8', '0.8 <= anc'],
            plot_axes = False,
            par_mar = [0, 0, 2, 0],
        ),
    ],
    duringOps = sim.pyTagger(passInfo, infoFields=['x', 'y', 'anc']),
    gen = 5,
)
#end_file

#begin_file log/histPlotter.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import histPlotter, qqPlotter, boxPlotter
import random
pop = sim.population([500], infoFields=['x', 'y', 'anc'])
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sim.sexSplitter())
def passInfo(fields):
    'Parental fields will be passed as x1, y1, anc1, x2, y2, anc2'
    x1, y1, anc1, x2, y2, anc2 = fields
    anc = (anc1 + anc2)/2.
    x = (x1 + x2)/2 + random.normalvariate(anc - 0.5, 0.1)
    y = (y1 + y2)/2 + random.normalvariate(0, 0.1)
    return x, y, anc

simu = sim.simulator(pop, sim.randomMating())
simu.evolve(
    initOps = [
        sim.initSex(),
        # random geographic location
        sim.initInfo(random.random, infoFields=['x', 'y']),
        # anc is 0 or 1
        sim.initInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    postOps = [
        histPlotter(infoFields='anc', 
            subPops=[(0,0), (0,1)], col_sp=['blue', 'red'],
            saveAs='log/histPlotter.png',
            main="!'Histogram of ancestry values at generation %d' % gen",
        ),
        qqPlotter(infoFields='anc', 
            subPops=[(0,0), (0,1)], col_sp=['blue', 'red'],
            saveAs='log/qqPlotter.png',
            main="!'QQ plot of ancestry values at generation %d' % gen",
        ),
        boxPlotter(infoFields='anc', 
            subPops=[(0,0), (0,1)],
            saveAs='whatever',
            dev_print_file='!"log/Gen%d.png" % gen',
            main="!'Boxplots of ancestry values at generation %d' % gen",
        ),
    ],
    duringOps = sim.pyTagger(passInfo, infoFields=['x', 'y', 'anc']),
    gen = 5,
)
#end_file

#begin_file log/getParam.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import types, simuOpt
options = [
    {'arg': 'r:',
     'longarg': 'rate=',
     'default': [0.01],
     'useDefault': True,
     'label': 'Recombination rate',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Recombination rate for each replicate. If a single value
            is given, it will be used for all replicates.''',
     'validate': simuOpt.valueListOf(simuOpt.valueBetween(0, 1))
    },
    {'longarg': 'rep=',
     'default': 5,
     'label': 'Number of replicates',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Number of replicates to simulate.',
     'validate': simuOpt.valueGT(0)
    }, 
    {'longarg': 'pop=',
     'default': 'CEU',
     'label': 'Initial sim.population',
     'allowedTypes': [types.StringType],
     'description': '''Use one of the HapMap sim.populations as the initial
            sim.population for this simulation. You can choose from:
            |YRI: 33 trios from the Yoruba people in Nigeria (Africa)
            |CEU: 30 trios from Utah with European ancestry (European)
            |CHB+JPT: 90 unrelated individuals from China and Japan (Asia)
            ''',
     'chooseOneOf': ['CEU', 'YRI', 'CHB+JPT'],
     'validate': simuOpt.valueOneOf(['CEU', 'YRI', 'CHB+JPT'])
    }
]
pars = simuOpt.simuParam(options, 'A demo simulation')
print pars.usage()
# You can manually feed parameters...
pars.processArgs(['--rep=10'])
pars.rep
#begin_ignore
import sys
oldArg = [x for x in sys.argv]
sys.argv.pop()
import os
if not os.path.isfile('figures/getParam.png'):
    print 'Run a GUI if getParam has not been runned'
else:
    sys.argv = ['getParam.py', '--rate=[0.25]', '--rep=5', '--pop="CEU"']
    simuOpt.setOptions(gui=False)

import simuPOP as sim
pars.processArgs(sys.argv)
#end_ignore
if not pars.getParam():
    sys.exit(1)

#begin_ignore
sys.argv = oldArg
#end_ignore
pars.saveConfig('sample.cfg')
# post-process parameters
pars.rate
pars.rep
pars.rate = pars.rate * pars.rep
# extract parameters as a dictionary or a list
pars.asDict()
pars.asList()
# Default value of parameter rep is changed
# additional attribute is added.
par1 = simuOpt.simuParam(options, # all parameters with default values
    rep=50,                     # default value of rep is changed
    additional=10               # derived parameters are added
)
# print all parameters except for derived ones.
print par1.asDict()
# All parameters are derived ...
par2 = simuOpt.simuParam(rep=50, pop='CEU', rate=[0.5])
print par2.asDict()
print par2.rep, par2.pop
#end_file

#begin_file log/paramFunc.py
import types, simuOpt
options = [
    simuOpt.param(
        arg = 'r:',
        longarg = 'rate=',
        default = [0.01],
        useDefault = True,
        label = 'Recombination rate',
        allowedTypes = [types.ListType, types.TupleType],
        description = '''Recombination rate for each replicate. If a single value
            is given, it will be used for all replicates.''',
        validate = simuOpt.valueListOf(simuOpt.valueBetween(0, 1))),
    simuOpt.param(
        longarg = 'rep=',
        default = 5,
        label = 'Number of replicates',
        allowedTypes = [types.IntType, types.LongType],
        description = 'Number of replicates to simulate.',
        validate = simuOpt.valueGT(0)), 
    simuOpt.param(
        longarg = 'pop=',
        default = 'CEU',
        label = 'Initial sim.population',
        allowedTypes = [types.StringType],
        description = '''Use one of the HapMap sim.populations as the initial
            sim.population for this simulation. You can choose from:
            |YRI: 33 trios from the Yoruba people in Nigeria (Africa)
            |CEU: 30 trios from Utah with European ancestry (European)
            |CHB+JPT: 90 unrelated individuals from China and Japan (Asia)
            ''',
        chooseOneOf = ['CEU', 'YRI', 'CHB+JPT'],
        validate = simuOpt.valueOneOf(['CEU', 'YRI', 'CHB+JPT']))
]
#begin_ignore
pars = simuOpt.simuParam(options, 'A demo simulation')
print pars.usage()
#end_ignore
#end_file


#begin_file log/reichDemo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import math
def demo_model(model, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    model: linear or exponential
    N0:   Initial sim.population size.
    N1:   Ending sim.population size.
    G0:   Length of burn-in stage.
    G1:   Length of sim.population expansion stage.
    '''
    def ins_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:
            return N1
    rate = (math.log(N1) - math.log(N0))/G1
    def exp_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    if model == 'instant':
        return ins_expansion
    elif model == 'exponential':
        return exp_expansion

# when needed, create a demographic function as follows
demo_func = demo_model('exponential', 1000, 100000, 500, 500)
# sim.population size at generation 700
print demo_func(700)
#end_file

#begin_file log/reichStat.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
class ne(sim.pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a sim.population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        sim.pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    #
    def calcNe(self, pop):
        sim.Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(freq[x]/sumFreq)**2 for x in freq.keys() if x != 0])
        # save the result to the sim.population.
        pop.dvars().ne = ne
        return True

def Ne(pop, loci):
    '''Function form of operator ne'''
    ne(loci).apply(pop)
    return pop.dvars().ne

pop = sim.population(100, loci=[10])
sim.InitByFreq(pop, [.2] * 5)
print Ne(pop, loci=[2, 4])
#end_file

#begin_file log/reichEvolve.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore

#begin_ignore
import math
def demo_model(model, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    model: linear or exponential
    N0:   Initial sim.population size.
    N1:   Ending sim.population size.
    G0:   Length of burn-in stage.
    G1:   Length of sim.population expansion stage.
    '''
    rate = (math.log(N1) - math.log(N0))/G1
    def ins_expansion(gen, pop=None):
        if gen < G0:
            return N0
        else:
            return N1
    
    def exp_expansion(gen, pop=None):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    
    if model == 'instant':
        return ins_expansion
    elif model == 'exponential':
        return exp_expansion

class ne(sim.pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a sim.population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        sim.pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    
    def calcNe(self, pop):
        sim.Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(freq[x]/sumFreq)**2 for x in freq.keys() if x != 0])
        
        # save the result to the sim.population.
        pop.dvars().ne = ne
        return True

#end_ignore

def simulate(model, N0, N1, G0, G1, spec, s, mu, k):
    '''Evolve a sim.population using given demographic model
    and observe the evolution of its allelic spectrum.
    model: type of demographic model.
    N0, N1, G0, G1: parameters of demographic model.
    spec: initial allelic spectrum, should be a list of allele
        frequencies for each allele.
    s: selection pressure.
    mu: mutation rate.
    k: k for the k-allele model
    '''
    demo_func = demo_model(model, N0, N1, G0, G1)
    simu = sim.simulator(
        sim.population(size=demo_func(0), loci=1, infoFields='fitness'),
        sim.randomMating(subPopSize=demo_func)
    )
    simu.evolve(
        initOps = [
            sim.initSex(),
            sim.initByFreq(loci=0, alleleFreq=spec)
        ],
        postOps = [
            sim.kamMutator(k=k, rates=mu),
            sim.maSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=[0], step=100),
            sim.pyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
                step=100),
        ],
        gen = G0 + G1
    )

simulate('instant', 1000, 10000, 500, 500, [0.9]+[0.02]*5, 0.01, 1e-4, 200)
#end_file

#begin_file log/simuCDCV.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
#!/usr/bin/env python
#
# Author:  Bo Peng
# Purpose: A real world example for simuPOP user's guide.
#
'''
Simulation the evolution of allelic spectra (number and frequencies
of alleles at a locus), under the influence of sim.population expansion,
mutation, and natural selection.
'''
import simuOpt
simuOpt.setOptions(quiet=True, alleleType='long')
import simuPOP as sim
import sys, types, os, math
options = [
    {'longarg': 'demo=',
     'default': 'instant',
     'label': 'Population growth model',
     'description': 'How does a sim.population grow from N0 to N1.',
     'chooseOneOf': ['instant', 'exponential'],
    },
    {'longarg': 'N0=',
     'default': 10000,
     'label': 'Initial sim.population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial sim.population size. This size will be maintained
                till the end of burnin stage''',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'N1=',
     'default': 100000,
     'label': 'Final sim.population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Ending sim.population size (after sim.population expansion)',
     'validate': simuOpt.valueGT(0)
    }, 
    {'longarg': 'G0=',
     'default': 500,
     'label': 'Length of burn-in stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of generations of the burn in stage.',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'G1=',
     'default': 1000,
     'label': 'Length of expansion stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of geneartions of the sim.population expansion stage',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'spec=',
     'default': [0.9] + [0.02]*5,
     'label': 'Initial allelic spectrum',
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''Initial allelic spectrum, should be a list of allele
            frequencies, for allele 0, 1, 2, ... respectively.''',
     'validate': simuOpt.valueListOf(simuOpt.valueBetween(0, 1)),
    },
    {'longarg': 's=',
     'default': 0.01,
     'label': 'Selection pressure',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': '''Selection coefficient for homozygtes (aa) genotype.
            A recessive selection model is used so the fitness values of
            genotypes AA, Aa and aa are 1, 1 and 1-s respectively.''',
     'validate': simuOpt.valueGT(-1),
    },
    {'longarg': 'mu=',
     'default': 1e-4,
     'label': 'Mutation rate',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': 'Mutation rate of a k-allele mutation model',
     'validate': simuOpt.valueBetween(0, 1),
    },
    {'longarg': 'k=',
     'default': 200,
     'label': 'sim.Maximum allelic sim.state',
     'allowedTypes': [types.IntType],
     'description': 'sim.Maximum allelic sim.state for a k-allele mutation model',
     'validate': simuOpt.valueGT(1),
    },
]

def demo_model(type, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    type: linear or exponential
    N0:   Initial sim.population size.
    N1:   Ending sim.population size.
    G0:   Length of burn-in stage.
    G1:   Length of sim.population expansion stage.
    '''
    rate = (math.log(N1) - math.log(N0))/G1
    def ins_expansion(gen, pop=None):
        if gen < G0:
            return N0
        else:
            return N1
    
    def exp_expansion(gen, pop=None):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    
    if type == 'instant':
        return ins_expansion
    elif type == 'exponential':
        return exp_expansion

class ne(sim.pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a sim.population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        sim.pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    
    def calcNe(self, pop):
        sim.Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc][1:]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(x/sumFreq)**2 for x in freq])
        # save the result to the sim.population.
        pop.dvars().ne = ne
        return True

def simuCDCV(model, N0, N1, G0, G1, spec, s, mu, k):
    '''Evolve a sim.population using given demographic model
    and observe the evolution of its allelic spectrum.
    model: type of demographic model.
    N0, N1, G0, G1: parameters of demographic model.
    spec: initial allelic spectrum, should be a list of allele
        frequencies for each allele.
    s: selection pressure.
    mu: mutation rate.
    k: maximum allele for the k-allele model
    '''
    demo_func = demo_model(model, N0, N1, G0, G1)
    print demo_func(0)
    simu = sim.simulator(
        sim.population(size=demo_func(0), loci=1, infoFields='fitness'),
        sim.randomMating(subPopSize=demo_func)
    )
    simu.evolve(
        initOps = [
            sim.initSex(),
            sim.initByFreq(loci=0, alleleFreq=spec)
        ],
        postOps = [
            sim.kamMutator(rate=mu, maxAllele=k),
            sim.maSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=0, step=100),
            sim.pyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
                step=100),
        ],
        gen = G0 + G1
    )
    return simu.extract(0)

if __name__ == '__main__':
    # get parameters
    par = simuOpt.simuParam(options, __doc__)
    if not par.getParam():
        sys.exit(1)
    
    if not sum(par.spec) == 1:
        print 'Initial allelic spectrum should add up to 1.'
        sys.exit(1)
    # save user input to a configuration file
    par.saveConfig('simuCDCV.cfg')
    #
    simuCDCV(*par.asList())

#begin_ignore
import os
if os.path.file('log/simuCDCV.py'):
    out = os.popen('python log/simuCDCV.py -h')
    hlp = open('log/simuCDCV.hlp', 'w')
    print >> hlp, out.read()
    hlp.close()
#end_ignore
#end_file



#begin_file log/ageStructured.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.GetRNG().setSeed(12345)
#end_ignore
import random
N = 10000
pop = sim.population(10000, loci=1, infoFields='age')
pop.setVirtualSplitter(sim.infoSplitter(field='age', cutoff=[20, 50, 75]))
def demoModel(gen, pop):
    '''A demographic model that keep a constant supply of new individuals'''
    # number of individuals that will die
    sim.Stat(pop, popSize=True, subPops=[(0,3)])
    # individuals that will be kept, plus some new guys.
    return pop.popSize() - pop.dvars().popSize + N / 75

def pene(geno, age, gen):
    'Define an age-dependent penetrance function'
    # the probability of getting disease increases with age
    return (0., 0.01*age[0], 0.01*age[0])[sum(geno)]

def outputStat(pop):
    'Calculate and output statistics'
    sim.Stat(pop, popSize=True, numOfAffected=True,
        subPops=[(0,0), (0,1), (0,2), (0,3)],
        vars=['popSize_sp', 'propOfAffected_sp'])
    for sp in range(3):
        print '%s: %.3f%% (size %d)' % (pop.subPopName((0,sp)),
            pop.dvars((0,sp)).propOfAffected * 100.,
            pop.dvars((0,sp)).popSize),
    print
    return True

simu = sim.simulator(pop,
    sim.heteroMating([
        # all individuals with age < 75 will be kept. Note that
        # cloneMating will keep individual sex, affection status and all
        # information fields (by default).
        sim.cloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
        # only individuals with age between 20 and 50 will mate and produce
        # offspring. The age of offspring will be zero.
        sim.randomMating(subPops=[(0,1)]),
    ], subPopSize=demoModel)
)

simu.evolve(
    initOps = [
        sim.initSex(),
        # random assign age
        sim.initInfo(lambda: random.randint(0, 75), infoFields='age'),
        # random genotype
        sim.initByFreq([0.5, 0.5]),
        sim.pyOutput('Prevalence of disease in each age group:\n'),
    ],
    # increase the age of everyone by 1 before mating.
    preOps = sim.infoExec('age += 1'),
    # number of individuals?
    postOps = [
        sim.pyPenetrance(func=pene, loci=0, paramFields='age'),
        sim.pyOperator(func=outputStat)
    ],
    gen = 10
)     
#end_file
