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
sim.getRNG().setSeed(12345)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=1000, loci=2)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([1, 2, 2, 1])
    ],
    matingScheme=sim.randomMating(ops=sim.Recombinator(rates=0.01)),
    postOps=[
        sim.stat(LD=[0, 1], step=10),
        sim.PyEval(r"'%.2f\n' % LD[0][1]", step=10),
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
sim.getRNG().setSeed(12345)
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
sim.getRNG().setSeed(12345)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=2, loci=[5, 6])
sim.initByFreq(pop, [0.2, 0.3, 0.5])
for ind in pop.individuals():
    for loc in range(pop.chromBegin(1), pop.chromEnd(1)):
        print ind.allele(loc),
    print 
#end_file

#begin_file log/defdict.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population([100]*2, loci=1)
sim.initByFreq(pop, [0, 0.2, 0.8], subPops=0)
sim.initByFreq(pop, [0.2, 0.8], subPops=1)
sim.Stat(pop, alleleFreq=0, vars=['alleleFreq_sp'])
for sp in range(2):
    print 'Subpop %d (with %d alleles): ' % (sp, len(pop.dvars(sp).alleleFreq[0])),
    for a in range(3):
        print '%.2f ' % pop.dvars(sp).alleleFreq[0][a],
    print

#end_file


#begin_file log/userFunc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(1000, loci=1, infoFields='smoking')
sim.initInfo(pop, lambda:random.randint(0,1), infoFields='smoking')
sim.initByFreq(pop, [0.3, 0.7])

# a penetrance function that depends on smoking
def func(geno, smoking):
    if smoking:
        return (geno[0]+geno[1])*0.4
    else:
        return (geno[0]+geno[1])*0.1

sim.pyPenetrance(pop, loci=0, func=func)
sim.Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected

#end_file

#begin_file log/withArgs.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(1000, loci=1, infoFields=('x', 'y'))
sim.initInfo(pop, lambda:random.randint(0,1), infoFields=('x', 'y'))
sim.initByFreq(pop, [0.3, 0.7])

# a penetrance function that depends on unknown information fields
def func(*fields):
    return 0.4*sum(fields)

# function withArgs tells PyPenetrance that func accepts fields x, y so that
# it will pass values at fields x and y to func.
sim.pyPenetrance(pop, loci=0, func=sim.withArgs(func, pop.infoFields()))
sim.Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#end_file

#begin_file log/genoStru.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2,5], ploidy=sim.HAPLODIPLOID, loci=[3, 5])
sim.initByFreq(pop, [0.3, 0.7])
sim.Dump(pop)
#end_file

#begin_file log/chromType.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=6, ploidy=2, loci=[3, 3, 6, 4, 4, 4],
    chromTypes=[sim.AUTOSOME]*2 + [sim.CHROMOSOME_X, sim.CHROMOSOME_Y] + [sim.CUSTOMIZED]*2)
sim.initByFreq(pop, [0.3, 0.7])
sim.Dump(pop, structure=False) # does not display genotypic structure information
#end_file

#begin_file log/infoField.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(10, loci=[20], ancGen=1,
    infoFields=['father_idx', 'mother_idx'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([0]*20+[1]*20)
    ],
    matingScheme=sim.randomMating(
        ops=[
            sim.Recombinator(rates=0.01),
            sim.ParentsTagger()
        ]
    ),
    gen = 1
)
pop.indInfo('mother_idx')  # mother of all offspring
ind = pop.individual(0)
mom = pop.ancestor(ind.mother_idx, 1)
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
sim.getRNG().setSeed(12345)
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
sim.getRNG().setSeed(12345)
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
#
geno.count(1)           # count
geno.index(2)           # index 
ind.setAllele(5, 3)    # change underlying genotype using setAllele
print geno              # geno is change
print geno             # but not geno
geno[2:5] = 4           # can use regular Python slice operation
print ind.genotype()
#end_file


#begin_file log/subPop.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
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
# remove subpopulations
pop.removeSubPops(1)
pop.subPopSizes()
#end_file

#begin_file log/subPopName.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
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
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[200, 400], loci=[30], infoFields='x')
# assign random information fields
sim.initSex(pop)
sim.initInfo(pop, lambda: random.randint(0, 3), infoFields='x')
# define a virtual splitter by sex
pop.setVirtualSplitter(sim.SexSplitter())
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 1])    # Size of VSP 0 in subpopulation 0
# define a virtual splitter by information field 'x'
pop.setVirtualSplitter(sim.InfoSplitter(field='x', values=[0, 1, 2, 3]))
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
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(10, loci=[2, 3], infoFields='Sex')
sim.initSex(pop)
pop.setVirtualSplitter(sim.SexSplitter())
# initialize male and females with different genotypes. 
sim.initByValue(pop, [[0]*5, [1]*5], subPops=([0, 0], [0, 1]))
# set Sex information field to 0 for all males, and 1 for all females
pop.setIndInfo([sim.MALE], 'Sex', [0, 0])
pop.setIndInfo([sim.FEMALE], 'Sex', [0, 1])
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
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[2000, 4000], loci=[30], infoFields='x')
# assign random information fields
sim.initSex(pop)
sim.initInfo(pop, lambda: random.randint(0, 3), infoFields='x')
#
# 1, use a combined splitter
pop.setVirtualSplitter(sim.CombinedSplitter(splitters = [
    sim.SexSplitter(),
    sim.InfoSplitter(field='x', values=[0, 1, 2, 3])
]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.MALE
pop.subPopSize([1, 4])    # individuals in sp 1 with value 2 at field x
#
# use a product splitter that defines additional VSPs by sex and info
pop.setVirtualSplitter(sim.ProductSplitter(splitters = [
    sim.SexSplitter(names=['M', 'F']),  # give a new set of names
    sim.InfoSplitter(field='x', values=[0, 1, 2, 3])
]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.MALE with value 1 in sp 0
pop.subPopSize([1, 5])    # sim.FEMALE with value 1 in sp 1
#
# use a combined splitter to join VSPs defined by a
# product splitter
pop.setVirtualSplitter(sim.CombinedSplitter([
    sim.ProductSplitter([
        sim.SexSplitter(),
        sim.InfoSplitter(field='x', values=[0, 1, 2, 3])])],
    vspMap = [[0,1,2], [4,5,6], [7]],
    names = ['Male x<=3', 'Female x<=3', 'Female x=4']))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # sim.MALE with value 0, 1, 2 at field x
pop.subPopSize([1, 1])    # sim.FEMALE with value 0, 1 or 2 at field x
#end_file


#begin_file log/accessIndividual.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
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
pop.setVirtualSplitter(sim.InfoSplitter(cutoff=[3, 7], field='x'))
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
sim.tagID(pop)
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
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[4, 6], loci=2, infoFields='x')
pop.setIndInfo([random.randint(0, 10) for x in range(10)], 'x')
pop.indInfo('x')
pop.setGenotype([0, 1, 2, 3], 0)
pop.genotype(0)
pop.setVirtualSplitter(sim.InfoSplitter(cutoff=[3], field='x'))
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
sim.getRNG().setSeed(12345)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(500, loci=1, ancGen=2)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme = sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0, begin=-3),
        sim.PyEval(r"'%.3f\n' % alleleFreq[0][0]", begin=-3)
    ],
    gen = 20
)
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
sim.getRNG().setSeed(12345)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(5, loci=[5], alleleNames=['A', 'T', 'C', 'G'])
sim.initByFreq(pop, [0.2, 0.3, 0.4, 0.1])
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
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[200, 200], loci=[5, 5], infoFields='age')
sim.initByValue(pop, range(10))
sim.initInfo(pop, lambda: random.randint(0,75), infoFields='age')
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[20, 60]))
# remove individuals
pop.removeIndividuals(indexes=range(0, 300, 10))
print pop.subPopSizes()
# remove subPopulation
pop.removeSubPops(1)
print pop.subPopSizes()
# remove virtual subpopulation (people with age between 20 and 60)
pop.removeSubPops([(0, 1)])
print pop.subPopSizes()
# extract another virtual subpopulation (people with age greater than 60)
pop1 = pop.extractSubPops([(0,2)])
sim.Dump(pop1, structure=False, max=10)
#end_file


#begin_file log/popVars.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
from pprint import pprint
pop = sim.population(100, loci=2)
sim.initByFreq(pop, [0.3, 0.7])
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
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=1), rep=5)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme = sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.TerminateIf('len(alleleFreq[0]) == 1')
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=5, chromNames=['chrom1'])
pop.dvars().name = 'my sim.population'
pop.save('sample.pop')
pop1 = sim.loadPopulation('sample.pop')
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=[20])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.8, 0.2])
    ],
    preOps=[
        sim.PyEval(r"'At the beginning of gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1])
    ],
    matingScheme = sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0, begin=80, step=10),
        sim.PyEval(r"'At the end of gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        sim.PyEval(r"'At the end of gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1])
    ],
    finalOps=sim.SavePopulation(output='sample.pop'),
    gen=100
)
#begin_ignore
import os
os.remove('sample.pop')
#end_ignore
#end_file

#begin_file log/replicate.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=[20]), 5)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.2, 0.8])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0, step=10),
        sim.PyEval('gen', step=10, reps=0),
        sim.PyEval(r"'\t%.2f' % alleleFreq[0][0]", step=10, reps=(0, 2, -1)),
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
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=1000, loci=2), rep=3)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([1, 2, 2, 1])
    ],
    matingScheme = sim.randomMating(ops=sim.Recombinator(rates=0.01)),
    postOps=[
        sim.stat(LD=[0, 1]),
        sim.PyEval(r"'%.2f\t' % LD[0][1]", step=20, output='>>LD.txt'),
        sim.pyOutput('\n', reps=-1, step=20, output='>>LD.txt'),
        sim.PyEval(r"'%.2f\t' % R2[0][1]", output='R2.txt'),
        sim.PyEval(r"'%.2f\t' % LD[0][1]", step=20, output="!'>>LD_%d.txt' % rep"),
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
sim.getRNG().setSeed(12345)
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
pop = sim.population(size=1000, loci=2)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([1, 2, 2, 1])
    ],
    matingScheme = sim.randomMating(ops=sim.Recombinator(rates=0.01)),
    postOps=[
        sim.stat(LD=[0, 1]),
        sim.PyEval(r"'LD: %d, %.2f' % (gen, LD[0][1])", step=20,
            output=logger.info),   # send LD to console and a logfile
        sim.PyEval(r"'R2: %d, %.2f' % (gen, R2[0][1])", step=20,
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=10000, loci=2)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([1, 2, 2, 1])
    ],
    matingScheme = sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(end=29),
        sim.Recombinator(rates=0.01, begin=30),
    ]),
    postOps=[
        sim.stat(LD=[0, 1]),
        sim.PyEval(r"'gen %d, LD: %.2f\n' % (gen, LD[0][1])", step=20)
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
sim.getRNG().setSeed(12345)
#end_ignore
def myPenetrance(geno):
    'A three-locus heterogeneity penetrance model'
    if sum(geno) < 2:
        return 0
    else:
        return sum(geno)*0.1

pop = sim.population(1000, loci=[20]*3)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.8, 0.2])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.PyPenetrance(func=myPenetrance, loci=[10, 30, 50]),
        sim.stat(numOfAffected=True),
        sim.PyEval(r"'%d: %d\n' % (gen, numOfAffected)")
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
sim.getRNG().setSeed(12345)
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
            sim.kamMutate(pop, k=2, rates=mu1, loci=[i])
        else:
            sim.kamMutate(pop, k=2, rates=mu2, loci=[i])
    return True

pop = sim.population(size=10000, loci=[2, 3])
pop.evolve(
    initOps=[ 
        sim.InitSex(),
        sim.InitByFreq([.99, .01], loci=[0, 2, 4]),
        sim.InitByFreq([.8, .2], loci=[1, 3])
    ],
    preOps=sim.pyOperator(func=dynaMutator, param=(.2, 1e-2, 1e-5)),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=range(5), step=10),
        sim.PyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
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
sim.getRNG().setSeed(12345)
#end_ignore
def rejectInd(off):
    'reject an individual if it off.allele(0) == 1'
    return off.allele(0) == 0

pop = sim.population(size=100, loci=1)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme=sim.randomMating(
        ops=[
            sim.MendelianGenoTransmitter(),
            sim.pyOperator(func=rejectInd)
        ]),
    gen = 1
)
# You should see no individual with allele 1 at locus 0, ploidy 0.
pop.genotype()[:20]
#end_file

#begin_file log/funcinitByFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
from simuPOP import InitByFreq, population
#begin_ignore
from simuPOP import getRNG
getRNG().setSeed(12345)
#end_ignore
def initByFreq(pop, *args, **kwargs):
    InitByFreq(*args, **kwargs).apply(pop)

pop = population(1000, loci=[2,3])
initByFreq(pop, [.2, .3, .5])
#end_file

#begin_file log/migrSize.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[500, 1000], infoFields='migrate_to')
pop.evolve(
    initOps=sim.InitSex(),
    preOps=sim.migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(popSize=True),
        sim.PyEval(r'"%s\n" % subPopSize')
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[500, 1000], infoFields='migrate_to')
pop.evolve(
    initOps=sim.InitSex(),
    preOps=[
        sim.migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        sim.stat(popSize=True),
        sim.PyEval(r'"%s\n" % subPopSize')
    ],
    matingScheme=sim.randomMating(subPopSize=[500, 1000]),
    postOps=[
        sim.stat(popSize=True),
        sim.PyEval(r'"%s\n" % subPopSize')
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
sim.getRNG().setSeed(12345)
#end_ignore
def demo(gen):
    return [500 + gen*10, 1000 + gen*10]

pop = sim.population(size=[500, 1000], infoFields='migrate_to')
pop.evolve(
    initOps=sim.InitSex(),
    preOps=sim.migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
    matingScheme=sim.randomMating(subPopSize=demo),
    postOps=[
        sim.stat(popSize=True),
        sim.PyEval(r'"%s\n" % subPopSize')
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
sim.getRNG().setSeed(12345)
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

pop = sim.population(1000)
pop.evolve(
    preOps=[
        sim.stat(popSize=True),
        sim.PyEval(r'"Gen %d:\t%s (before mating)\t" % (gen, subPopSize)')
    ],
    matingScheme=sim.RandomSelection(subPopSize=demo),
    postOps=[
        sim.stat(popSize=True),
        sim.PyEval(r'"%s (after mating)\n" % subPopSize')
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
sim.getRNG().setSeed(12345)
#end_ignore
def checkNumOffspring(numOffspring):
    '''Check the number of offspring for each family using
       information field father_idx
    '''
    pop = sim.population(size=[30], infoFields=['father_idx', 'mother_idx'])
    pop.evolve(
        initOps=sim.InitSex(),
        matingScheme=sim.randomMating(ops=[
            sim.MendelianGenoTransmitter(),
            sim.ParentsTagger(),
            ],
            numOffspring=numOffspring),
        gen=1)
    # get the parents of each offspring
    parents = [(x, y) for x, y in zip(pop.indInfo('mother_idx'),
        pop.indInfo('father_idx'))]
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
checkNumOffspring(numOffspring=2)
# Case 2: Use a Python function
import random
def func(gen):
    return random.randint(5, 8)

checkNumOffspring(numOffspring=func)
# Case 3: A geometric distribution
checkNumOffspring(numOffspring=(sim.GEOMETRIC_DISTRIBUTION, 0.3))
# Case 4: A Possition distribution
checkNumOffspring(numOffspring=(sim.POISSON_DISTRIBUTION, 3))
# Case 5: A Binomial distribution
checkNumOffspring(numOffspring=(sim.BINOMIAL_DISTRIBUTION, 0.1, 10))
# Case 6: A uniform distribution
checkNumOffspring(numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 6))
#end_file

#begin_file log/sexMode.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
def checkSexMode(ms):
    '''Check the assignment of sex to offspring'''
    pop = sim.population(size=[40])
    pop.evolve(initOps=sim.InitSex(), matingScheme=ms, gen=1)
    # return individual sex as a string
    return ''.join([ind.sexChar() for ind in pop.individuals()])

# Case 1: sim.NO_SEX (all male, sim.randomMating will not continue)
checkSexMode(sim.randomMating(sexMode=sim.NO_SEX))
# Case 2: sim.RANDOM_SEX (sim.Male/Female with probability 0.5)
checkSexMode(sim.randomMating(sexMode=sim.RANDOM_SEX))
# Case 3: sim.PROB_OF_MALES (Specify probability of male)
checkSexMode(sim.randomMating(sexMode=(sim.PROB_OF_MALES, 0.8)))
# Case 4: sim.NUM_OF_MALES (Specify number of male in each family)
checkSexMode(sim.randomMating(numOffspring=3, sexMode=(sim.NUM_OF_MALES, 1)))
# Case 5: sim.NUM_OF_FEMALES (Specify number of female in each family)
checkSexMode(sim.randomMating(
    numOffspring=(sim.UNIFORM_DISTRIBUTION, 4, 6),
    sexMode=(sim.NUM_OF_FEMALES, 2))
)
#end_file

#begin_file log/monogamous.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(20, infoFields=['father_idx', 'mother_idx'])
pop.evolve(
    initOps=sim.InitSex(sex=(sim.MALE, sim.FEMALE)),
    matingScheme=sim.MonogamousMating(
        numOffspring=2,
        sexMode=(sim.NUM_OF_MALES, 1),
        ops=[
            sim.MendelianGenoTransmitter(),
            sim.ParentsTagger(),
        ],
    ),
    gen = 5
)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, infoFields=['father_idx', 'mother_idx'])
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.PolygamousMating(polySex=sim.MALE, polyNum=2,
        ops=[sim.ParentsTagger(),
            sim.MendelianGenoTransmitter()],
    ),
    gen = 5
)
[int(ind.father_idx) for ind in pop.individuals()][:20]
[int(ind.mother_idx) for ind in pop.individuals()][:20]
#end_file

#begin_file log/RandomSelection.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, ploidy=1, loci=[5, 5], ancGen=1,
    infoFields='parent_idx')
pop.evolve(
    initOps=sim.InitByFreq([0.3, 0.7]),
    matingScheme=sim.RandomSelection(ops=[
        sim.ParentsTagger(infoFields='parent_idx'),
        sim.CloneGenoTransmitter(),
    ]),
    gen = 5
)
ind = pop.individual(0)
par = pop.ancestor(ind.parent_idx, 1)
print ind.sex(), ind.genotype()
print par.sex(), par.genotype()
#end_file

#begin_file log/AlphaMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=5, 
    infoFields=['father_idx', 'mother_idx', 'fitness'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    preOps=sim.MaSelector(loci=0, fitness=[0.8, 0.8, 1]),
    matingScheme=sim.AlphaMating(alphaSex=sim.MALE, alphaNum=2,
        ops=[sim.MendelianGenoTransmitter(), sim.ParentsTagger()]),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.PyEval(r'"%.2f\n" % alleleFreq[0][1]', step=5)
    ],
    gen = 20,
)
[int(ind.father_idx) for ind in pop.individuals()][:10]
[int(ind.mother_idx) for ind in pop.individuals()][:10]
#end_file

#begin_file log/HaplodiploidMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(10, ploidy=sim.HAPLODIPLOID, loci=[5, 5],
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(sim.SexSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([0]*10, subPops=[(0, 0)]),
        sim.InitByValue([1]*10+[2]*10, subPops=[(0, 1)])
    ],
    preOps=sim.dumper(structure=False),
    matingScheme=sim.HaplodiploidMating(
        ops=[sim.HaplodiploidGenoTransmitter(), sim.ParentsTagger()]),
    postOps=sim.dumper(structure=False),
    gen = 1
)
#end_file

#begin_file log/SelfMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(20, loci=8)
# every chromosomes are different. :-)
for idx, ind in enumerate(pop.individuals()):
    ind.setGenotype([idx*2], 0)
    ind.setGenotype([idx*2+1], 1)

pop.evolve(
    matingScheme=sim.SelfMating(ops=sim.Recombinator(rates=0.01)),
    gen = 1
)
sim.Dump(pop, width=3, structure=False, max=10)
#end_file

#begin_file log/heteroMatingSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000, 1000], loci=2,
    infoFields=['father_idx', 'mother_idx'])
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.heteroMating([
        sim.randomMating(numOffspring=2, subPops=0,
            ops=[sim.MendelianGenoTransmitter(), sim.ParentsTagger()]
        ),
        sim.randomMating(numOffspring=4, subPops=1,
            ops=[sim.MendelianGenoTransmitter(), sim.ParentsTagger()]
        )
    ]),
    gen=10
)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000], loci=2,
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(sim.ProportionSplitter([0.2, 0.8]))
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.heteroMating(matingSchemes=[
        sim.SelfMating(subPops=[(0, 0)],
            ops=[sim.SelfingGenoTransmitter(), sim.ParentsTagger()]
        ),
        sim.randomMating(subPops=[(0, 1)],
            ops=[sim.SelfingGenoTransmitter(), sim.ParentsTagger()]
        )
    ]),
    gen = 10
)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000], loci=2,
    infoFields='mark')
pop.setVirtualSplitter(sim.RangeSplitter([[0, 500], [200, 1000]]))
def markOff(param):
    '''define a Python during mating operator that marks
       individual information field 'mark'
    '''
    def func(off, param):
        off.mark = param
        return True
    return sim.pyOperator(func=func, param=param)

pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.heteroMating([
        sim.randomMating(subPops=0, weight=-0.5, ops=[markOff(0), sim.MendelianGenoTransmitter()]),
        sim.randomMating(subPops=[(0, 0)], weight=2, ops=[markOff(1), sim.MendelianGenoTransmitter()]),
        sim.randomMating(subPops=[(0, 1)], weight=3, ops=[markOff(2), sim.MendelianGenoTransmitter()])
    ]),
    gen = 10
)
marks = list(pop.indInfo('mark'))
marks.count(0.)
marks.count(1.)
marks.count(2.)
#end_file
#begin_file log/simulator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=10)
# five copies of the same population
simu = sim.simulator(pop, rep=5)
simu.numRep()
# evolve for ten generations and save the populations
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.3, 0.7])
    ],
    matingScheme=sim.randomMating(),
    finalOps=sim.SavePopulation('!"pop%d.pop"%rep'),
    gen=10
)
# load the population and create another simulator
simu = sim.simulator([sim.loadPopulation('pop%d.pop' % x) for x in range(5)])
# continue to evolve
simu.evolve(
    matingScheme=sim.randomMating(),
    gen=10
)
# print out allele frequency
for pop in simu.populations():
    sim.Stat(pop, alleleFreq=0)
    print '%.2f' % pop.dvars().alleleFreq[0][0],

print
# get a population
pop = simu.extract(0)
simu.numRep()
#begin_ignore
import os
for x in range(5):
    os.remove('pop%d.pop' % x)

#end_ignore
#end_file


#begin_file log/simuGen.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(50, loci=[10], ploidy=1),
    rep=3)
simu.evolve(gen = 5)
simu.dvars(0).gen
simu.evolve(
    initOps=[sim.InitByFreq([0.5, 0.5])],
    matingScheme=sim.RandomSelection(),
    postOps=[
        sim.stat(alleleFreq=5),
        sim.IfElse('alleleNum[5][0] == 0',
            sim.PyEval(r"'Allele 0 is lost in rep %d at gen %d\n' % (rep, gen)")),
        sim.IfElse('alleleNum[5][0] == 50',
            sim.PyEval(r"'Allele 0 is fixed in rep %d at gen %d\n' % (rep, gen)")),
        sim.TerminateIf('len(alleleNum[5]) == 1'),
    ],
)
[simu.dvars(x).gen for x in range(3)]
#end_file

#begin_file log/describe.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim

def outputStat(pop):
    'Calculate and output statistics, ignored'
    return True

# describe this evolutionary process
print sim.describe(
    initOps=[
        sim.InitSex(),
        sim.InitInfo(lambda: random.randint(0, 75), infoFields='age'),
        sim.InitByFreq([0.5, 0.5]),
        sim.IdTagger(),
        sim.pyOutput('Prevalence of disease in each age group:\n'),
    ],
    preOps=sim.InfoExec('age += 1'),
    matingScheme=sim.heteroMating([
        sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
        sim.randomMating(ops=[
            sim.IdTagger(),
            sim.Recombinator(intensity=1e-4)
        ], subPops=[(0,1)]),
    ]),
    postOps=[
        sim.MaPenetrance(loci=0, penetrance=[0.01, 0.1, 0.3]),
        sim.pyOperator(func=outputStat)
    ],
    gen = 100,
    numRep = 3
)     
#end_file

#begin_file log/twoStage.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
# First stage: use the standard random mating scheme, do not use any
# information field for efficiency considerations.
pop = sim.population(500, loci=[10])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme=sim.randomMating(),
    gen = 50
)
# Second stage: track parents and produce more offspring per mating
# event. In preparation for sim.pedigree ascertainment.
pop.addInfoFields(['ind_id', 'father_id', 'mother_id'])
pop.setAncestralDepth(1)
pop.evolve(
    initOps=sim.IdTagger(),
    matingScheme=sim.randomMating(numOffspring=2, ops=[
        sim.IdTagger(),
        sim.PedigreeTagger(),
        sim.MendelianGenoTransmitter(),
    ]),
    postOps=sim.MaPenetrance(loci=0, penetrance=(0.2, 0.4, 0.5)),
    gen = 5
)
# Sample affected sibpairs
from simuPOP.sampling import drawAffectedSibpairSample
sample = drawAffectedSibpairSample(pop, families=5)
[int(ind.father_id) for ind in sample.individuals()]
#end_file

#begin_file log/changeStru.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
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
simu = sim.simulator(sim.population(1000, loci=[]), rep=3)
simu.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.randomMating(),
    postOps=sim.pyOperator(func=mutator, param=(10000, 2e-6)),
    gen = 200
)
for pop in simu.populations():
    print pop.totNumLoci(), pop.lociPos()

#end_file


#begin_file log/locateRelative.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, ancGen=2, infoFields=['ind_id', 'father_id', 'mother_id'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
    ],
    matingScheme=sim.randomMating(
        numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
        ops=[
            sim.MendelianGenoTransmitter(),
            sim.IdTagger(),
            sim.PedigreeTagger()
        ],
    ),
    gen = 5
)
ped = sim.pedigree(pop)
offFields = ['off%d' % x for x in range(4)]
grandOffFields = ['grandOff%d' % x for x in range(5)]
ped.addInfoFields(['spouse'] + offFields + grandOffFields)
# only look spouse for fathers...
ped.locateRelatives(sim.OUTBRED_SPOUSE, ['spouse'], sex=sim.FEMALE_ONLY)
ped.locateRelatives(sim.COMMON_OFFSPRING, ['spouse'] + offFields)
# trace offspring of offspring
ped.traceRelatives([offFields, offFields], resultFields=grandOffFields)
# 
IDs = ped.individualsWithRelatives(grandOffFields)
# check on ID.
grandFather = IDs[0]
grandMother = ped.indByID(grandFather).spouse
# some ID might be invalid.
children = [ped.indByID(grandFather).info(x) for x in offFields]
childrenSpouse = [ped.indByID(x).spouse for x in children if x >= 1]
childrenParents = [ped.indByID(x).father_id for x in children if x >= 1] \
    + [ped.indByID(x).mother_id for x in children if x >= 1]
grandChildren = [ped.indByID(grandFather).info(x) for x in grandOffFields]
grandChildrenParents = [ped.indByID(x).father_id for x in grandChildren if x >= 1] \
    + [ped.indByID(x).mother_id for x in grandChildren if x >= 1]

def idString(IDs):
    uniqueIDs = list(set(IDs))
    uniqueIDs.sort()
    return ', '.join(['%d' % x for x in uniqueIDs if x >= 1])

print '''GrandParents: %d, %d
Children: %s
Children's spouses: %s
Children's parents: %s
GrandChildren: %s
GrandChildren's parents: %s ''' % \
(grandFather, grandMother, idString(children), idString(childrenSpouse),
    idString(childrenParents), idString(grandChildren), idString(grandChildrenParents))

#end_file


#begin_file log/InitSex.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000, 1000])
sim.initSex(pop, maleFreq=0.3, subPops=0)
sim.initSex(pop, sex=[sim.MALE, sim.FEMALE, sim.FEMALE], subPops=1)
sim.Stat(pop, numOfMales=True, vars='numOfMales_sp')
print pop.dvars(0).numOfMales
print pop.dvars(1).numOfMales
#end_file

#begin_file log/InitByFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2, 3], loci=[5, 7])
sim.initByFreq(pop, alleleFreq=[[.2, .8], [.8, .2]])
sim.Dump(pop, structure=False)
#end_file

#begin_file log/InitByFreqIdenticalInds.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2, 3], loci=[5, 7])
sim.initByFreq(pop, alleleFreq=[.2, .8], identicalInds=True)
sim.Dump(pop, structure=False)
#end_file

#begin_file log/InitByValue.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2, 3], loci=[5, 7])
sim.initByValue(pop, [1]*5 + [2]*7 + [3]*5 +[4]*7)
sim.Dump(pop, structure=False)
#end_file

#begin_file log/InitByValueProp.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[6, 8], loci=[5, 7])
pop.setVirtualSplitter(sim.SexSplitter())
# initialize sex and the first two loci
sim.initSex(pop)
sim.initByValue(pop, loci=range(5), value=range(10))
# initialize all males
sim.initByValue(pop, loci=range(5, 12), value=[2]*7,
    subPops=[(0, 0), (1, 0)])
# initialize females by proportion
sim.initByValue(pop, loci=range(5, 12), ploidy=1, value=[[3]*7, [4]*7],
    subPops=[(0, 1), (1, 1)], proportions=[0.4, 0.6])
sim.Dump(pop, structure=False)
#end_file

#begin_file log/InitInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=[5], loci=[2], infoFields=['sex', 'age'])
pop.setVirtualSplitter(sim.SexSplitter())
sim.initSex(pop)
sim.initInfo(pop, 0, subPops=[(0,0)], infoFields='sex')
sim.initInfo(pop, 1, subPops=[(0,1)], infoFields='sex')
sim.initInfo(pop, lambda: random.randint(20, 70), infoFields='age')
sim.Dump(pop, structure=False)
#end_file

#begin_file log/dumper.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[10, 10], loci=[20, 30], infoFields='gen',
    ancGen=-1)
pop.setVirtualSplitter(sim.SexSplitter())
pop1 = pop.clone()
sim.initByFreq(pop, [0]*20 + [0.1]*10)
pop.setIndInfo(1, 'gen')
sim.initByFreq(pop1, [0]*50 + [0.1]*10)
pop1.setIndInfo(2, 'gen')
pop.push(pop1)
sim.Dump(pop, width=3, loci=[5, 6, 30], subPops=([0, 0], [1, 1]),
    max=10, structure=False, ancGen=-1)
#end_file

#begin_file log/SavePopulation.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=2),
    rep=5)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.2, 0.8])
    ],
    matingScheme=sim.randomMating(),
    postOps=sim.SavePopulation(output="!'snapshot_%d_%d.pop' % (rep, gen)",
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


#begin_file log/IfElseFixed.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=1000, loci=1)
verbose = True
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5]),
    ],
    matingScheme=sim.randomMating(),
    postOps=sim.IfElse(verbose,
        ifOps=[
            sim.stat(alleleFreq=0),
            sim.PyEval(r"'Gen: %3d, allele freq: %.3f\n' % (gen, alleleFreq[0][1])",
                step=5)
        ],
        begin=10),
    gen = 30
)
#end_file


#begin_file log/IfElse.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=1000, loci=1),
    rep=4)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5]),
        sim.PyExec('below40, above60 = 0, 0')
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.IfElse('alleleFreq[0][1] < 0.4',
            sim.PyExec('below40 += 1')),
        sim.IfElse('alleleFreq[0][1] > 0.6',
            sim.PyExec('above60 += 1')),
        sim.IfElse('len(alleleFreq[0]) == 1',
            sim.PyExec('stoppedAt = gen')),
        sim.TerminateIf('len(alleleFreq[0]) == 1')
    ]
)
for pop in simu.populations():
    print 'Overall: %4d, below 40%%: %4d, above 60%%: %4d' % \
        (pop.dvars().stoppedAt, pop.dvars().below40, pop.dvars().above60)

#end_file

#begin_file log/TerminateIf.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(
    sim.population(size=100, loci=1),
    rep=10)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5]),
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.TerminateIf('len(alleleFreq[0]) == 1', stopAll=True)
    ]
)

#end_file


#begin_file log/Pause.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100), rep=10)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme=sim.randomMating(),
    postOps=[sim.Pause(stopOnKeyStroke=str(x), reps=x) for x in range(10)],
    gen = 100
)
#end_file

#begin_file log/TicToc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(10000, loci=[100]*5), rep=2)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.1, 0.9])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.TicToc(step=50, reps=-1),
    ],
    gen = 101
)
#end_file

#begin_file log/PyExec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(100, loci=1),
    rep=2)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.2, 0.8]),
        sim.PyExec('traj=[]')
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.PyExec('traj.append(alleleFreq[0][1])'),
    ],
    gen=5
)
# print trajectory
print ', '.join(['%.3f' % x for x in simu.dvars(0).traj])
#end_file

#begin_file log/PyEval.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=1,
    infoFields=['mother_idx', 'father_idx'])
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.ParentsTagger(),
    ]),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.PyEval(r'"gen %d, #father %d, #mother %d\n"' \
            ' % (gen, numFather, numMother)',
            stmts="numFather = len(set(pop.indInfo('father_idx')))\n"
                "numMother = len(set(pop.indInfo('mother_idx')))",
            exposePop='pop')
    ],
    gen=3
)
#end_file

#begin_file log/InfoEval.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(20, loci=1, infoFields='a')
pop.setVirtualSplitter(sim.InfoSplitter('a', cutoff=[3]))
sim.initByFreq(pop, [0.2, 0.8])
pop.setIndInfo([random.randint(2, 5) for x in range(20)], 'a')
sim.infoEval(pop, 'a', subPops=[(0, 0)]);print
sim.infoEval(pop, 'ind.allele(0, 0)', exposeInd='ind');print
# use sim.population variables
pop.dvars().b = 5
sim.infoEval(pop, '"%d " % (a+b)', usePopVars=True);print
#end_file

#begin_file log/InfoExec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=1, infoFields=['a', 'b', 'c'])
sim.initByFreq(pop, [0.2, 0.8])
sim.infoExec(pop, 'a=1')
print pop.indInfo('a')[:10]
sim.infoExec(pop, 'b=ind.sex()', exposeInd='ind')
print pop.indInfo('b')[:10]
sim.infoExec(pop, 'c=a+b')
print pop.indInfo('c')[:10]
pop.dvars().d = 5
sim.infoExec(pop, 'c+=d', usePopVars=True)
print pop.indInfo('c')[:10]
#end_file

#begin_file log/migrateByProb.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000]*3, infoFields='migrate_to')
pop.evolve(
    initOps=sim.InitSex(),
    preOps=sim.migrator(rate=[
            [0, 0.1, 0.1],
            [0, 0, 0.1],
            [0, 0.1, 0]
        ]), 
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(popSize=True),
        sim.PyEval('subPopSize'),
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000]*3, infoFields='migrate_to')
pop.evolve(
    initOps=sim.InitSex(),
    preOps=sim.migrator(rate=[[0.1], [0.2]],
            mode=sim.BY_PROPORTION,
            subPops=[1, 2],
            toSubPops=[3]),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(popSize=True),
        sim.PyEval('subPopSize'),
        sim.pyOutput('\n')
    ],
    gen = 5
)        
#
pop.evolve(
    preOps=sim.migrator(rate=[[50, 50], [100, 50]],
            mode=sim.BY_COUNTS,
            subPops=[3, 2],
            toSubPops=[2, 1]),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(popSize=True),
        sim.PyEval('subPopSize'),
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000]*2, infoFields='migrate_to')
pop.setVirtualSplitter(sim.SexSplitter())
pop.evolve(
    # 500 males and 500 females
    initOps=sim.InitSex(sex=[sim.MALE, sim.FEMALE]),
    preOps=[
        sim.migrator(rate=[
            [0, 0.10],
            [0, 0.05],
            ],
            mode = sim.BY_PROPORTION,
            subPops=[(0, 0), (0, 1)]),
        sim.stat(popSize=True, numOfMales=True, vars='numOfMales_sp'),
        sim.PyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMales'], subPopSize[0], "
            "subPop[1]['numOfMales'], subPopSize[1])"),
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(popSize=True, numOfMales=True, vars='numOfMales_sp'),
        sim.PyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMales'], subPopSize[0], "
            "subPop[1]['numOfMales'], subPopSize[1])"),
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population([10]*2, infoFields='migrate_to')
pop.setIndInfo([0, 1, 2, 3]*5, 'migrate_to')
sim.Migrate(pop, mode=sim.BY_IND_INFO)
pop.subPopSizes()
#end_file

#begin_file log/splitBySize.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000)
pop.evolve(
    preOps=[
        sim.SplitSubPops(subPops=0, sizes=[300, 300, 400], at=2),
        sim.stat(popSize=True),
        sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    matingScheme=sim.RandomSelection(),
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
sim.getRNG().setSeed(12345)
#end_ignore
def demo(gen, pop):
    if gen < 2:
        return 1000 + 100 * gen
    else:
        return [x + 50 * gen for x in pop.subPopSizes()]

pop = sim.population(1000)
pop.evolve(
    preOps=[
        sim.SplitSubPops(subPops=0, proportions=[.5]*2, at=2),
        sim.stat(popSize=True),
        sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    matingScheme=sim.RandomSelection(subPopSize=demo),
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
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population([1000]*3, subPopNames=['a', 'b', 'c'], infoFields='x')
pop.setIndInfo([random.randint(0, 3) for x in range(1000)], 'x')
print pop.subPopSizes()
print pop.subPopNames()
sim.splitSubPops(pop, subPops=[0, 2], infoFields=['x'])
print pop.subPopSizes()
print pop.subPopNames()
#end_file

#begin_file log/MergeSubPops.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population([500]*2)
pop.evolve(
    preOps=[
        sim.MergeSubPops(subPops=[0, 1], at=3),
        sim.stat(popSize=True),
        sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    matingScheme=sim.RandomSelection(),
    gen = 5
)
#end_file

#begin_file log/ResizeSubPops.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population([500]*2)
pop.evolve(
    preOps=[
        sim.ResizeSubPops(proportions=(1.5, 2), at=3),
        sim.stat(popSize=True),
        sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    matingScheme=sim.RandomSelection(),
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
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=[1000], loci=[100]),
    rep=2)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([0]*100 + [1]*100)
    ],
    matingScheme=sim.randomMating(ops = [
        sim.Recombinator(rates=0.01, reps=0),
        sim.Recombinator(rates=[0.01]*10, loci=range(50, 60), reps=1),
    ]),
    postOps=[
        sim.stat(LD=[[40, 55], [60, 70]]),
        sim.PyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[60][70])'),
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000], loci=3, lociPos=[0, 1, 1.1])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([0]*3 + [1]*3)
    ],
    matingScheme=sim.randomMating(ops=sim.Recombinator(intensity=0.01)),
    postOps=[
        sim.stat(LD=[[0, 1], [1, 2]]),
        sim.PyEval(r'"%.3f\t%.3f\n" % (LD_prime[0][1], LD_prime[1][2])', step=10)
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
sim.getRNG().setSeed(12345)
#end_ignore
simu = sim.simulator(sim.population(size=[1000], loci=[100]),
    rep=2)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([0]*100 + [1]*100)
    ],
    matingScheme=sim.randomMating(ops=[
        sim.Recombinator(rates=0.01, loci=50, reps=0),
        sim.Recombinator(rates=0.01, loci=50, reps=1, convMode=(sim.NUM_MARKERS, 1, 10)),
    ]),
    postOps=[
        sim.stat(LD=[[40, 55], [40, 70]]),
        sim.PyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[40][70])'),
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=[1000, 2000], infoFields='ind_id')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
    ],
    matingScheme=sim.randomMating(ops = [
        sim.IdTagger(),
        sim.Recombinator(rates=0.001, output='>>rec.log', infoFields='ind_id')]),
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

#begin_file log/MatrixMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=1)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.2, 0.3, 0.5])
    ],
    preOps=sim.MatrixMutator(rate = [
            [0, 1e-5, 1e-5],
            [1e-4, 0, 1e-4],
            [1e-3, 1e-3, 0]
        ]),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0, step=100),
        sim.PyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(3)]) + '\n'",
            step=100),
    ],
    gen=1000
)
#end_file

#begin_file log/KamMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=1*3)
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.KamMutator(k=5, rates=[1e-2, 1e-3], loci=[0, 1]),
        sim.stat(alleleFreq=range(3), step=100),
        sim.PyEval(r"', '.join(['%.3f' % alleleFreq[x][0] for x in range(3)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/SnpMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=[1, 1], infoFields='fitness')
pop.evolve(
    initOps=sim.InitSex(),
    preOps=[
        sim.SnpMutator(u=0.001),
        sim.MaSelector(loci=0, fitness=[1, 0.99, 0.98]),
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=[0, 1], step=100),
        sim.PyEval(r"'%.3f\t%.3f\n' % (alleleFreq[0][1], alleleFreq[1][1])",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/AcgtMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=1,
    alleleNames=['A', 'C', 'G', 'T'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([.1, .1, .1, .7])
    ],
    matingScheme=sim.randomMating(),
    preOps=[
        sim.AcgtMutator(rate=[1e-4, 0.5], model='K80'),
        sim.stat(alleleFreq=0, step=100),
        sim.PyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(4)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/SmmMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=1000, loci=[1, 1])
pop.evolve(
    # all start from allele 50
    initOps=[
        sim.InitSex(),
        sim.InitByFreq( [0]*50 + [1])
    ],
    matingScheme=sim.randomMating(),
    preOps=[
        sim.SmmMutator(rates=1e-3, loci=0),
        sim.SmmMutator(rates=1e-3, incProb=0.6, loci=1,
            mutStep=(sim.GEOMETRIC_DISTRIBUTION, 0.2)),
    ],
    gen=100
)
# count the average number tandem repeats at both loci
cnt0 = cnt1 = 0
for ind in pop.individuals():
    cnt0 += ind.allele(0, 0) + ind.allele(0, 1)
    cnt1 += ind.allele(1, 0) + ind.allele(1, 1)

print 'Average number of repeats at two loci are %.2f and %.2f.' % \
    (cnt0/2000., cnt1/2000.)
#end_file

#begin_file log/PyMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
def incAllele(allele):
    return allele + random.randint(1, 5)

pop = sim.population(size=1000, loci=[20])
pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.randomMating(),
    postOps=sim.PyMutator(func=incAllele, rates=[1e-4, 1e-3],
            loci=[2, 10]),
    gen = 1000
)
# count the average number tandem repeats at both loci
def avgAllele(pop, loc):
    ret = 0
    for ind in pop.individuals():
        ret += ind.allele(loc, 0) + ind.allele(loc, 1)
    return ret / (pop.popSize() * 2.)

print 'Average number of repeats at two loci are %.2f and %.2f.' % \
    (avgAllele(pop, 2), avgAllele(pop, 10))
#end_file

#begin_file log/MixedMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(5000, loci=[1, 1])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([50, 50])
    ],
    preOps=[
        # the first locus uses a pure stepwise mutation model
        sim.SmmMutator(rates=0.001, loci=0),
        # the second locus uses a mixed model
        sim.MixedMutator(rates=0.001, loci=1, mutators=[        
            sim.KamMutator(rates=1, k=100),
            sim.SmmMutator(rates=1)
        ], prob=[0.1, 0.9])],
    matingScheme=sim.randomMating(),
    gen = 20
)
# what alleles are there?
geno0 = []
geno1 = []
for ind in pop.individuals():
    geno0.extend([ind.allele(0, 0), ind.allele(0, 1)])
    geno1.extend([ind.allele(1, 0), ind.allele(1, 1)])

print 'Locus 0 has alleles', ', '.join([str(x) for x in set(geno0)])
print 'Locus 1 has alleles', ', '.join([str(x) for x in set(geno1)])
#end_file

#begin_file log/ContextMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(5000, loci=[3, 3])
pop.evolve(
    # initialize locus by 0, 0, 0, 1, 0, 1
    initOps=[
        sim.InitSex(),
        sim.InitByValue([1, 1], loci=[3, 5])
    ],
    preOps=[
        sim.ContextMutator(mutators=[
            sim.SnpMutator(u=0.1),
            sim.SnpMutator(u=1),
            ],
            contexts=[(0, 0), (1, 1)],
            loci=[1, 4],
            rates=0.01
        ),
        sim.stat(alleleFreq=[1, 4], step=5),
        sim.PyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    matingScheme=sim.randomMating(),
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
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(5000, loci=[3, 3])
def contextMut(allele, context):
    if context == [0, 0]:
        if allele == 0 and random.random() < 0.1:
            return 1
    elif context == [1, 1]:
        if allele == 0:
            return 1
    # do not mutate
    return allele

pop.evolve(
    # initialize locus by 0, 0, 0, 1, 0, 1
    initOps=[
        sim.InitSex(),
        sim.InitByValue([1, 1], loci=[3, 5])
    ],
    preOps=[
        sim.PyMutator(func=contextMut, context=1,
            loci=[1, 4],  rates=0.01
        ),
        #sim.SnpMutator(u=0.01, v= 0.01, loci=[1, 4]),
        sim.stat(alleleFreq=[1, 4], step=5),
        sim.PyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    matingScheme=sim.randomMating(),
    gen = 20
)
#end_file

#begin_file log/PointMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=1, infoFields='fitness')
pop.evolve(
    initOps=sim.pyOutput('Introducing alleles at generation'),
    preOps=sim.MaSelector(loci=0, wildtype=0, fitness=[1, 1.05, 1.1]),
    matingScheme=sim.RandomSelection(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.IfElse('alleleNum[0][1] == 0', ifOps=[
            sim.PyEval(r"' %d' % gen"),
            sim.PointMutator(inds=0, loci=0, allele=1),
        ]),
        sim.IfElse('alleleFreq[0][1] > 0.05', ifOps=[
            sim.PyEval(r"'.\nTerminate at generation %d at allele freq %.3f.\n'" +
                " % (gen, alleleFreq[0][1])"),
            sim.TerminateIf('True'),
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
sim.getRNG().setSeed(12345)
#end_ignore
def fragileX(geno):
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
pop.setVirtualSplitter(sim.AffectionSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([50, 50])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        # determine affection sim.status for each offspring (duringMating)
        sim.PyPenetrance(func=fragileX, loci=0),
        # unaffected offspring, mutation rate is high to save some time
        sim.SmmMutator(rates=1e-3, loci=1),
        # unaffected offspring, mutation rate is high to save some time
        sim.SmmMutator(rates=1e-3, loci=0, subPops=[(0, 0)]),
        # affected offspring have high probability of mutating upward
        sim.SmmMutator(rates=1e-2, loci=0, subPops=[(0, 1)],
           incProb=0.7, mutStep=3),
        # number of affected
        sim.pyOperator(func=avgAllele, step=20),
        sim.PyEval(r"'Gen: %3d #Aff: %d AvgRepeat: %.2f (unaff), %.2f (aff), %.2f (unrelated)\n'"
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[2000], loci=1)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0]*4 + [0.1, 0.2, 0.3, 0.4])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.KamMutator(k=4, rates=1e-4, mapIn=[0]*4 + range(4),
            mapOut=[4, 5, 6, 7]),
        sim.stat(alleleFreq=0, step=100),
        sim.PyEval(r"', '.join(['%.2f' % alleleFreq[0][x] for x in range(8)]) + '\n'",
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
sim.getRNG().setSeed(12345)
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
            loc = sim.getRNG().randGeometric(rate)
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
pop.evolve(
    initOps=sim.InitSex(),
    preOps=[
        # recombine in a 10Mb region at rate 1e-8
        sim.pyOperator(func=infSitesMutate, param=(1, 10000000, 1e-8)),
    ],
    matingScheme=sim.randomMating(),
    gen = 100
)
# now, we get a sim.population. Let us have a look at the 'alleles'.
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population([5000]*3, loci=5)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(structure=range(5), subPops=(0, 1), suffix='_01', step=40),
        sim.stat(structure=range(5), subPops=(1, 2), suffix='_12', step=40),
        sim.stat(structure=range(5), subPops=(0, 2), suffix='_02', step=40),
        sim.stat(structure=range(5), step=40),
        sim.PyEval(r"'Fst=%.3f (pairwise: %.3f %.3f %.3f)\n' % (F_st, F_st_01, F_st_12, F_st_02)",
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(10000, loci=1)
pop.setVirtualSplitter(sim.CombinedSplitter(
    [sim.SexSplitter(), sim.AffectionSplitter()]))
sim.initSex(pop)
sim.initByFreq(pop, [0.2, 0.8])
sim.maPenetrance(pop, loci=0, penetrance=[0.1, 0.2, 0.5])
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
sim.Stat(pop, numOfMales=True)
print pop.dvars().numOfMales
# Count the number of affected and unaffected male individual
sim.Stat(pop, numOfMales=True, subPops=[(0, 2), (0, 3)], vars='numOfMales_sp')
print pop.dvars((0,2)).numOfMales, pop.dvars((0,3)).numOfMales
# or number of affected male and females
sim.Stat(pop, numOfAffected=True, subPops=[(0, 0), (0, 1)], vars='numOfAffected_sp')
print pop.dvars((0,0)).numOfAffected, pop.dvars((0,1)).numOfAffected
# These can also be done using a sim.ProductSplitter...
pop.setVirtualSplitter(sim.ProductSplitter(
    [sim.SexSplitter(), sim.AffectionSplitter()]))
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(10000, loci=1)
pop.setVirtualSplitter(sim.AffectionSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(loci=0, alleleFreq=[0.8, 0.2])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.MaPenetrance(penetrance=[0.1, 0.4, 0.6], loci=0),
        sim.stat(alleleFreq=0, subPops=[(0, 0), (0, 1)],
            vars=['alleleFreq', 'alleleFreq_sp']),
        sim.PyEval(r"'Gen: %d, freq: %.2f, freq (aff): %.2f, freq (unaff): %.2f\n' % " + \
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=[1, 1, 1], chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y])
sim.initByFreq(pop, [0.01, 0.05, 0.94])
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=1)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(heteroFreq=0, step=10),
        sim.PyEval(r"'Gen: %d, HeteroFreq: %.2f\n' % (gen, heteroFreq[0])", step=20)
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
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import viewVars
pop = sim.population(100, loci=3)
sim.initByFreq(pop, [0.2, 0.4, 0.4], loci=0)
sim.initByFreq(pop, [0.2, 0.8], loci=2)
sim.Stat(pop, genoFreq=[0, 1, 2], haploFreq=[0, 1, 2],
    vars=['genoNum', 'haploFreq'])
viewVars(pop.vars(), gui=False)
#end_file

#begin_file log/statInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population([500], infoFields='anc')
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sim.InfoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
#
pop.evolve(
    initOps=[
        sim.InitSex(),
        # anc is 0 or 1
        sim.InitInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.InheritTagger(mode=sim.MEAN, infoFields='anc')
    ]),
    postOps=[
        sim.stat(popSize=True, meanOfInfo='anc', varOfInfo='anc',
            subPops=[(0,x) for x in range(5)]),
        sim.PyEval(r"'Anc: %.2f (%.2f), #inds: %s\n' %" + \
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population([1000]*2, loci=3)
sim.initByFreq(pop, [0.2, 0.8], subPops=0)
sim.initByFreq(pop, [0.8, 0.2], subPops=1)
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
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import *
from simuPOP.sampling import drawCaseControlSample
def assoTest(pop):
    'Draw case-control sample and apply association tests'
    sample = drawCaseControlSample(pop, cases=500, controls=500)
    sim.Stat(sample, association=(0, 2), vars=['Allele_ChiSq_p', 'Geno_ChiSq_p', 'Armitage_p'])
    print 'Allele test: %.2e, %.2e, Geno test: %.2e, %.2e, Trend test: %.2e, %.2e' \
        % (sample.dvars().Allele_ChiSq_p[0], sample.dvars().Allele_ChiSq_p[2],
        sample.dvars().Geno_ChiSq_p[0], sample.dvars().Geno_ChiSq_p[2],
        sample.dvars().Armitage_p[0], sample.dvars().Armitage_p[2])
    return True

pop = sim.population(size=100000, loci=3)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([[0]*3, [1]*3], proportions=[0.5, 0.5])
    ],
    matingScheme=sim.randomMating(ops=sim.Recombinator(loci=[0, 1], rates=[0.01, 0.005])),
    postOps=[
        sim.MaPenetrance(loci=1, penetrance=[0.1, 0.2, 0.4]),
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
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import migrIslandRates
simu = sim.simulator(sim.population([5000]*3, loci=10, infoFields='migrate_to'),
    rep=2)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    preOps=sim.migrator(rate=migrIslandRates(0.01, 3), reps=1),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(structure=range(10), step=40),
        sim.PyEval("'Fst=%.3f (rep=%d without migration) ' % (F_st, rep)", step=40, reps=0),
        sim.PyEval("'Fst=%.3f (rep=%d with migration) ' % (F_st, rep)", step=40, reps=1),
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population([1000], loci=1)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([[0,0], [0, 1], [1,1]], proportions=[0.4, 0.4, 0.2])
    ],
    preOps=[
        sim.stat(HWE=0, genoFreq=0),
        sim.PyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
            'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])'),
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(HWE=0, genoFreq=0),
        sim.PyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
            'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])'),
    ],
    gen = 1
)
#end_file

#begin_file log/InheritTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[1000]*10, loci=1, infoFields='x')
# tag the first individual of each subpopulation.
for sp in range(pop.numSubPop()):
    pop.individual(0, sp).x = 1

pop.evolve(
    initOps=sim.InitSex(),
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.InheritTagger(mode=sim.MAXIMUM, infoFields='x'),
    ]),
    postOps=[
        sim.stat(sumOfInfo='x', vars=['sumOfInfo_sp']),
        sim.PyEval(r'", ".join(["%3d" % subPop[i]["sumOfInfo"]["x"] for i in range(10)])+"\n"'),
    ],
    gen = 5
)
#end_file

#begin_file log/SummaryTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
sim.IdTagger().reset(1)
#end_ignore
pop = sim.population(1000, loci=1, infoFields=['fitness', 'avgFitness'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5]),
    ],
    preOps=sim.MaSelector(loci=0, wildtype=0, fitness=[1, 0.99, 0.95]),
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.SummaryTagger(mode=sim.MEAN, infoFields=['fitness', 'avgFitness']),
    ]),
    postOps=[
        sim.stat(alleleFreq=0, meanOfInfo='avgFitness', step=10),
        sim.PyEval(r"'gen %d: allele freq: %.3f, average fitness of parents: %.3f\n' % "
            "(gen, alleleFreq[0][1], meanOfInfo['avgFitness'])", step=10)
    ],
    gen = 50,
)
#end_file


#begin_file log/IdTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
sim.IdTagger().reset(1)
#end_ignore
pop = sim.population(10, infoFields='ind_id', ancGen=1)
pop.evolve(
    initOps=sim.IdTagger(),
    matingScheme=sim.RandomSelection(ops=[
        sim.CloneGenoTransmitter(),
        sim.IdTagger(),
    ]),
    gen = 1
)
print [int(ind.ind_id) for ind in pop.individuals()]
pop.useAncestralGen(1)
print [int(ind.ind_id) for ind in pop.individuals()]
sim.tagID(pop) # re-assign ID
print [int(ind.ind_id) for ind in pop.individuals()]
#end_file

#begin_file log/PedigreeTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
sim.IdTagger().reset(1)
#end_ignore
pop = sim.population(100, infoFields=['ind_id', 'father_id', 'mother_id'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger()
    ],
    matingScheme=sim.randomMating(ops=[
        sim.IdTagger(),
        sim.PedigreeTagger(output=">>sim.pedigree.txt"),
        sim.MendelianGenoTransmitter()]
    ),
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

#begin_file log/PyTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
def randomMove(x, y):
    '''Pass parental information fields to offspring'''
    # shift right with high concentration of alleles... 
    off_x = random.normalvariate((x[0]+x[1])/2., 0.1)
    off_y = random.normalvariate((y[0]+y[1])/2., 0.1)
    return off_x, off_y

pop = sim.population(1000, loci=[1], infoFields=['x', 'y'])
pop.setVirtualSplitter(sim.GenotypeSplitter(loci=0, alleles=[[0, 0], [0,1], [1, 1]]))
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5]),
        sim.InitInfo(random.random, infoFields=['x', 'y'])
    ],
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.PyTagger(func=randomMove),
    ]),
    postOps=[
        sim.stat(minOfInfo='x', maxOfInfo='x'),
        sim.PyEval(r"'Range of x: %.2f, %.2f\n' % (minOfInfo['x'], maxOfInfo['x'])")
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(1000, loci=[1], infoFields=['aff', 'numOfAff'])
# define virtual subpopulations by affection sim.status
pop.setVirtualSplitter(sim.AffectionSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5]),
    ],
    preOps=[
        # get affection sim.status for parents
        sim.MaPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4]),
        # set 'aff' of parents
        sim.InfoExec('aff = ind.affected()', exposeInd='ind'),
    ],
        # get number of affected parents for each offspring and store in numOfAff
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.SummaryTagger(mode=sim.SUMMATION, infoFields=['aff', 'numOfAff'])]),
    postOps=[
        # get affection sim.status for offspring
        sim.MaPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4]),
        # calculate mean 'numOfAff' of offspring, for unaffected and affected subpopulations.
        sim.stat(meanOfInfo='numOfAff', subPops=[(0,0), (0,1)], vars=['meanOfInfo_sp']),
        # print mean number of affected parents for unaffected and affected offspring.
        sim.PyEval(r"'sim.MEAN number of affected parents: %.2f (unaff), %.2f (aff)\n' % "
            "(subPop[(0,0)]['meanOfInfo']['numOfAff'], subPop[(0,1)]['meanOfInfo']['numOfAff'])")
    ],
    gen = 5
)

#end_file



#begin_file log/MapPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=2000, loci=2)
sim.initByFreq(pop, [.2, .8])
sim.mapPenetrance(pop, loci=0,
    penetrance={(0,0):0, (0,1):.2, (1,1):.3})
sim.Stat(pop, genoFreq=0, numOfAffected=1, vars='genoNum')
# number of affected individuals
pop.dvars().numOfAffected
# which should be roughly (#01 + #10) * 0.2 + #11 * 0.3
(pop.dvars().genoNum[0][(0,1)] + pop.dvars().genoNum[0][(1,0)]) * 0.2 \
+ pop.dvars().genoNum[0][(1,1)] * 0.3
#end_file

#begin_file log/MaPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(5000, loci=3)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.9] + [0.02]*5)
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.MaPenetrance(loci=0, penetrance=(0.01, 0.2, 0.3)),
        sim.stat(numOfAffected=True, vars='propOfAffected'),
        sim.PyEval(r"'Gen: %d Prevalence: %.1f%%\n' % (gen, propOfAffected*100)"),
    ],
    gen = 5
)
#end_file

#begin_file log/MlPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(5000, loci=3)
sim.initByFreq(pop, [0.2]*5)
# the multi-loci penetrance
sim.mlPenetrance(pop, mode=sim.MULTIPLICATIVE,
    ops = [sim.MaPenetrance(loci=loc,
        penetrance=[0, 0.3, 0.6]) for loc in range(3)])
# count the number of affected individuals.
sim.Stat(pop, numOfAffected=True)
pop.dvars().numOfAffected
#end_file

#begin_file log/PyPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=2000, loci=[1]*2, infoFields=['p', 'smoking'])
pop.setVirtualSplitter(sim.InfoSplitter(field='smoking', values=[0,1]))
# the second parameter gen can be used for varying selection pressure
def penet(geno, smoking):
    #     BB     Bb      bb
    # AA  0.01   0.01    0.01
    # Aa  0.01   0.03    0.03
    # aa  0.01   0.03    0.05
    #
    # geno is (A1 A2 B1 B2)
    if geno[0] + geno[1] == 1 and geno[2] + geno[3] != 0:
        v = 0.03   # case of AaBb
    elif geno[0] + geno[1] == 2 and geno[2] + geno[3] == 1:
        v = 0.03   # case of aaBb
    elif geno[0] + geno[1] ==2 and geno[2] + geno[3] == 2:
        v = 0.05   # case of aabb
    else:                
        v = 0.01   # other cases
    if smoking:
        return v * 2
    else:
        return v

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(alleleFreq=[.5, .5]),
        sim.pyOutput('Calculate prevalence in smoker and non-smokers'),
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        # set smoking status randomly
        sim.InitInfo(lambda : random.randint(0,1), infoFields='smoking'),
        # assign affection status
        sim.PyPenetrance(loci=[0, 1], func=penet),
        sim.stat(numOfAffected=True, subPops=[(0,0),(0,1)], 
            vars='propOfAffected_sp', step=20),
        sim.PyEval(r"'Non-smoker: %.2f%%\tSmoker: %.2f%%\n' % "
            "(subPop[(0,0)]['propOfAffected']*100, subPop[(0,1)]['propOfAffected']*100)",
            step=20)
    ],
    gen = 50
)

#end_file



#begin_file log/PyQuanTrait.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
pop = sim.population(size=5000, loci=2, infoFields=['qtrait1', 'qtrait2', 'age'])
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[40]))
def qtrait(geno, age):
    'Return two traits that depends on genotype and age'
    return random.normalvariate(age * sum(geno), 10), random.randint(0, 10*sum(geno))

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.2, 0.8]),
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        # use random age for simplicity
        sim.InitInfo(lambda:random.randint(20, 75), infoFields='age'),
        sim.PyQuanTrait(loci=(0,1), func=qtrait, infoFields=['qtrait1', 'qtrait2']),
        sim.stat(meanOfInfo=['qtrait1'], subPops=[(0,0), (0,1)],
            vars='meanOfInfo_sp'),
        sim.PyEval(r"'Mean of trait1: %.3f (age < 40), %.3f (age >=40)\n' % "
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(4000, loci=1, infoFields='fitness')
simu = sim.simulator(pop, rep=3)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    preOps=sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):0.98, (1,1):0.97}),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0, step=10),
        sim.PyEval("'Gen:%3d ' % gen", reps=0, step=10),
        sim.PyEval(r"'%.3f\t' % alleleFreq[0][1]", step=10),
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(10000, loci=1)
simu = sim.simulator(pop, rep=3)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):0.98, (1,1):0.97}),
    ]),
    postOps=[
        sim.stat(alleleFreq=0, step=10),
        sim.PyEval("'Gen:%3d ' % gen", reps=0, step=10),
        sim.PyEval(r"'%.3f\t' % alleleFreq[0][1]", step=10),
        sim.pyOutput('\n', reps=-1, step=10)
    ],
    gen = 50
)
#end_file

#begin_file log/MapSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=1000, loci=1, infoFields='fitness')
s1 = .1
s2 = .2
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(alleleFreq=[.2, .8])
    ],
    preOps=sim.MapSelector(loci=0, fitness={(0,0):1-s1, (0,1):1, (1,1):1-s2}),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.PyEval(r"'%.4f\n' % alleleFreq[0][0]", step=100)
    ],
    gen=301
)

#end_file


#begin_file log/MaSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=1000, loci=1, infoFields='fitness')
s1 = .1
s2 = .2
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(alleleFreq=[.2] * 5)
    ],
    preOps=sim.MaSelector(loci=0, fitness=[1-s1, 1, 1-s2]),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.PyEval(r"'%.4f\n' % alleleFreq[0][0]", step=100)
    ],
    gen = 301)
#end_file

#begin_file log/MaSelectorHaploid.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=10000, ploidy=1, loci=[1,1], infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(alleleFreq=[.5, .5])
    ],
    # fitness values for AB, Ab, aB and ab
    preOps=sim.MaSelector(loci=[0,1], fitness=[1, 1, 1, 0.95]),
    matingScheme=sim.RandomSelection(),
    postOps=[
        sim.stat(haploFreq=[0, 1], step=25),
        sim.PyEval(r"'%.3f\t%.3f\t%.3f\t%.3f\n' % (haploFreq[(0,1)][(0,0)],"
                "haploFreq[(0,1)][(0,1)], haploFreq[(0,1)][(1,0)],"
                "haploFreq[(0,1)][(1,1)])", step=25)
    ],
    gen = 100
)
#end_file

#begin_file log/MlSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=10000, loci=2, infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(alleleFreq=[.5, .5])
    ],
    preOps=[
        sim.MlSelector([
            sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}),
            sim.MapSelector(loci=1, fitness={(0,0):1, (0,1):0.9, (1,1):.8}),
            ], mode = sim.ADDITIVE, reps=0),
        sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}, reps=1),
        sim.MapSelector(loci=1, fitness={(0,0):1, (0,1):0.9, (1,1):.8}, reps=2)
    ],
    matingScheme=sim.randomMating(),
    postOps=[
         sim.stat(alleleFreq=[0,1]),
         sim.PyEval(r"'REP %d:\t%.3f\t%.3f\t' % (rep, alleleFreq[0][1], alleleFreq[1][1])"),
         sim.pyOutput('\n', reps=-1),
    ],
    gen = 5
)
#end_file

#begin_file log/PySelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
import random
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=2000, loci=[1]*2, infoFields=['fitness', 'smoking'])
s1 = .02
s2 = .03
# the second parameter gen can be used for varying selection pressure
def sel(geno, smoking):
    #     BB  Bb   bb
    # AA  1   1    1
    # Aa  1   1-s1 1-s2
    # aa  1   1    1-s2
    #
    # geno is (A1 A2 B1 B2)
    if geno[0] + geno[1] == 1 and geno[2] + geno[3] == 1:
        v = 1 - s1  # case of AaBb
    elif geno[2] + geno[3] == 2:
        v = 1 - s2  # case of ??bb
    else:                
        v = 1       # other cases
    if smoking:
        return v * 0.9
    else:
        return v

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(alleleFreq=[.5, .5])
    ],
    preOps=sim.PySelector(loci=[0, 1], func=sel),
    matingScheme=sim.randomMating(),
    postOps=[
        # set smoking status randomly
        sim.InitInfo(lambda : random.randint(0,1), infoFields='smoking'),
        sim.stat(alleleFreq=[0, 1], step=20),
        sim.PyEval(r"'%.4f\t%.4f\n' % (alleleFreq[0][1], alleleFreq[1][1])", step=20)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=2000, loci=1, infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(alleleFreq=[.5, .5])
    ],
    preOps=[
        sim.MaPenetrance(loci=0, penetrance=[0.01, 0.1, 0.2]),
        sim.stat(numOfAffected=True, step=25, vars='propOfAffected'),
        sim.PyEval(r"'Percent of affected: %.3f\t' % propOfAffected", step=50),
        sim.InfoExec('fitness = not ind.affected()', exposeInd='ind')
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.PyEval(r"'%.4f\n' % alleleFreq[0][1]", step=50)
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(size=[5000, 5000], loci=1, infoFields='fitness')
pop.setVirtualSplitter(sim.SexSplitter())
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq(alleleFreq=[.5, .5])
    ],
    preOps=[
        sim.MaSelector(loci=0, fitness=[1, 1, 0.98], subPops=[(0,0), (1,1)]),
        sim.MaSelector(loci=0, fitness=[1, 1, 1], subPops=[(0,1), (1,0)]),
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=[0], subPops=[(0,0), (0,1), (1,0), (1,1)],
            vars='alleleFreq_sp', step=50),
        sim.PyEval(r"'%.4f\t%.4f\t%.4f\t%.4f\n' % "
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
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import trajectory, simulateForwardTrajectory

traj = simulateForwardTrajectory(N=[2000, 4000], fitness=[1, 0.99, 0.98],
    beginGen=0, endGen=100, beginFreq=[0.2, 0.3],
    endFreq=[[0.1, 0.11], [0.2, 0.21]])
traj.plot('log/forwardTrajectory.png', plot_ylim=[0, 0.5], col_sp=['red', 'blue'],
    plot_main='Simulated trajectory (forward-time)')
pop = sim.population(size=[2000, 4000], loci=10, infoFields='fitness')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.8, 0.2], subPops=0),
        sim.InitByFreq([0.7, 0.3], subPops=1),
        sim.pyOutput('Sp0: loc2\tloc5\tSp1: loc2\tloc5\n'),
    ],
    matingScheme=sim.ControlledRandomMating(
        ops=[sim.Recombinator(rates=0.01)],
        loci=5, alleles=1, freqFunc=traj.func()),
    postOps=[
        sim.stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=20),
        sim.PyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
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
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import trajectory, simulateBackwardTrajectory
from math import exp
def Nt(gen):
    'An exponential sim.population growth demographic model.'
    return int((10**4) * exp(.00115 * gen))

def fitness(gen, sp):
    'sim.CONSTANT positive selection pressure.'
    return [1, 1.01, 1.02]

# simulate a trajectory backward in time, from generation 1000
traj = simulateBackwardTrajectory(N=Nt, fitness=fitness, nLoci=2,
     endGen=1000, endFreq=[0.1, 0.2])
traj.plot('log/backTrajectory.png', plot_ylim=[0, 0.3], plot_xlim=[0, 1000],
    col_loc=['red', 'blue'], plot_main='Simulated trajectory (backward-time)')
print 'Trajectory simulated with length %s ' % len(traj.traj)
pop = sim.population(size=Nt(0), loci=[1]*2)
# save trajectory function in the sim.population's local namespace
# so that the sim.PyEval operator can access it.
pop.dvars().traj = traj.func()
pop.evolve(
    initOps=[sim.InitSex()],
    preOps=traj.mutators(loci=[0, 1]),
    matingScheme=sim.ControlledRandomMating(loci=[0, 1], alleles=[1, 1],
        subPopSize=Nt, freqFunc=traj.func()),
    postOps=[
        sim.stat(alleleFreq=[0, 1], begin=500, step=100),
        sim.PyEval(r"'%4d: %.3f (exp: %.3f), %.3f (exp: %.3f)\n' % (gen, alleleFreq[0][1],"
            "traj(gen)[0], alleleFreq[1][1], traj(gen)[1])",
            begin=500, step=100)
    ],
    gen=1001  # evolve 1001 generations to reach the end of generation 1000
)
#end_file


#begin_file log/simuProgress.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import simuProgress
pop = sim.population(10000, loci=[10], infoFields='index')
prog = simuProgress('Setting individual genotype...\n', pop.popSize(), gui=False)
for idx in range(pop.popSize()):
    # do something to each individaul
    pop.individual(idx).index = idx
    prog.update(idx + 1)

#end_file


#begin_file log/viewVars.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True, gui=False)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import viewVars
pop = sim.population([1000, 2000], loci=3)
sim.initByFreq(pop, [0.2, 0.4, 0.4], loci=0)
sim.initByFreq(pop, [0.2, 0.8], loci=2)
sim.Stat(pop, genoFreq=[0, 1, 2], haploFreq=[0, 1, 2],
    alleleFreq=range(3),
    vars=['genoFreq', 'genoNum', 'haploFreq', 'alleleNum_sp'])
viewVars(pop.vars())
#end_file


#begin_file log/varPlotter.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.plotter import varPlotter
pop = sim.population(size=1000, loci=2)
simu = sim.simulator(pop, rep=3)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByValue([1, 2, 2, 1])
    ],
    matingScheme=sim.randomMating(ops=sim.Recombinator(rates=0.01)),
    postOps=[
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
sim.getRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import varPlotter
pop = sim.population(size=1000, loci=1*4)
simu = sim.simulator(pop, rep=3)
simu.evolve(
    initOps=[sim.InitSex()] +
        [sim.InitByFreq([0.1*(x+1), 1-0.1*(x+1)], loci=x) for x in range(4)],
    matingScheme=sim.randomMating(),
    postOps=[
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
sim.getRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import varPlotter
pop = sim.population(size=1000, loci=1*4)
simu = sim.simulator(pop, rep=3)
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
    initOps=[sim.InitSex()]+
        [sim.InitByFreq([0.1*(x+1), 1-0.1*(x+1)], loci=x) for x in range(4)],
    matingScheme=sim.randomMating(),
    postOps=[
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
sim.getRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import scatterPlotter
import random
pop = sim.population([500], infoFields=['x', 'y', 'anc'])
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sim.InfoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
#
def passInfo(x, y, anc):
    'Parental fields will be passed as tuples'
    off_anc = (anc[0] + anc[1])/2.
    off_x = (x[0] + x[1])/2 + random.normalvariate(off_anc - 0.5, 0.1)
    off_y = (y[0] + y[1])/2 + random.normalvariate(0, 0.1)
    return off_x, off_y, off_anc

pop.evolve(
    initOps=[
        sim.InitSex(),
        # random geographic location
        sim.InitInfo(random.random, infoFields=['x', 'y']),
        # anc is 0 or 1
        sim.InitInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.PyTagger(passInfo)]),
    postOps=[
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
sim.getRNG().setSeed(12345)
#end_ignore
import simuPOP as sim
from simuPOP.plotter import histPlotter, qqPlotter, boxPlotter
import random
pop = sim.population([500], infoFields=['x', 'y', 'anc'])
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sim.SexSplitter())

def passInfo(x, y, anc):
    'Parental fields will be passed as tuples'
    off_anc = (anc[0] + anc[1])/2.
    off_x = (x[0] + x[1])/2 + random.normalvariate(off_anc - 0.5, 0.1)
    off_y = (y[0] + y[1])/2 + random.normalvariate(0, 0.1)
    return off_x, off_y, off_anc

pop.evolve(
    initOps=[
        sim.InitSex(),
        # random geographic location
        sim.InitInfo(random.random, infoFields=['x', 'y']),
        # anc is 0 or 1
        sim.InitInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    matingScheme=sim.randomMating(ops=[
        sim.MendelianGenoTransmitter(),
        sim.PyTagger(passInfo)]),
    postOps=[
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
sim.getRNG().setSeed(12345)
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


#begin_file log/randomSample.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.sampling import drawRandomSample
pop = sim.population([2000]*5, loci=1)
# sample from the whole population
sample = drawRandomSample(pop, sizes=500)
print sample.subPopSizes()
# sample from each subpopulation
sample = drawRandomSample(pop, sizes=[100]*5)
print sample.subPopSizes()
#end_file


#begin_file log/caseControlSample.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
from simuPOP.sampling import drawCaseControlSamples
pop = sim.population([10000], loci=5)
sim.initByFreq(pop, [0.2, 0.8])
sim.maPenetrance(pop, loci=2, penetrance=[0.11, 0.15, 0.20])
# draw multiple case control sample
samples = drawCaseControlSamples(pop, cases=500, controls=500, numOfSamples=5)
for sample in samples:
    sim.Stat(sample, association=range(5))
    print ', '.join(['%.6f' % sample.dvars().Allele_ChiSq_p[x] for x in range(5)])

#end_file

#begin_file log/plotPedigree.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12347)
#end_ignore
from simuPOP.sampling import IndexToID, PlotPedigree
pop = sim.population(size=15, loci=5, infoFields=['father_idx', 'mother_idx'], ancGen=2)
pop.evolve(
    preOps=[
        sim.InitSex(),
        sim.InitByFreq([0.7, 0.3]),
    ],
    matingScheme=sim.randomMating(numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
        ops=[sim.MendelianGenoTransmitter(), sim.ParentsTagger()]),
    postOps=sim.MaPenetrance(loci=3, penetrance=(0.1, 0.4, 0.7)),
    gen = 5
)
IndexToID(pop, reset=True)
# three information fields were added
print pop.infoFields()
# save this population for future use
pop.save('log/pedigree.pop')
# draw pedigree
PlotPedigree(pop, filename='log/pedigree.png')
#end_file


#begin_file log/sampleAffectedSibpair.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12347)
#end_ignore
from simuPOP.sampling import drawAffectedSibpairSample, PlotPedigree
pop = sim.loadPopulation('log/pedigree.pop')
sample = drawAffectedSibpairSample(pop, families=2)
PlotPedigree(sample, filename='log/affectedSibpair.png')
#end_file


#begin_file log/sampleNuclearFamily.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12347)
#end_ignore
from simuPOP.sampling import drawNuclearFamilySample, PlotPedigree
pop = sim.loadPopulation('log/pedigree.pop')
sample = drawNuclearFamilySample(pop, families=2, numOffspring=(2,4),
    affectedParents=(1,2), affectedOffspring=(1, 3))
PlotPedigree(sample, filename='log/nuclerFamily.png')
#end_file



#begin_file log/sampleThreeGenFamily.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12347)
#end_ignore
from simuPOP.sampling import drawThreeGenFamilySample, PlotPedigree
pop = sim.loadPopulation('log/pedigree.pop')
sample = drawThreeGenFamilySample(pop, families=2, numOffspring=(1, 3),
    pedSize=(8, 15), numOfAffected=(2, 5))
PlotPedigree(sample, filename='log/threeGenFamily.png')
#end_file


#begin_file log/combinedSampling.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12347)
#end_ignore
from simuPOP.sampling import drawCombinedSample, affectedSibpairSampler, nuclearFamilySampler, PlotPedigree
pop = sim.loadPopulation('log/pedigree.pop')
sample = drawCombinedSample(pop, samplers = [
    affectedSibpairSampler(families=1),
    nuclearFamilySampler(families=1, numOffspring=(2,4), affectedParents=(1,2), affectedOffspring=(1,3))
    ])
PlotPedigree(sample, filename='log/combinedSampling.png')
#end_file


#begin_file log/samplingVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
# create an age-structured population with a disease
import random
pop = sim.population(10000, loci=10, infoFields='age')
sim.initByFreq(pop, [0.3, 0.7])
sim.initInfo(pop, lambda: random.randint(0, 70), infoFields='age')
pop.setVirtualSplitter(sim.InfoSplitter(cutoff=(40, 60), field='age'))
sim.maPenetrance(pop, loci=5, penetrance=(0.1, 0.2, 0.3))
#
from simuPOP.sampling import drawCaseControlSample
sample = drawCaseControlSample(pop, cases=500, controls=500, subPops=[(0,1)])
ageInSample = sample.indInfo('age')
print min(ageInSample), max(ageInSample)
#end_file

#begin_file log/samplingSeparateVSPs.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
# create an age-structured population with a disease
import random
pop = sim.population(10000, loci=10, infoFields='age')
sim.initByFreq(pop, [0.3, 0.7])
sim.initInfo(pop, lambda: random.randint(0, 70), infoFields='age')
pop.setVirtualSplitter(sim.InfoSplitter(cutoff=(20, 40), field='age'))
# different age group has different penetrance
sim.maPenetrance(pop, loci=5, penetrance=(0.1, 0.2, 0.3), subPops=[(0,1)])
sim.maPenetrance(pop, loci=5, penetrance=(0.2, 0.4, 0.6), subPops=[(0,2)])
# count the number of affected individuals in each group
sim.Stat(pop, numOfAffected=True, subPops=[(0,1), (0,2)], vars='numOfAffected_sp')
print pop.dvars((0,1)).numOfAffected, pop.dvars((0,2)).numOfAffected
#
from simuPOP.sampling import drawRandomSample
sample = drawRandomSample(pop, sizes=[500, 500], subPops=[(0,1), (0,2)])
# virtual subpopulations are rearranged to different subpopulations.
print sample.subPopSizes()
#end_file


#begin_file log/reichDemo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
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
    def ins_expansion(gen):
        if gen < G0:
            return N0
        else:
            return N1
    rate = (math.log(N1) - math.log(N0))/G1
    def exp_expansion(gen):
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
sim.getRNG().setSeed(12345)
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
sim.initByFreq(pop, [.2] * 5)
print Ne(pop, loci=[2, 4])
#end_file

#begin_file log/reichEvolve.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
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
    def ins_expansion(gen):
        if gen < G0:
            return N0
        else:
            return N1
    
    def exp_expansion(gen):
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
    pop = sim.population(size=demo_func(0), loci=1, infoFields='fitness')
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitByFreq(loci=0, alleleFreq=spec)
        ],
        matingScheme=sim.randomMating(subPopSize=demo_func),
        postOps=[
            sim.KamMutator(k=k, rates=mu),
            sim.MaSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=[0], step=100),
            sim.PyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
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
     'label': 'sim.MAXIMUM allelic sim.state',
     'allowedTypes': [types.IntType],
     'description': 'sim.MAXIMUM allelic sim.state for a k-allele mutation model',
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
    def ins_expansion(gen):
        if gen < G0:
            return N0
        else:
            return N1
    
    def exp_expansion(gen):
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
    pop = sim.population(size=demo_func(0), loci=1, infoFields='fitness')
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitByFreq(loci=0, alleleFreq=spec)
        ],
        matingScheme=sim.randomMating(subPopSize=demo_func),
        postOps=[
            sim.KamMutator(rate=mu, maxAllele=k),
            sim.MaSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=0, step=100),
            sim.PyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
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
if os.path.isfile('log/simuCDCV.py'):
    out = os.popen('python log/simuCDCV.py -h')
    hlp = open('log/simuCDCV.hlp', 'w')
    print >> hlp, out.read()
    hlp.close()

#end_ignore
#end_file


#begin_file log/randomSeed.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
def simulate():
    pop = sim.population(1000, loci=10, infoFields='age')
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitByFreq([0.5, 0.5]),
            sim.InitInfo(lambda: random.randint(0, 10), infoFields='age')
        ],
        matingScheme=sim.randomMating(),
        finalOps=sim.stat(alleleFreq=0),
        gen=100
    )
    return pop.dvars().alleleFreq[0][0]

seed = sim.getRNG().seed()
random.seed(seed)
print '%.4f' % simulate()
# will yield different result
print '%.4f' % simulate()
sim.setRNG(seed=seed)
random.seed(seed)
# will yield identical result because the same seeds are used
print '%.4f' % simulate()
#end_file


#begin_file log/debug.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
# redirect system stderr
import sys
debugOutput = open('debug.txt', 'w')
old_stderr = sys.stderr
sys.stderr = debugOutput
# start simulation
simu = sim.simulator(sim.population(100, loci=1), rep=5)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.1, 0.9])
    ],
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=0),
        sim.IfElse('alleleNum[0][0] == 0',
            ifOps=[
                # the is None part makes the function return True
                sim.pyOperator(lambda : sim.turnOnDebug("DBG_MUTATOR") is None),
                sim.PointMutator(loci=0, allele=0, inds=0),
            ],
            elseOps=sim.pyOperator(lambda : sim.turnOffDebug("DBG_MUTATOR") is None)),
    ],
    gen = 100
)
# replace standard stdandard error
sys.stderr = old_stderr
debugOutput.close()
print ''.join(open('debug.txt').readlines()[:5])
#begin_ignore
sim.turnOffDebug("DBG_MUTATOR")
#end_ignore
#end_file



#begin_file log/newOperator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
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
                sim.kamMutate(pop, k=2, rates=self.mu1, loci=[i])
            else:
                sim.kamMutate(pop, k=2, rates=self.mu2, loci=[i])
        return True

pop = sim.population(size=10000, loci=[2, 3])
pop.evolve(
    initOps=[ 
        sim.InitSex(),
        sim.InitByFreq([.99, .01], loci=[0, 2, 4]),
        sim.InitByFreq([.8, .2], loci=[1, 3])
    ],
    preOps=dynaMutator(cutoff=.2, mu1=1e-2, mu2=1e-5),
    matingScheme=sim.randomMating(),
    postOps=[
        sim.stat(alleleFreq=range(5), step=10),
        sim.PyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)          
#end_file

#begin_file log/randomMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
from simuPOP import *
#begin_ignore
getRNG().setSeed(12345)
#end_ignore
def randomMating(numOffspring=1., sexMode=RANDOM_SEX,
        ops=MendelianGenoTransmitter(), subPopSize=[],
        subPops=ALL_AVAIL, weight=0, selectionField='fitness'):
    'A basic diploid sexual random mating scheme.'
    return homoMating(
        chooser=RandomParentsChooser(True, selectionField),
        generator=OffspringGenerator(ops, numOffspring, sexMode),
        subPopSize=subPopSize,
        subPops=subPops,
        weight=weight)

#end_file

#begin_file log/sequentialSelfing.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(100, loci=5*3, infoFields='parent_idx')
pop.evolve(
    initOps=[sim.InitByFreq([0.2]*5)],
    preOps=sim.dumper(structure=False, max=5),
    matingScheme=sim.homoMating(
        sim.SequentialParentChooser(),
        sim.OffspringGenerator(ops=[
            sim.SelfingGenoTransmitter(),
            sim.ParentsTagger(infoFields='parent_idx'),
        ])
    ),
    postOps=sim.dumper(structure=False, max=5),
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
sim.getRNG().setSeed(12345)
#end_ignore
def traj(gen):
    return [0.5 + gen * 0.01]

pop = sim.population(1000, loci=[10]*2)
# evolve the sim.population while keeping allele frequency 0.5
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme=sim.homoMating(sim.RandomParentChooser(),
        sim.ControlledOffspringGenerator(loci=5,
            alleles=[0], freqFunc=traj,
            ops = sim.SelfingGenoTransmitter())),
    postOps=[
        sim.stat(alleleFreq=[5, 15]),
        sim.PyEval(r'"%.2f\t%.2f\n" % (alleleFreq[5][0], alleleFreq[15][0])')
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
sim.getRNG().setSeed(12345)
#end_ignore
pop = sim.population(10, loci=[5]*5,
    # one autosome, two sex chromosomes, and two mitochondrial chromosomes
    chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y] + [sim.CUSTOMIZED]*2,
    infoFields=['father_idx', 'mother_idx'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.4] + [0.2]*3)
    ],
    matingScheme=sim.randomMating(ops= [
        sim.Recombinator(rates=0.1),
        sim.MitochondrialGenoTransmitter(),
        sim.ParentsTagger()
    ]),
    postOps=sim.dumper(structure=False),
    gen = 2
)
#end_file

#begin_file log/sexSpecificRec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
from simuPOP import *
#begin_ignore
getRNG().setSeed(12345)
#end_ignore
class sexSpecificRecombinator(pyOperator):
    def __init__(self, intensity=0, rates=0, loci=[], convMode=NO_CONVERSION,
            maleIntensity=0, maleRates=0, maleLoci=[], maleConvMode=NO_CONVERSION,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.Recombinator = Recombinator(rates, intensity, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = Recombinator(maleRates, maleIntensity,
            maleLoci, maleConvMode)
        #
        self.initialized = False
        #
        pyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)
    #
    def transmitGenotype(self, pop, off, dad, mom):
        # Recombinators need to be initialized. Basically, they cache some
        # population properties to speed up genotype transmission.
        if not self.initialized:
            self.Recombinator.initialize(pop)
            self.maleRecombinator.initialize(pop)
            self.initialized = True
        # Form the first homologous copy of offspring.
        self.Recombinator.transmitGenotype(mom, off, 0)
        self.maleRecombinator.transmitGenotype(dad, off, 1)
        return True

pop = population(10, loci=[15]*2, infoFields=['father_idx', 'mother_idx'])
pop.evolve(
    initOps=[
        InitSex(),
        InitByFreq([0.4] + [0.2]*3)
    ],
    matingScheme=randomMating(ops=[
        sexSpecificRecombinator(rates=0.1, maleRates=0),
        ParentsTagger()
    ]),
    postOps=dumper(structure=False),
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
sim.getRNG().setSeed(12345)
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

#begin_file log/PyParentsChooser.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
from random import randint
def randomChooser(pop, sp):
    males = []
    females = []
    # identify males and females in each social rank
    for rank in range(3):
        males.append([x for x in pop.individuals(sp) \
            if x.sex() == sim.MALE and x.rank == rank])
        females.append([x for x in pop.individuals(sp) \
            if x.sex() == sim.FEMALE and x.rank == rank])
    #
    while True:
        # choose a rank randomly
        rank = int(pop.individual(randint(0, pop.subPopSize(sp) - 1), sp).rank)
        yield males[rank][randint(0, len(males[rank]) - 1)], \
            females[rank][randint(0, len(females[rank]) - 1)]

def setRank(off, dad):
    'The rank of offspring can increase or drop to zero randomly'
    off.rank = (dad.rank + randint(-1, 1)) % 3

pop = sim.population(size=[1000, 2000], loci=1, infoFields='rank')
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitInfo(lambda : randint(0, 2), infoFields='rank')
    ],
    matingScheme=sim.homoMating(
        sim.PyParentsChooser(randomChooser),
        sim.OffspringGenerator(ops=sim.MendelianGenoTransmitter())
    ),
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
sim.getRNG().setSeed(12345)
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
        [x for x in range(pop.popSize()) if pop.individual(x).sex() == sim.MALE],
        [x for x in range(pop.popSize()) if pop.individual(x).sex() == sim.FEMALE])
    while True:
        # return indexes of parents repeatedly
        yield pc.chooseParents()

pop = sim.population(100, loci=1)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitByFreq([0.5, 0.5])
    ],
    matingScheme=sim.homoMating(sim.PyParentsChooser(parentsChooser),
        sim.OffspringGenerator(ops=sim.MendelianGenoTransmitter())),
    gen = 100
)
#end_file



#begin_file log/ageStructured.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
#end_ignore
import simuPOP as sim
#begin_ignore
sim.getRNG().setSeed(12345)
#end_ignore
import random
N = 10000
pop = sim.population(N, loci=1, infoFields=['age', 'ind_id', 'father_id', 'mother_id'])
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[20, 50, 75]))
def demoModel(gen, pop):
    '''A demographic model that keep a constant supply of new individuals'''
    # number of individuals that will die
    sim.Stat(pop, popSize=True, subPops=[(0,3)])
    # individuals that will be kept, plus some new guys.
    return pop.popSize() - pop.dvars().popSize + N / 75

def pene(geno, age, ind):
    'Define an age-dependent penetrance function'
    # this disease does not occur in children
    if age < 16:
        return 0
    # if an individual is already affected, keep so
    if ind.affected():
        return 1
    # the probability of getting disease increases with age
    return (0., 0.001*age, 0.001*age)[sum(geno)]

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


pop.evolve(
    initOps=[
        sim.InitSex(),
        # random assign age
        sim.InitInfo(lambda: random.randint(0, 75), infoFields='age'),
        # random genotype
        sim.InitByFreq([0.5, 0.5]),
        # assign an unique ID to everyone.
        sim.IdTagger(),
        sim.pyOutput('Prevalence of disease in each age group:\n'),
    ],
    # increase the age of everyone by 1 before mating.
    preOps=sim.InfoExec('age += 1'),
    matingScheme=sim.heteroMating([
        # all individuals with age < 75 will be kept. Note that
        # CloneMating will keep individual sex, affection status and all
        # information fields (by default).
        sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
        # only individuals with age between 20 and 50 will mate and produce
        # offspring. The age of offspring will be zero.
        sim.randomMating(ops=[
            sim.IdTagger(),                   # give new born an ID
            sim.PedigreeTagger(),             # track parents of each individual
            sim.MendelianGenoTransmitter(),   # transmit genotype
        ],
        numOffspring=(sim.UNIFORM_DISTRIBUTION, 1, 3),
        subPops=[(0,1)]),],
        subPopSize=demoModel),
    # number of individuals?
    postOps=[
        sim.PyPenetrance(func=pene, loci=0),
        sim.pyOperator(func=outputStat, step=20)
    ],
    gen = 200
)

# draw two pedigrees from the last age-structured population
from simuPOP import sampling
sample = sampling.drawNuclearFamilySample(pop, families=2, numOffspring=(2,3),
    affectedParents=(1,2), affectedOffspring=(1,3))
sim.Dump(sample)

#end_file
