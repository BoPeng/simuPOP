#!/usr/bin/env python

#begin_file log/importSimuPOP.py
import simuOpt
simuOpt.setOptions(optimized=False, alleleType='long', quiet=False)
from simuPOP import *
#end_file

#begin_file log/standard.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
#expect_error
from simuPOP import *
pop = population(10, loci=2)
pop.locusPos(10)
pop.individual(20).setAllele(1, 0)
#end_file

#begin_file log/simpleExample.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=1000, loci=2)
simu = simulator(pop, randomMating(ops=recombinator(rates=0.01)), rep=3)
simu.evolve(
    preOps = [
        initSex(),
        initByValue([1, 2, 2, 1])
    ],
    ops = [
        stat(LD=[0, 1]),
        pyEval(r"'%.2f\t' % LD[0][1]", step=10),
        pyOutput('\n', reps=-1, step=10)
    ],
    gen=100
)
#end_file

#begin_file log/help.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
help(population.addInfoFields)
#end_file

#begin_file log/absIndex.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[10, 20], loci=[5, 7])
print pop.chromLocusPair(7)
print pop.absLocusIndex(1, 1)
print pop.absIndIndex(2, 1)
print pop.subPopIndPair(25)
#end_file

#begin_file log/iterator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=2, loci=[5, 6])
InitByFreq(pop, [0.2, 0.3, 0.5])
for ind in pop.individuals():
    for loc in range(pop.chromBegin(1), pop.chromEnd(1)):
        print ind.allele(loc),
    print 
#end_file

#begin_file log/carray.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=2, loci=[3, 4])
InitByFreq(pop, [0.3, 0.5, 0.2])
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population([100]*2, loci=1)
InitByFreq(pop, [0, 0.2, 0.8], subPops=0)
InitByFreq(pop, [0.2, 0.8], subPops=1)
Stat(pop, alleleFreq=0, vars=['alleleFreq_sp'])
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2, 3], ploidy=2, loci=[5, 10],
    lociPos=range(0, 5) + range(0, 20, 2), chromNames=['Chr1', 'Chr2'],
    alleleNames=['A', 'C', 'T', 'G'])
# access genotypic information from the population
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
pop = population(loci=[2, 3], lociPos=[3, 1, 1, 3, 2],
    lociNames=['loc%d' % x for x in range(5)])
pop.lociPos()
pop.lociNames()
#end_file

#begin_file log/haplodiploid.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2,5], ploidy=Haplodiploid, loci=[3, 5])
InitByFreq(pop, [0.3, 0.7])
Dump(pop)
#end_file

#begin_file log/chromType.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=6, ploidy=2, loci=[3, 3, 6, 4, 4, 4],
    chromTypes=[Autosome]*2 + [ChromosomeX, ChromosomeY] + [Customized]*2)
InitByFreq(pop, [0.3, 0.7])
Dump(pop, structure=False) # does not display genotypic structure information
#end_file

#begin_file log/infoField.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(10, loci=[20], ancGen=1,
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, randomMating(ops=recombinator(rates=0.01)))
simu.evolve(
    preOps = [
        initSex(),
        initByValue([0]*20+[1]*20)
    ],
    ops = parentsTagger(),
    gen = 1
)
pop = simu.extract(0)
pop.indInfo('mother_idx')  # mother of all offspring
ind = pop.individual(0)
mom = pop.ancestor(ind.intInfo('mother_idx'), 1)
print ind.genotype(0)
print mom.genotype(0)
print mom.genotype(1)
#end_file

#begin_file log/individual.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population([5, 4], loci=[2, 5], infoFields='x')
# get an individual
ind = pop.individual(3)
ind.ploidy()            # access to genotypic structure
ind.numChrom()
ind.affected()
ind.setAffected(True)   # access affection status,
ind.sex()               # sex,
ind.setInfo(4, 'x')     # and information fields
ind.info('x')
ind.intInfo(0)          # obtain the value of 'x' as an integer.
#end_file

#begin_file log/individual_genotype.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population([2, 1], loci=[2, 5])
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[3, 4, 5], ploidy=1, loci=1, infoFields='x')
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
Dump(pop, width=2, structure=False)
#end_file

#begin_file log/subPopName.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[3, 4, 5], subPopNames=['x', 'y', 'z'])
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population(size=[200, 400], loci=[30], infoFields='x')
# assign random information fields
InitSex(pop)
InitInfo(pop, lambda: random.randint(0, 3), infoFields='x')
# define a virtual splitter by sex
pop.setVirtualSplitter(sexSplitter())
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 1])    # Size of VSP 0 in subpopulation 0
# define a virtual splitter by information field 'x'
pop.setVirtualSplitter(infoSplitter(field='x', values=[0, 1, 2, 3]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # Size of VSP 0 in subpopulation 0
pop.subPopSize([1, 0])    # Size of VSP 0 in subpopulation 1
#end_file

#begin_file log/virtualSubPop.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population(10, loci=[2, 3], infoFields='Sex')
InitSex(pop)
pop.setVirtualSplitter(sexSplitter())
# initialize male and females with different genotypes. Set initSex
# to False because this operator will by default also initialize sex.
InitByValue(pop, [[0]*5, [1]*5], subPops=([0, 0], [0, 1]))
# set Sex information field to 0 for all males, and 1 for all females
pop.setIndInfo([Male], 'Sex', [0, 0])
pop.setIndInfo([Female], 'Sex', [0, 1])
# Print individual genotypes, followed by values at information field Sex
Dump(pop, structure=False)
#end_file


#begin_file log/advancedVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population(size=[2000, 4000], loci=[30], infoFields='x')
# assign random information fields
InitSex(pop)
InitInfo(pop, lambda: random.randint(0, 3), infoFields='x')
#
# 1, use a combined splitter
pop.setVirtualSplitter(combinedSplitter(splitters = [
    sexSplitter(),
    infoSplitter(field='x', values=[0, 1, 2, 3])
]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # Male
pop.subPopSize([1, 4])    # individuals in sp 1 with value 2 at field x
#
# use a product splitter that defines additional VSPs by sex and info
pop.setVirtualSplitter(productSplitter(splitters = [
    sexSplitter(),
    infoSplitter(field='x', values=[0, 1, 2, 3])
]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # Male with value 1 in sp 0
pop.subPopSize([1, 5])    # Female with value 1 in sp 1
#
# use a combined splitter to join VSPs defined by a
# product splitter
pop.setVirtualSplitter(combinedSplitter([
    productSplitter([
        sexSplitter(),
        infoSplitter(field='x', values=[0, 1, 2, 3])])],
    vspMap = [[0,1,2], [4,5,6], [7]]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # Male with value 0, 1, 2 at field x
pop.subPopSize([1, 1])    # Female with value 0, 1 or 2 at field x
#end_file


#begin_file log/accessIndividual.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
# create a population with two generations. The current generation has values
# 0-9 at information field x, the parental generation has values 10-19.
pop = population(size=[5, 5], loci=[2, 3], infoFields='x', ancGen=1)
pop.setIndInfo(range(11, 20), 'x')
pop1 = pop.clone()
pop1.setIndInfo(range(10), 'x')
pop.push(pop1)
#
ind = pop.individual(5)       # using absolute index
ind.info('x')
# use a for loop, and relative index
for idx in range(pop.subPopSize(1)):
    print pop.individual(idx, 1).info('x'),

# It is usually easier to use an iterator
for ind in pop.individuals(1):
    print ind.info('x'),

# Access individuals in VSPs
pop.setVirtualSplitter(infoSplitter(cutoff=[3, 7], field='x'))
for ind in pop.individuals([1, 1]):
    print ind.info('x'),

# Access individuals in ancetral generations
pop.ancestor(5, 1).info('x')        # absolute index
pop.ancestor(0, 1, 1).info('x')     # relative index
# Or make ancestral generation the current generation and use 'individual'
pop.useAncestralGen(1)
pop.individual(5).info('x')         # absolute index
pop.individual(0, 1).info('x')      # relative index
# 'ancestor' can still access the 'present' (generation 0) generation
pop.ancestor(5, 0).info('x')
# access individual by ID
pop.addInfoFields('ind_id')
TagID(pop)
[ind.intInfo('ind_id') for ind in pop.individuals()]
# access individual by ID. Note that individual 12 is in the parental generation
pop.indByID(12).info('x')
#end_file

#begin_file log/batchAccess.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population(size=[4, 6], loci=2, infoFields='x')
pop.setIndInfo([random.randint(0, 10) for x in range(10)], 'x')
pop.indInfo('x')
pop.setGenotype([0, 1, 2, 3], 0)
pop.genotype(0)
pop.setVirtualSplitter(infoSplitter(cutoff=[3], field='x'))
pop.setGenotype([0])    # clear all values
pop.setGenotype([5, 6, 7], [1, 1])
pop.indInfo('x', 1)
pop.genotype(1)
#end_file

#begin_file log/popInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(10)
pop.setInfoFields(['a', 'b'])
pop.addInfoFields('c')
pop.addInfoFields(['d', 'e'])
pop.infoFields()
#
cIdx = pop.infoIdx('c')
eIdx = pop.infoIdx('e')
# information fields can be accessed in batch mode
pop.setIndInfo([1], cIdx)
# as well as individually.
for ind in pop.individuals():
    ind.setInfo(ind.info(cIdx) + 1, eIdx)

print pop.indInfo(eIdx)
#end_file

#begin_file log/ancestralPop.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(500, loci=1), randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [
        # start recording ancestral generations at generation 18
        setAncestralDepth(2, at=[-2]),
        stat(alleleFreq=0, begin=-3),
        pyEval(r"'%.3f\n' % alleleFreq[0][0]", begin=-3)
    ],
    gen = 20
)
pop = simu.population(0)
# start from current generation
for i in range(pop.ancestralGens(), -1, -1):
  pop.useAncestralGen(i)
  Stat(pop, alleleFreq=0)
  print '%d   %.3f' % (i, pop.dvars().alleleFreq[0][0])

# restore to the current generation  
pop.useAncestralGen(0)  
#end_file

#begin_file log/addRemoveLoci.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(10, loci=3, chromNames=['chr1'])
# 1 1 1, 
pop.setGenotype([1])
# 1 1 1, 0 0 0
pop.addChrom(lociPos=[0.5, 1, 2], lociNames=['rs1', 'rs2', 'rs3'],
    chromName='chr2')
pop1 = population(10, loci=3, chromNames=['chr3'],
    lociNames=['rs4', 'rs5', 'rs6'])
# 2 2 2,
pop1.setGenotype([2])
# 1 1 1, 0 0 0, 2 2 2
pop.addChromFrom(pop1)
# 1 1 1, 0 0 0, 2 0 2 2 0
pop.addLoci(chrom=[2, 2], pos=[1.5, 3.5], lociNames=['rs7', 'rs8'])
# 1 1 1, 0 0 0, 2 0 2 0
pop.removeLoci([8])
Dump(pop)
#end_file

#begin_file log/recode.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(5, loci=[5], alleleNames=['A', 'T', 'C', 'G'])
InitByFreq(pop, [0.2, 0.3, 0.4, 0.1])
Dump(pop, structure=False)
print pop.genotype()
pop.recodeAlleles([0, 3, 1, 2], alleleNames=['A', 'C', 'G', 'T'])
Dump(pop, structure=False)
print pop.genotype()
#end_file

#begin_file log/extract.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population(size=[10, 10], loci=[5, 5], infoFields=['x', 'y'])
InitByValue(pop, range(10))
pop.setIndInfo([-1]*4 + [0]*3 + [-1]*3 + [2]*4 + [-1]*3 + [1]*4, 'x')
pop1 = pop.extract(field='x', loci=[1, 2, 3, 6, 7], infoFields='x')
Dump(pop1, structure=False)
#end_file

#begin_file log/popVars.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from pprint import pprint
pop = population(100, loci=2)
InitByFreq(pop, [0.3, 0.7])
print pop.vars()    # No variable now
pop.dvars().myVar = 21
print pop.vars()
Stat(pop, popSize=1, alleleFreq=0)
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=1),
    randomMating(), 5)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [
        stat(alleleFreq=0),
        terminateIf('len(alleleFreq[0]) == 1')
    ]
)
#end_file

#begin_file log/savePop.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(100, loci=5, chromNames=['chrom1'])
pop.dvars().name = 'my population'
pop.save('sample.pop')
pop1 = LoadPopulation('sample.pop')
pop1.chromName(0)
pop1.dvars().name
#begin_ignore
import os
os.remove('sample.pop')
#end_ignore
#end_file

#begin_file log/stageAndGen.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=[20]), randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.2, 0.8])
    ],
    ops = [
        stat(alleleFreq=0, begin=80, step=10),
        pyEval(r"'After gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        pyEval(r"'Around gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1], stage=PrePostMating)
    ],
    postOps = [savePopulation(output='sample.pop')],
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=[20]), randomMating())
simu.evolve(
    preOps = initByFreq([0.2, 0.8]),
    ops = [
        stat(alleleFreq=0, begin=80, step=10),
        pyEval(r"'After gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        pyEval(r"'Around gen %d: alleleFreq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1], stage=PrePostMating)
    ],
    postOps = [savePopulation(output='sample.pop')],
    gen=100,
    dryrun = True
)
#end_file

#begin_file log/replicate.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=[20]), randomMating(), 5)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.2, 0.8])
    ],
    ops = [
        stat(alleleFreq=0, step=10),
        pyEval('gen', step=10, reps=0),
        pyEval(r"'\t%.2f' % alleleFreq[0][0]", step=10, reps=(0, 2, -1)),
        pyOutput('\n', step=10, reps=-1)
    ],
    gen=30,
)
#end_file

#begin_file log/output.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(size=1000, loci=2),
    randomMating(ops=recombinator(rates=0.01)), rep=3)
simu.evolve(
    preOps = [
        initSex(),
        initByValue([1, 2, 2, 1])
    ],
    ops = [
        stat(LD=[0, 1]),
        pyEval(r"'%.2f\t' % LD[0][1]", step=20, output='>>LD.txt'),
        pyOutput('\n', reps=-1, step=20, output='>>LD.txt'),
        pyEval(r"'%.2f\t' % R2[0][1]", output='R2.txt'),
        pyEval(r"'%.2f\t' % LD[0][1]", step=20, output="!'>>LD_%d.txt' % rep"),
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
from simuPOP import *
GetRNG().setSeed(12345)
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
simu = simulator(population(size=1000, loci=2),
    randomMating(ops=recombinator(rates=0.01)))
simu.evolve(
    preOps = [
        initSex(),
        initByValue([1, 2, 2, 1])
    ],
    ops = [
        stat(LD=[0, 1]),
        pyEval(r"'LD: %d, %.2f' % (gen, LD[0][1])", step=20,
            output=logger.info),   # send LD to console and a logfile
        pyEval(r"'R2: %d, %.2f' % (gen, R2[0][1])", step=20,
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(size=10000, loci=2),
    randomMating(ops=[
        mendelianGenoTransmitter(end=29),
        recombinator(rates=0.01, begin=30),
    ])
)
simu.evolve(
    preOps = [
        initSex(),
        initByValue([1, 2, 2, 1])
    ],
    ops = [
        stat(LD=[0, 1]),
        pyEval(r"'gen %d, LD: %.2f\n' % (gen, LD[0][1])", step=20)
    ],
    gen=100
)
#end_file

#begin_file log/hybrid.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def myPenetrance(geno):
    'A three-locus heterogeneity penetrance model'
    if sum(geno) < 2:
        return 0
    else:
        return sum(geno)*0.1

simu = simulator(population(1000, loci=[20]*3), randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.8, 0.2])
    ],
    ops = [
        pyPenetrance(func=myPenetrance, loci=[10, 30, 50]),
        stat(numOfAffected=True),
        pyEval(r"'%d: %d\n' % (gen, numOfAffected)")
    ],
    gen = 5
)
#end_file

#begin_file log/pyOperator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def dynaMutator(pop, param):
    '''This mutator mutates commom loci with low mutation rate and rare
    loci with high mutation rate, as an attempt to raise allele frequency
    of rare loci to an higher level.'''
    # unpack parameter
    (cutoff, mu1, mu2) = param;
    Stat(pop, alleleFreq=range(pop.totNumLoci()))
    for i in range(pop.totNumLoci()):
        # Get the frequency of allele 1 (disease allele)
        if pop.dvars().alleleFreq[i][1] < cutoff:
            KamMutate(pop, k=2, rates=mu1, loci=[i])
        else:
            KamMutate(pop, k=2, rates=mu2, loci=[i])
    return True

simu = simulator(population(size=10000, loci=[2, 3]),
    randomMating())
simu.evolve(
    preOps = [ 
        initSex(),
        initByFreq([.99, .01], loci=[0, 2, 4]),
        initByFreq([.8, .2], loci=[1, 3])
    ],
    ops = [ 
        pyOperator(func=dynaMutator, param=(.2, 1e-2, 1e-5), stage=PreMating),
        stat(alleleFreq=range(5), step=10),
        pyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)                
#end_file

#begin_file log/pyDuringMatingOperator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def rejectInd(off):
    'reject an individual if it off.allele(0) == 1'
    return off.allele(0) == 0

simu = simulator(population(size=100, loci=1),
    randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [ 
        pyOperator(func=rejectInd, stage=DuringMating, offspringOnly=True),
    ],
    gen = 1
)
# You should see no individual with allele 1 at locus 0, ploidy 0.
simu.population(0).genotype()[:20]
#end_file

#begin_file log/newOperator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
class dynaMutator(pyOperator):
    '''This mutator mutates commom loci with low mutation rate and rare
    loci with high mutation rate, as an attempt to raise allele frequency
    of rare loci to an higher level.'''
    def __init__(self, cutoff, mu1, mu2, *args, **kwargs):
        self.cutoff = cutoff
        self.mu1 = mu1
        self.mu2 = mu2
        pyOperator.__init__(self, func=self.mutate, *args, **kwargs)
    #
    def mutate(self, pop):
        Stat(pop, alleleFreq=range(pop.totNumLoci()))
        for i in range(pop.totNumLoci()):
            # Get the frequency of allele 1 (disease allele)
            if pop.dvars().alleleFreq[i][1] < self.cutoff:
                KamMutate(pop, k=2, rates=self.mu1, loci=[i])
            else:
                KamMutate(pop, k=2, rates=self.mu2, loci=[i])
        return True

simu = simulator(population(size=10000, loci=[2, 3]),
    randomMating())
simu.evolve(
    preOps = [ 
        initSex(),
        initByFreq([.99, .01], loci=[0, 2, 4]),
        initByFreq([.8, .2], loci=[1, 3])
    ],
    ops = [ 
        dynaMutator(cutoff=.2, mu1=1e-2, mu2=1e-5, stage=PreMating),
        stat(alleleFreq=range(5), step=10),
        pyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)          
#end_file

#begin_file log/funcInitByFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def InitByFreq(pop, *args, **kwargs):
    initByFreq(*args, **kwargs).apply(pop)

pop = population(1000, loci=[2,3])
InitByFreq(pop, [.2, .3, .5])
#end_file

#begin_file log/migrSize.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=[500, 1000], infoFields='migrate_to'),
    randomMating())

simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end_file

#begin_file log/migrFixedSize.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=[500, 1000], infoFields='migrate_to'),
    randomMating(subPopSize=[500, 1000]))
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True, stage=PrePostMating),
        pyEval(r'"%s\n" % subPopSize', stage=PrePostMating)
    ],
    gen = 3
)
#end_file

#begin_file log/demoFunc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def demo(gen, oldSize=[]):
    return [500 + gen*10, 1000 + gen*10]

simu = simulator(
    population(size=[500, 1000], infoFields='migrate_to'),
    randomMating(subPopSize=demo))
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end_file

#begin_file log/numOff.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def checkNumOffspring(ms):
    '''Check the number of offspring for each family using
       information field father_idx
    '''
    simu = simulator(
        population(size=[30], infoFields=['father_idx', 'mother_idx']),
        matingScheme=ms)
    simu.evolve(
        preOps = initSex(),
        ops=[parentsTagger()],
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
checkNumOffspring(randomMating(numOffspring=2))
# Case 2: Use a Python function
import random
def func(gen):
    return random.randint(5, 8)

checkNumOffspring(randomMating(numOffspring=func))
# Case 3: A geometric distribution
checkNumOffspring(randomMating(numOffspring=(GeometricDistribution, 0.3)))
# Case 4: A Possition distribution
checkNumOffspring(randomMating(numOffspring=(PoissonDistribution, 3)))
# Case 5: A Binomial distribution
checkNumOffspring(randomMating(numOffspring=(BinomialDistribution, 0.1, 10)))
# Case 6: A uniform distribution
checkNumOffspring(randomMating(numOffspring=(UniformDistribution, 2, 6)))
#end_file

#begin_file log/sexMode.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def checkSexMode(ms):
    '''Check the assignment of sex to offspring'''
    simu = simulator(
        population(size=[40]),
        matingScheme=ms)
    simu.evolve(preOps = initSex(), ops=[], gen=1)
    # return individual sex as a string
    return ''.join([ind.sexChar() for ind in simu.population(0).individuals()])

# Case 1: NoSex (all male, randomMating will not continue)
checkSexMode(randomMating(sexMode=NoSex))
# Case 2: RandomSex (Male/Female with probability 0.5)
checkSexMode(randomMating(sexMode=RandomSex))
# Case 3: ProbOfMale (Specify probability of male)
checkSexMode(randomMating(sexMode=(ProbOfMale, 0.8)))
# Case 4: NumOfMale (Specify number of male in each family)
checkSexMode(randomMating(numOffspring=3, sexMode=(NumOfMale, 1)))
# Case 5: NumOfFemale (Specify number of female in each family)
checkSexMode(randomMating(
    numOffspring=(UniformDistribution, 4, 6),
    sexMode=(NumOfFemale, 2))
)
#end_file

#begin_file log/monogamous.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(20, infoFields=['father_idx', 'mother_idx']),
    monogamousMating(numOffspring=2, sexMode=(NumOfMale, 1)))
simu.evolve(
    preOps = [initSex(sex=(Male, Female))],
    ops = [parentsTagger()],
    gen = 5
)
pop = simu.extract(0)
[ind.sex() for ind in pop.individuals()]
[ind.intInfo('father_idx') for ind in pop.individuals()]
[ind.intInfo('mother_idx') for ind in pop.individuals()]
# count the number of distinct parents
len(set(pop.indInfo('father_idx')))
len(set(pop.indInfo('mother_idx')))
#end_file

#begin_file log/polygamous.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, infoFields=['father_idx', 'mother_idx']),
    polygamousMating(polySex=Male, polyNum=2))
simu.evolve(
    preOps = initSex(),
    ops = [parentsTagger()],
    gen = 5
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals()][:20]
[ind.intInfo('mother_idx') for ind in pop.individuals()][:20]
#end_file

#begin_file log/randomSelection.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, ploidy=1, loci=[5, 5], ancGen=1,
    infoFields='parent_idx'),
    randomSelection())
simu.evolve(
    preOps = [initByFreq([0.3, 0.7])],
    ops = [parentsTagger(infoFields='parent_idx')],
    gen = 5
)
pop = simu.extract(0)
ind = pop.individual(0)
par = pop.ancestor(ind.intInfo('parent_idx'), 1)
print ind.sex(), ind.genotype()
print par.sex(), par.genotype()
#end_file

#begin_file log/alphaMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(1000, loci=5, 
    infoFields=['father_idx', 'mother_idx', 'fitness']),
    alphaMating(alphaSex=Male, alphaNum=2))
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [parentsTagger(),
        maSelector(loci=0, fitness=[0.8, 0.8, 1]),
        stat(alleleFreq=0),
        pyEval(r'"%.2f\n" % alleleFreq[0][1]', step=5)
    ],
    gen = 20,
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals()][:10]
[ind.intInfo('mother_idx') for ind in pop.individuals()][:10]
#end_file

#begin_file log/haplodiploidMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(10, ploidy=Haplodiploid, loci=[5, 5],
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(sexSplitter())
simu = simulator(pop, haplodiploidMating())
simu.evolve(
    preOps = [
        initSex(),
        initByValue([0]*10, subPops=[(0, 0)]),
        initByValue([1]*10+[2]*10, subPops=[(0, 1)])
    ],
    ops = [parentsTagger(),
        dumper(structure=False, stage=PrePostMating)],
    gen = 1
)
#end_file

#begin_file log/selfMating.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(20, loci=8)
# every chromosomes are different. :-)
for idx, ind in enumerate(pop.individuals()):
    ind.setGenotype([idx*2], 0)
    ind.setGenotype([idx*2+1], 1)

simu = simulator(pop, selfMating(ops=recombinator(rates=0.01)))
simu.evolve(
    ops = [],
    gen = 1
)
Dump(simu.population(0), width=3, structure=False, max=10)
#end_file

#begin_file log/heteroMatingSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[1000, 1000], loci=2,
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, heteroMating(
    [randomMating(numOffspring=2, subPop=0),
     randomMating(numOffspring=4, subPop=1)
    ])
)
simu.evolve(
    preOps = initSex(),
    ops= [parentsTagger()],
    gen=10
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals(0)][:10]
[ind.intInfo('father_idx') for ind in pop.individuals(1)][:10]
#end_file

#begin_file log/heteroMatingVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[1000], loci=2,
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(proportionSplitter([0.2, 0.8]))
simu = simulator(pop, heteroMating(
    matingSchemes = [
        selfMating(subPop=(0, 0)),
        randomMating(subPop=(0, 1))
    ])
)
simu.evolve(
    preOps = initSex(),
    ops= [parentsTagger()],
    gen = 10
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals(0)][:15]
[ind.intInfo('mother_idx') for ind in pop.individuals(0)][:15]
#end_file

#begin_file log/heteroMatingWeight.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[1000], loci=2,
    infoFields='mark')
pop.setVirtualSplitter(rangeSplitter([[0, 500], [200, 1000]]))
def markOff(param):
    '''define a Python during mating operator that marks
       individual information field 'mark'
    '''
    def func(off, param):
        off.setInfo(param, 'mark')
        return True
    return pyOperator(func=func, param=param, stage=DuringMating,
        offspringOnly=True)

simu = simulator(pop, heteroMating(
    matingSchemes = [
        randomMating(subPop=0, weight=-0.5, ops=[markOff(0), mendelianGenoTransmitter()]),
        randomMating(subPop=(0, 0), weight=2, ops=[markOff(1), mendelianGenoTransmitter()]),
        randomMating(subPop=(0, 1), weight=3, ops=[markOff(2), mendelianGenoTransmitter()])
    ])
)
simu.evolve(
    preOps = initSex(),
    ops= [],
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def randomMating(numOffspring = 1., sexMode = RandomSex,
        ops = mendelianGenoTransmitter(), subPopSize = [],
        subPop = (), weight = 0, selectionField = 'fitness'):
    'A basic diploid sexual random mating scheme.'
    return homoMating(
        chooser = randomParentsChooser(True, selectionField),
        generator = offspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)
#end_file

#begin_file log/sequentialSelfing.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=5*3, infoFields='parent_idx'),
    homoMating(sequentialParentChooser(),
    offspringGenerator(ops=selfingGenoTransmitter())))
simu.evolve(
    preOps = [initByFreq([0.2]*5)],
    ops = [
        parentsTagger(infoFields='parent_idx'),
        dumper(structure=False, stage=PrePostMating, max=5)],
    gen = 1
)
#end_file

#begin_file log/controlledOffGenerator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def traj(gen):
    return [0.5 + gen * 0.01]

simu = simulator(population(1000, loci=[10]*2),
    homoMating(randomParentChooser(),
        controlledOffspringGenerator(loci=5,
            alleles=[0], freqFunc=traj,
            ops = selfingGenoTransmitter()))
)
# evolve the population while keeping allele frequency 0.5
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [stat(alleleFreq=[5, 15]),
        pyEval(r'"%.2f\t%.2f\n" % (alleleFreq[5][0], alleleFreq[15][0])')],
    gen = 5
)
#end_file

#begin_file log/mitochondrial.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(10, loci=[5]*5,
    # one autosome, two sex chromosomes, and two mitochondrial chromosomes
    chromTypes=[Autosome, ChromosomeX, ChromosomeY] + [Customized]*2,
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, randomMating(ops= [
    recombinator(rates=0.1),
    mitochondrialGenoTransmitter()]))
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.4] + [0.2]*3)
    ],
    ops = [
        parentsTagger(),
        dumper(structure=False),
    ],
    gen = 2
)
#end_file

#begin_file log/sexSpecificRec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
class sexSpecificRecombinator(pyOperator):
    def __init__(self, intensity=0, rates=0, loci=[], convMode=NoConversion,
            maleIntensity=0, maleRates=0, maleLoci=[], maleConvMode=NoConversion,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.recombinator = recombinator(rates, intensity, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = recombinator(maleRates, maleIntensity,
            maleLoci, maleConvMode)
        #
        self.initialized = False
        #
        pyOperator.__init__(self, func=self.transmitGenotype,
            stage=DuringMating, *args, **kwargs)
    #
    def transmitGenotype(self, pop, off, dad, mom):
        # Recombinators need to be initialized. Basically, they cache some
        # population properties to speed up genotype transmission.
        if not self.initialized:
            self.recombinator.initialize(pop)
            self.maleRecombinator.initialize(pop)
            self.initialized = True
        # Form the first homologous copy of offspring.
        self.recombinator.transmitGenotype(mom, off, 0)
        self.maleRecombinator.transmitGenotype(dad, off, 1)
        return True

pop = population(10, loci=[15]*2, infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, randomMating(
    ops = sexSpecificRecombinator(rates=0.1, maleRates=0)))
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.4] + [0.2]*3)
    ],
    ops=[
        parentsTagger(),
        dumper(structure=False),
    ],
    gen = 2
)
#end_file

#begin_file log/infoChooser.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(100, loci=[10],
    infoFields=['father_idx', 'mother_idx', 'sibling'])
pop.setVirtualSplitter(sexSplitter())
def locate_sibling(pop):
    '''The population is arranged as MFMFMFMF... where MF are siblings, so the
    sibling of males are 1, 3, 5, .. and the slibling of females are 0, 2, 4, ...
    '''
    pop.setIndInfo([2*x+1 for x in range(pop.popSize()/2)], 'sibling', (0, 0))
    pop.setIndInfo([2*x for x in range(pop.popSize()/2)], 'sibling', (0, 1))

simu = simulator(pop, consanguineousMating(func=locate_sibling, infoFields='sibling',
    numOffspring=2, sexMode=(NumOfMale, 1)))
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.2, 0.8])
    ],
    ops = [
        parentsTagger(),
        dumper(structure=False, max=6, at=[-1])
    ],
    gen = 2
)
#end_file

#begin_file log/generator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from random import randint
def randomChooser(pop, sp):
    males = []
    females = []
    # identify males and females in each social rank
    for rank in range(3):
        males.append([x for x in pop.individuals(sp) \
            if x.sex() == Male and x.info('rank') == rank])
        females.append([x for x in pop.individuals(sp) \
            if x.sex() == Female and x.info('rank') == rank])
    #
    while True:
        # choose a rank randomly
        rank = pop.individual(randint(0, pop.subPopSize(sp) - 1), sp).intInfo('rank')
        yield males[rank][randint(0, len(males[rank]) - 1)], \
            females[rank][randint(0, len(females[rank]) - 1)]

def setRank(pop, dad, mom, off):
    'The rank of offspring can increase or drop to zero randomly'
    off.setInfo((dad.info('rank') + randint(-1, 1)) % 3, 'rank')

pop = population(size=[1000, 2000], loci=1, infoFields='rank')
simu = simulator(pop, homoMating(
    pyParentsChooser(randomChooser),
    offspringGenerator(ops=mendelianGenoTransmitter())))
simu.evolve(
    preOps = [
        initSex(),
        initInfo(lambda : randint(0, 2), infoFields='rank')
    ],
    ops = [],
    gen = 5
)    
#end_file


#begin_file log/cppParentChooser.py
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore

# The class myParentsChooser is defined in module myParentsChooser
from myParentsChooser import myParentsChooser
def parentsChooser(pop, sp):
    'How to call a C++ level parents chooser.'
    # create an object with needed information (such as x, y) ...
    pc = myParentsChooser(
        [x for x in range(pop.popSize()) if pop.individual(x).sex() == Male],
        [x for x in range(pop.popSize()) if pop.individual(x).sex() == Female])
    while True:
        # return indexes of parents repeatedly
        yield pc.chooseParents()

pop = population(100, loci=1)
simu = simulator(pop,
    homoMating(pyParentsChooser(parentsChooser),
    offspringGenerator(ops=mendelianGenoTransmitter()))
)
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [],
    gen = 100
)
#end_file

#begin_file log/simuGen.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(50, loci=[10], ploidy=1),
    randomSelection(), rep=3)
simu.evolve(ops = [], gen = 5)
simu.gen()
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        stat(alleleFreq=5),
        ifElse('alleleNum[5][0] == 0',
            pyEval(r"'Allele 0 is lost in rep %d at gen %d\n' % (rep, gen)")),
        ifElse('alleleNum[5][0] == 50',
            pyEval(r"'Allele 0 is fixed in rep %d at gen %d\n' % (rep, gen)")),
        terminateIf('len(alleleNum[5]) == 1'),
    ],
)
simu.gen()
#end_file

#begin_file log/twoStage.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
# First stage: use the standard random mating scheme, do not use any
# information field for efficiency considerations.
simu = simulator(population(500, loci=[10]), randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [],
    gen = 50
)
# Second stage: track parents and produce more offspring per mating
# event. In preparation for pedigree ascertainment.
for pop in simu.populations():
    pop.addInfoFields(['father_idx', 'mother_idx'])
    pop.setAncestralDepth(1)

simu.setMatingScheme(randomMating(numOffspring=2))
simu.evolve(
    ops = [
        parentsTagger(),
        maPenetrance(loci=0, penetrance=(0.2, 0.4, 0.5))
    ],
    gen = 5
)
# Sample affected sibpairs
pop = simu.extract(0)
sample = AffectedSibpairSample(pop, size=5)[0]
[ind.intInfo('father_idx') for ind in sample.individuals()]
#end_file

#begin_file log/changeStru.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
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

# The populations start with no loci at all.
simu = simulator(population(1000, loci=[]), randomMating(), rep=3)
simu.evolve(
    preOps = initSex(),
    ops = [pyOperator(func=mutator, param=(10000, 2e-6))],
    gen = 200
)
for pop in simu.populations():
    print pop.totNumLoci(), pop.lociPos()
#end_file

#begin_file log/simuFunc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=[5, 10], infoFields='x'),
    randomMating(), rep=5)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.4, 0.6])
    ],
    ops=[],
    gen=10
)
# clone
cloned = simu.clone()
# save and load, using a different mating scheme
simu.save("sample.sim")
loaded = LoadSimulator("sample.sim", randomMating(numOffspring=2))
# 
simu.numRep()
loaded.numRep()
for pop1,pop2 in zip(cloned.populations(), loaded.populations()):
    assert pop1 == pop2

# continue to evolve
simu.evolve(ops=[], gen=10)
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[1000, 1000])
InitSex(pop, maleFreq=0.3, subPops=0)
InitSex(pop, sex=[Male, Female, Female], subPops=1)
Stat(pop, numOfMale=True, vars='numOfMale_sp')
print pop.dvars(0).numOfMale
print pop.dvars(1).numOfMale
#end_file

#begin_file log/initByFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2, 3], loci=[5, 7])
InitByFreq(pop, alleleFreq=[[.2, .8], [.8, .2]])
Dump(pop, structure=False)
#end_file

#begin_file log/initByFreqIdenticalInds.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2, 3], loci=[5, 7])
InitByFreq(pop, alleleFreq=[.2, .8], identicalInds=True)
Dump(pop, structure=False)
#end_file

#begin_file log/initByValue.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2, 3], loci=[5, 7])
InitByValue(pop, [1]*5 + [2]*7 + [3]*5 +[4]*7)
Dump(pop, structure=False)
#end_file

#begin_file log/initByValueProp.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[6, 8], loci=[5, 7])
pop.setVirtualSplitter(sexSplitter())
# initialize sex and the first two loci
InitSex(pop)
InitByValue(pop, loci=range(5), value=range(10))
# initialize all males
InitByValue(pop, loci=range(5, 12), value=[2]*7,
    subPops=[(0, 0), (1, 0)])
# initialize females by proportion
InitByValue(pop, loci=range(5, 12), ploidy=1, value=[[3]*7, [4]*7],
    subPops=[(0, 1), (1, 1)], proportions=[0.4, 0.6])
Dump(pop, structure=False)
#end_file

#begin_file log/initInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population(size=[5], loci=[2], infoFields=['sex', 'age'])
pop.setVirtualSplitter(sexSplitter())
InitSex(pop)
InitInfo(pop, 0, subPops=[(0,0)], infoFields='sex')
InitInfo(pop, 1, subPops=[(0,1)], infoFields='sex')
InitInfo(pop, lambda: random.randint(20, 70), infoFields='age')
Dump(pop, structure=False)
#end_file

#begin_file log/dumper.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[10, 10], loci=[20, 30], infoFields='gen',
    ancGen=-1)
pop.setVirtualSplitter(sexSplitter())
pop1 = pop.clone()
InitByFreq(pop, [0]*20 + [0.1]*10)
pop.setIndInfo(1, 'gen')
InitByFreq(pop1, [0]*50 + [0.1]*10)
pop1.setIndInfo(2, 'gen')
pop.push(pop1)
Dump(pop, width=3, loci=[5, 6, 30], subPops=([0, 0], [1, 1]),
    max=10, structure=False, ancGen=-1)
#end_file

#begin_file log/savePopulation.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=2),
    randomMating(), rep=5)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.2, 0.8])
    ],
    ops = [
        savePopulation(output="!'snapshot_%d_%d.pop' % (rep, gen)",
            step = 10),
        ],
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, infoFields=['father_idx', 'mother_idx']),
    randomMating(), rep=5)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.3, 0.7])
    ],
    ops = [
        setAncestralDepth(2, at=-2),
        parentsTagger(begin=-2)
    ],
    gen = 100
)
pop = simu.population(3)
print pop.ancestralGens()
print pop.ancestor(10, 1).info('father_idx')
#end_file

#begin_file log/ifElse.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=1000, loci=1),
    randomMating(), rep=4)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5]),
        pyExec('below40, above60 = 0, 0')
    ],
    ops = [
        stat(alleleFreq=0),
        ifElse('alleleFreq[0][1] < 0.4',
            pyExec('below40 += 1')),
        ifElse('alleleFreq[0][1] > 0.6',
            pyExec('above60 += 1')),
        ifElse('len(alleleFreq[0]) == 1',
            pyExec('stoppedAt = gen')),
        terminateIf('len(alleleFreq[0]) == 1')
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=100, loci=1),
    randomMating(), rep=10)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5]),
    ],
    ops = [
        stat(alleleFreq=0),
        terminateIf('len(alleleFreq[0]) == 1', stopAll=True)
    ]
)

#end_file

#begin_file log/debug.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=1), randomMating(), rep=5)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.1, 0.9])
    ],
    ops = [
        stat(alleleFreq=0),
        ifElse('alleleNum[0][0] == 0',
            ifOps = [
                turnOnDebug("DBG_MUTATOR"),
                pointMutator(loci=0, allele=0, inds=0),
            ],
            elseOps = turnOffDebug("DBG_MUTATOR")),
    ],
    gen = 100
)
#begin_ignore
TurnOffDebug("DBG_MUTATOR")
#end_ignore
#end_file

#begin_file log/pause.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100), randomMating(), rep=10)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [pause(stopOnKeyStroke=str(x), reps=x) for x in range(10)],
    gen = 100
)
#end_file

#begin_file log/ticToc.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(10000, loci=[100]*5), randomMating(), rep=2)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.1, 0.9])
    ],
    ops = [
        stat(alleleFreq=0),
        ticToc(step=50, reps=-1),
    ],
    gen = 101
)
#end_file

#begin_file log/pyExec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(100, loci=1),
    randomMating(), rep=2)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.2, 0.8]),
        pyExec('traj=[]')
    ],
    ops = [
        stat(alleleFreq=0),
        pyExec('traj.append(alleleFreq[0][1])'),
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(1000, loci=1,
    infoFields=['mother_idx', 'father_idx']),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        stat(alleleFreq=0),
        parentsTagger(),
        pyEval(r'"gen %d, #father %d, #mother %d\n"' \
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population(20, loci=1, infoFields='a')
pop.setVirtualSplitter(infoSplitter('a', cutoff=[3]))
InitByFreq(pop, [0.2, 0.8])
pop.setIndInfo([random.randint(2, 5) for x in range(20)], 'a')
InfoEval(pop, 'a', subPops=[(0, 0)]);print
InfoEval(pop, 'ind.allele(0, 0)', exposeInd='ind');print
# use population variables
pop.dvars().b = 5
InfoEval(pop, '"%d " % (a+b)', usePopVars=True);print
#end_file

#begin_file log/infoExec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(100, loci=1, infoFields=['a', 'b', 'c'])
InitByFreq(pop, [0.2, 0.8])
InfoExec(pop, 'a=1')
print pop.indInfo('a')[:10]
InfoExec(pop, 'b=ind.sex()', exposeInd='ind')
print pop.indInfo('b')[:10]
InfoExec(pop, 'c=a+b')
print pop.indInfo('c')[:10]
pop.dvars().d = 5
InfoExec(pop, 'c+=d', usePopVars=True)
print pop.indInfo('c')[:10]
#end_file

#begin_file log/migrateByProb.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=[1000]*3, infoFields='migrate_to'),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[
            [0, 0.1, 0.1],
            [0, 0, 0.1],
            [0, 0.1, 0]
        ]),
        stat(popSize=True),
        pyEval('subPopSize'),
        pyOutput('\n')
    ],
    gen = 5
)        
#end_file

#begin_file log/migrateByPropAndCount.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=[1000]*3, infoFields='migrate_to'),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[[0.1], [0.2]],
            mode=ByProportion,
            subPops=[1, 2],
            toSubPops=[3]),
        stat(popSize=True),
        pyEval('subPopSize'),
        pyOutput('\n')
    ],
    gen = 5
)        
#
simu.evolve(
    ops = [
        migrator(rate=[[50, 50], [100, 50]],
            mode=ByCounts,
            subPops=[3, 2],
            toSubPops=[2, 1]),
        stat(popSize=True),
        pyEval('subPopSize'),
        pyOutput('\n')
    ],
    gen = 5
)        
#end_file

#begin_file log/migrateVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[1000]*2, infoFields='migrate_to')
pop.setVirtualSplitter(sexSplitter())
simu = simulator(pop, randomMating())
simu.evolve(
    # 500 males and 500 females
    preOps = initSex(sex=[Male, Female]),
    ops = [
        migrator(rate=[
            [0, 0.10],
            [0, 0.05],
            ],
            mode = ByProportion,
            subPops=[(0, 0), (0, 1)]),
        stat(popSize=True, numOfMale=True, stage=PrePostMating, vars='numOfMale_sp'),
        pyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMale'], subPopSize[0], "
            "subPop[1]['numOfMale'], subPopSize[1])", stage=PrePostMating),
    ],
    gen = 2
)   
#end_file

#begin_file log/manualMigration.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population([10]*2, infoFields='migrate_to')
pop.setIndInfo([0, 1, 2, 3]*5, 'migrate_to')
Migrate(pop, mode=ByIndInfo)
pop.subPopSizes()
#end_file

#begin_file log/splitBySize.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(1000), randomSelection())
simu.evolve(
    ops = [
        splitSubPops(subPops=0, sizes=[300, 300, 400], at=2),
        stat(popSize=True),
        pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 4
)
#end_file

#begin_file log/splitByProp.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
def demo(gen, oldSize=[]):
    if gen < 2:
        return 1000 + 100 * gen
    else:
        return [x + 50 * gen for x in oldSize]

simu = simulator(population(1000),
    randomSelection(subPopSize=demo))
simu.evolve(
    ops = [
        splitSubPops(subPops=0, proportions=[.5]*2, at=2),
        stat(popSize=True),
        pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 4
)
#end_file

#begin_file log/splitByInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population([1000]*3, subPopNames=['a', 'b', 'c'], infoFields='x')
pop.setIndInfo([random.randint(0, 3) for x in range(1000)], 'x')
print pop.subPopSizes()
print pop.subPopNames()
SplitSubPops(pop, subPops=[0, 2], infoFields=['x'])
print pop.subPopSizes()
print pop.subPopNames()
#end_file

#begin_file log/mergeSubPops.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population([500]*2),
    randomSelection())
simu.evolve(
    ops = [
        mergeSubPops(subPops=[0, 1], at=3),
        stat(popSize=True),
        pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 5
)
#end_file

#begin_file log/resizeSubPops.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population([500]*2),
    randomSelection())
simu.evolve(
    ops = [
        resizeSubPops(proportions=(1.5, 2), at=3),
        stat(popSize=True),
        pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 5
)
#end_file

#begin_file log/recRate.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(size=[1000], loci=[100]),
    randomMating(ops = [
        recombinator(rates=0.01, reps=0),
        recombinator(rates=[0.01]*10, loci=range(50, 60), reps=1),
    ]), rep=2)
simu.evolve(
    preOps = [
        initSex(),
        initByValue([0]*100 + [1]*100)
    ],
    ops = [
        stat(LD=[[40, 55], [60, 70]]),
        pyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[60][70])'),
        pyOutput('\n', reps=-1)
    ],
    gen = 5
)
#end_file

#begin_file log/recIntensity.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(size=[1000], loci=3, lociPos=[0, 1, 1.1]),
    randomMating(ops=recombinator(intensity=0.01)))
simu.evolve(
    preOps = [
        initSex(),
        initByValue([0]*3 + [1]*3)
    ],
    ops = [
        stat(LD=[[0, 1], [1, 2]]),
        pyEval(r'"%.3f\t%.3f\n" % (LD_prime[0][1], LD_prime[1][2])', step=10)
    ],
    gen = 50
)
#end_file

#begin_file log/conversion.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(size=[1000], loci=[100]),
    randomMating(ops=[
        recombinator(rates=0.01, loci=50, reps=0),
        recombinator(rates=0.01, loci=50, reps=1, convMode=(NumMarkers, 1, 10)),
    ]), rep=2)
simu.evolve(
    preOps = [
        initSex(),
        initByValue([0]*100 + [1]*100)
    ],
    ops = [
        stat(LD=[[40, 55], [40, 70]]),
        pyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[40][70])'),
        pyOutput('\n', reps=-1)
    ],
    gen = 5
)
#end_file

#begin_file log/trackRec.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(1000, loci=[1000, 2000], infoFields='ind_id')
simu = simulator(pop, randomMating(ops = [
    idTagger(),
    recombinator(rates=0.001, output='>>rec.log', infoFields='ind_id')])
)
simu.evolve(
    preOps = [
        initSex(),
        idTagger(),
    ],
    ops = [],
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(size=[2000], loci=1),
    randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.2, 0.3, 0.5])
    ],
    ops = [
        matrixMutator(rate = [
            [0, 1e-5, 1e-5],
            [1e-4, 0, 1e-4],
            [1e-3, 1e-3, 0]
        ]),
        stat(alleleFreq=0, step=100),
        pyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(3)]) + '\n'",
            step=100),
    ],
    gen=1000
)
#end_file

#begin_file log/kamMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2000], loci=1*3)
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        kamMutator(k=5, rates=[1e-2, 1e-3], loci=[0, 1]),
        stat(alleleFreq=range(3), step=100),
        pyEval(r"', '.join(['%.3f' % alleleFreq[x][0] for x in range(3)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/snpMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2000], loci=[1, 1], infoFields='fitness')
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        snpMutator(u=0.001),
        maSelector(loci=0, fitness=[1, 0.99, 0.98]),
        stat(alleleFreq=[0, 1], step=100),
        pyEval(r"'%.3f\t%.3f\n' % (alleleFreq[0][1], alleleFreq[1][1])",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/acgtMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2000], loci=1,
    alleleNames=['A', 'C', 'G', 'T'])
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([.1, .1, .1, .7])
    ],
    ops = [
        acgtMutator(rate=[1e-4, 0.5], model='K80'),
        stat(alleleFreq=0, step=100),
        pyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(4)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end_file

#begin_file log/smmMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(size=1000, loci=[1, 1]), randomMating())
simu.evolve(
    # all start from allele 50
    preOps = [
        initSex(),
        initByFreq( [0]*50 + [1])
    ],
    ops = [
        smmMutator(rates=1e-3, loci=0),
        smmMutator(rates=1e-3, incProb=0.6, loci=1,
            mutStep=(GeometricDistribution, 0.2)),
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
def incAllele(allele):
    return allele + random.randint(1, 5)

simu = simulator(population(size=1000, loci=[20]),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        pyMutator(func=incAllele, rates=[1e-4, 1e-3],
            loci=[2, 10])
    ],
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(5000, loci=[1, 1]),
    randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByValue([50, 50])
    ],
    ops = [
        # the first locus uses a pure stepwise mutation model
        smmMutator(rates=0.001, loci=0),
        # the second locus uses a mixed model
        mixedMutator(rates=0.001, loci=1, mutators=[        
            kamMutator(rates=1, k=100),
            smmMutator(rates=1)
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population(5000, loci=[3, 3]),
    randomMating())
simu.evolve(
    # initialize locus by 0, 0, 0, 1, 0, 1
    preOps = [
        initSex(),
        initByValue([1, 1], loci=[3, 5])
    ],
    ops = [
        contextMutator(mutators=[
            snpMutator(u=0.1),
            snpMutator(u=1),
            ],
            contexts=[(0, 0), (1, 1)],
            loci=[1, 4],
            rates=0.01
        ),
        stat(alleleFreq=[1, 4], step=5),
        pyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    gen = 20
)
#end_file

#begin_file log/pyContextMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
simu = simulator(population(5000, loci=[3, 3]),
    randomMating())
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
    preOps = [
        initSex(),
        initByValue([1, 1], loci=[3, 5])
    ],
    ops = [
        pyMutator(func=contextMut, context=1,
            loci=[1, 4],  rates=0.01
        ),
        #snpMutator(u=0.01, v= 0.01, loci=[1, 4]),
        stat(alleleFreq=[1, 4], step=5),
        pyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    gen = 20
)
#end_file

#begin_file log/pointMutator.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(1000, loci=1, infoFields='fitness')
simu = simulator(pop, randomSelection())
simu.evolve(
    preOps = pyOutput('Introducing alleles at generation'),
    ops = [
        stat(alleleFreq=0),
        ifElse('alleleNum[0][1] == 0', ifOps=[
            pyEval(r"' %d' % gen"),
            pointMutator(inds=0, loci=0, allele=1),
        ]),
        maSelector(loci=0, wildtype=0, fitness=[1, 1.05, 1.1]),
        ifElse('alleleFreq[0][1] > 0.05', ifOps=[
            pyEval(r"'.\nTerminate at generation %d at allele freq %.3f.\n'" +
                " % (gen, alleleFreq[0][1])"),
            terminateIf('True'),
        ])
    ],
)
#end_file

#begin_file log/mutatorVSP.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
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
    'Get average allele by affection status.'
    Stat(pop, alleleFreq=(0,1), subPops=[(0,0), (0,1)],
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

pop = population(10000, loci=[1, 1])
pop.setVirtualSplitter(affectionSplitter())
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByValue([50, 50])
    ],
    ops = [
        # determine affection status for each offspring (duringMating)
        pyPenetrance(func=fragileX, loci=0),
        # unaffected offspring, mutation rate is high to save some time
        smmMutator(rates=1e-3, loci=1),
        # unaffected offspring, mutation rate is high to save some time
        smmMutator(rates=1e-3, loci=0, subPops=[(0, 0)]),
        # affected offspring have high probability of mutating upward
        smmMutator(rates=1e-2, loci=0, subPops=[(0, 1)],
           incProb=0.7, mutStep=3),
        # number of affected
        pyOperator(func=avgAllele, step=20),
        pyEval(r"'Gen: %3d #Aff: %d AvgRepeat: %.2f (unaff), %.2f (aff), %.2f (unrelated)\n'"
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2000], loci=1)
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0]*4 + [0.1, 0.2, 0.3, 0.4])
    ],
    ops = [
        kamMutator(k=4, rates=1e-4, mapIn=[0]*4 + range(4),
            mapOut=[4, 5, 6, 7]),
        stat(alleleFreq=0, step=100),
        pyEval(r"', '.join(['%.2f' % alleleFreq[0][x] for x in range(8)]) + '\n'",
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
from simuPOP import *
#begin_ignore
GetRNG().setSeed(12345)
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
            loc = GetRNG().randGeometric(rate)
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

pop = population(size=[2000], loci=[100])
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
    ],
    ops = [
        # recombine in a 10Mb region at rate 1e-8
        pyOperator(func=infSitesMutate, param=(1, 10000000, 1e-8)),
    ],
    gen = 100
)
# now, we get a population. Let us have a look at the 'alleles'.
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population([5000]*3, loci=5), randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [
        stat(structure=range(5), subPops=(0, 1), suffix='_01', step=40),
        stat(structure=range(5), subPops=(1, 2), suffix='_12', step=40),
        stat(structure=range(5), subPops=(0, 2), suffix='_02', step=40),
        stat(structure=range(5), step=40),
        pyEval(r"'Fst=%.3f (pairwise: %.3f %.3f %.3f)\n' % (F_st, F_st_01, F_st_12, F_st_02)",
            step=40),
    ],
    gen = 200
)
#end_file

#begin_file log/statCount.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(10000, loci=1)
pop.setVirtualSplitter(combinedSplitter(
    [sexSplitter(), affectionSplitter()]))
InitSex(pop)
InitByFreq(pop, [0.2, 0.8])
MaPenetrance(pop, loci=0, penetrance=[0.1, 0.2, 0.5])
# Count population size
Stat(pop, popSize=True, subPops=[(0, 0), (0, 2)])
# popSize is the size of two VSPs, does not equal to total population size.
# Because two VSPs overlap (all males and all unaffected), popSize can be
# greater than real population size.
print pop.dvars().subPopSize, pop.dvars().popSize
# print popSize of each virtual subpopulation.
Stat(pop, popSize=True, subPops=[(0, 0), (0, 2)], vars='popSize_sp')
# Note the two ways to access variable in (virtual) subpopulations.
print pop.dvars((0,0)).popSize, pop.dvars().subPop[(0,2)]['popSize']
# Count number of male (should be the same as the size of VSP (0,0).
Stat(pop, numOfMale=True)
print pop.dvars().numOfMale
# Count the number of affected and unaffected male individual
Stat(pop, numOfMale=True, subPops=[(0, 2), (0, 3)], vars='numOfMale_sp')
print pop.dvars((0,2)).numOfMale, pop.dvars((0,3)).numOfMale
# or number of affected male and females
Stat(pop, numOfAffected=True, subPops=[(0, 0), (0, 1)], vars='numOfAffected_sp')
print pop.dvars((0,0)).numOfAffected, pop.dvars((0,1)).numOfAffected
# These can also be done using a productSplitter...
pop.setVirtualSplitter(productSplitter(
    [sexSplitter(), affectionSplitter()]))
Stat(pop, popSize=True, subPops=[(0, x) for x in range(4)])
# counts for male unaffected, male affected, female unaffected and female affected
print pop.dvars().subPopSize
#end_file

#begin_file log/statAlleleFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(10000, loci=1)
pop.setVirtualSplitter(affectionSplitter())
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq(loci=0, alleleFreq=[0.8, 0.2])
    ],
    ops = [
        maPenetrance(penetrance=[0.1, 0.4, 0.6], loci=0),
        stat(alleleFreq=0, subPops=[(0, 0), (0, 1)],
            vars=['alleleFreq', 'alleleFreq_sp']),
        pyEval(r"'Gen: %d, freq: %.2f, freq (aff): %.2f, freq (unaff): %.2f\n' % " + \
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(100, loci=[1, 1, 1], chromTypes=[Autosome, ChromosomeX, ChromosomeY])
InitByFreq(pop, [0.01, 0.05, 0.94])
Stat(pop, genoFreq=[0, 1])
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(100, loci=1)
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [
        stat(heteroFreq=0, step=10),
        pyEval(r"'Gen: %d, HeteroFreq: %.2f\n' % (gen, heteroFreq[0])", step=20)
    ],
    gen = 100
)
#end_file

#begin_file log/statHaploFreq.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(100, loci=3)
InitByFreq(pop, [0.2, 0.4, 0.4], loci=0)
InitByFreq(pop, [0.2, 0.8], loci=2)
Stat(pop, genoFreq=[0, 1, 2], haploFreq=[0, 1, 2],
    vars=['genoNum', 'haploFreq'])
utils.ViewVars(pop.vars(), gui=False)
#end_file

#begin_file log/statInfo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
pop = population([500], infoFields='anc')
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(infoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
#
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        # anc is 0 or 1
        initInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    ops = [
        inheritTagger(mode=Mean, infoFields='anc'),
        stat(popSize=True, meanOfInfo='anc', varOfInfo='anc',
            subPops=[(0,x) for x in range(5)]),
        pyEval(r"'Anc: %.2f (%.2f), #inds: %s\n' %" + \
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population([1000]*2, loci=3)
InitByFreq(pop, [0.2, 0.8], subPops=0)
InitByFreq(pop, [0.8, 0.2], subPops=1)
Stat(pop, LD=[[0, 1, 0, 0], [1, 2]],
    vars=['LD', 'LD_prime', 'R2', 'LD_ChiSq', 'LD_ChiSq_p', 'CramerV',
        'LD_prime_sp', 'LD_ChiSq_p_sp'])
from pprint import pprint
pprint(pop.vars())
#end_file

#begin_file log/statAssociation.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import *
def assoTest(pop):
    'Draw case-control sample and apply association tests'
    sample = CaseControlSample(pop, cases=500, controls=500)[0]
    Stat(sample, association=(0, 2), vars=['Allele_ChiSq_p', 'Geno_ChiSq_p', 'Armitage_p'])
    print 'Allele test: %.2e, %.2e, Geno test: %.2e, %.2e, Trend test: %.2e, %.2e' \
        % (sample.dvars().Allele_ChiSq_p[0], sample.dvars().Allele_ChiSq_p[2],
        sample.dvars().Geno_ChiSq_p[0], sample.dvars().Geno_ChiSq_p[2],
        sample.dvars().Armitage_p[0], sample.dvars().Armitage_p[2])
    return True

simu = simulator(population(size=100000, loci=3),
    randomMating(ops=recombinator(loci=[0, 1], rates=[0.01, 0.005])))
simu.evolve(
    preOps = [
        initSex(),
        initByValue([[0]*3, [1]*3], proportions=[0.5, 0.5])
    ],
    ops = [
        maPenetrance(loci=1, penetrance=[0.1, 0.2, 0.4]),
        pyOperator(func=assoTest, step=20),
    ],
    gen = 100
)
#end_file

#begin_file log/statStructure.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import MigrIslandRates
simu = simulator(population([5000]*3, loci=10, infoFields='migrate_to'),
    randomMating(), rep=2)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5])
    ],
    ops = [
        migrator(rate=MigrIslandRates(0.01, 3), reps=1),
        stat(structure=range(10), step=40),
        pyEval("'Fst=%.3f (rep=%d without migration) ' % (F_st, rep)", step=40, reps=0),
        pyEval("'Fst=%.3f (rep=%d with migration) ' % (F_st, rep)", step=40, reps=1),
        pyOutput('\n', reps=-1, step=40)
    ],
    gen = 200
)
#end_file

#begin_file log/statHWE.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(population([1000], loci=1), randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByValue([[0,0], [0, 1], [1,1]], proportions=[0.4, 0.4, 0.2])
    ],
    ops = [
        stat(HWE=0, genoFreq=0, stage=PrePostMating),
        pyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
            'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])',
            stage=PrePostMating),
    ],
    gen = 1
)
#end_file

#begin_file log/inheritTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[1000]*10, loci=1, infoFields='x')
# tag the first individual of each subpopulation.
for sp in range(pop.numSubPop()):
    pop.individual(0, sp).setInfo(1, 'x')

simu = simulator(pop, randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        inheritTagger(mode=Maximum, infoFields='x'),
        stat(sumOfInfo='x', vars=['sumOfInfo_sp']),
        pyEval(r'", ".join(["%3d" % subPop[i]["sumOfInfo"]["x"] for i in range(10)])+"\n"'),
    ],
    gen = 5
)
#end_file

#begin_file log/summaryTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
idTagger().reset(0)
#end_ignore
simu = simulator(
    population(1000, loci=1, infoFields=['fitness', 'avgFitness']),
    randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5]),
    ],
    ops = [
        maSelector(loci=0, wildtype=0, fitness=[1, 0.99, 0.95]),
        summaryTagger(mode=Mean, infoFields=['fitness', 'avgFitness']),
        stat(alleleFreq=0, meanOfInfo='avgFitness', step=10),
        pyEval(r"'gen %d: allele freq: %.3f, average fitness of parents: %.3f\n' % "
            "(gen, alleleFreq[0][1], meanOfInfo['avgFitness'])", step=10)
    ],
    gen = 50,
)
#end_file


#begin_file log/idTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
idTagger().reset(0)
#end_ignore
simu = simulator(
    population(10, infoFields='ind_id', ancGen=1),
    randomSelection())
simu.evolve(
    preOps = idTagger(),
    ops = idTagger(),
    gen = 1
)
pop = simu.extract(0)
print [ind.intInfo('ind_id') for ind in pop.individuals()]
pop.useAncestralGen(1)
print [ind.intInfo('ind_id') for ind in pop.individuals()]
TagID(pop) # re-assign ID
print [ind.intInfo('ind_id') for ind in pop.individuals()]
#end_file

#begin_file log/pedigreeTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
idTagger().reset(0)
#end_ignore
simu = simulator(
    population(100, infoFields=['ind_id', 'father_id', 'mother_id']),
    randomMating())
simu.evolve(
    preOps = [
        initSex(),
        idTagger()
    ],
    ops = [
        idTagger(),
        pedigreeTagger(output=">>pedigree.txt")
    ],
    gen = 100
)
ped = open('pedigree.txt')
print ''.join(ped.readlines()[100:105])
#begin_ignore
ped.close()
import os
os.remove('pedigree.txt')
#end_ignore
#end_file

#begin_file log/pyTagger.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import random
def randomMove(values):
    '''Pass parental information fields to offspring'''
    x1, y1, x2, y2 = values
    # shift right with high concentration of alleles... 
    x = random.normalvariate((x1+x2)/2., 0.1)
    y = random.normalvariate((y1+y2)/2., 0.1)
    return (x, y)

pop = population(1000, loci=[1], infoFields=['x', 'y'])
pop.setVirtualSplitter(genotypeSplitter(loci=0, alleles=[[0, 0], [0,1], [1, 1]]))
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5]),
        initInfo(random.random, infoFields=['x', 'y'])
    ],
    ops = [
        pyTagger(func=randomMove, infoFields=['x', 'y']),
        stat(minOfInfo='x', maxOfInfo='x'),
        pyEval(r"'Range of x: %.2f, %.2f\n' % (minOfInfo['x'], maxOfInfo['x'])")
    ],
    gen = 5
)

#end_file

#begin_file log/otherTagging.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(1000, loci=[1], infoFields=['aff', 'numAff'])
# define virtual subpopulations by affection status
pop.setVirtualSplitter(affectionSplitter())
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.5, 0.5]),
    ],
    ops = [
        # get affection status for both parents and offspring
        maPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4], stage=PrePostMating),
        # set 'aff' of parents
        infoExec('aff = ind.affected()', exposeInd='ind', stage=PreMating),
        # get number of affected parents for each offspring and store in numAff
        summaryTagger(mode=Summation, infoFields=['aff', 'numAff']),
        # calculate mean 'numAff' of offspring, for unaffected and affected subpopulations.
        stat(meanOfInfo='numAff', subPops=[(0,0), (0,1)], vars=['meanOfInfo_sp']),
        # print mean number of affected parents for unaffected and affected offspring.
        pyEval(r"'Mean number of affected parents: %.2f (unaff), %.2f (aff)\n' % "
            "(subPop[(0,0)]['meanOfInfo']['numAff'], subPop[(0,1)]['meanOfInfo']['numAff'])")
    ],
    gen = 5
)

#end_file


#begin_file log/forwardTrajectory.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import trajectory, ForwardTrajectory

traj = ForwardTrajectory(N=[2000, 4000], fitness=[1, 0.99, 0.98],
    beginGen=0, endGen=100, beginFreq=[0.2, 0.3],
    endFreq=[[0.1, 0.11], [0.2, 0.21]])
traj.plot('forwardTrajectory.png', plot_ylim=[0, 0.5], col_sp=['red', 'blue'],
    plot_main='Simulated trajectory (forward-time)')
simu = simulator(
    population(size=[2000, 4000], loci=10, infoFields='fitness'),
    controlledRandomMating(
        ops=[recombinator(rates=0.01)],
        loci=5, alleles=1, freqFunc=traj.func())
)
simu.evolve(
    preOps = [
        initSex(),
        initByFreq([0.8, 0.2], subPops=0),
        initByFreq([0.7, 0.3], subPops=1),
        pyOutput('Sp0: loc2\tloc5\tSp1: loc2\tloc5\n'),
    ],
    ops = [
        stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=20),
        pyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP.utils import trajectory, BackwardTrajectory
from math import exp
def Nt(gen, oldSize=[]):
    'An exponential population growth demographic model.'
    return int((10**4) * exp(.00115 * gen))

def fitness(gen, sp):
    'Constant positive selection pressure.'
    return [1, 1.01, 1.02]

# simulate a trajectory backward in time, from generation 1000
traj = BackwardTrajectory(N=Nt, fitness=fitness, nLoci=2,
     endGen=1000, endFreq=[0.1, 0.2])
traj.plot('backTrajectory.png', plot_ylim=[0, 0.3], plot_xlim=[0, 1000],
    col_loc=['red', 'blue'], plot_main='Simulated trajectory (backward-time)')
print 'Trajectory simulated with length %s ' % len(traj.traj)
pop = population(size=Nt(0), loci=[1]*2)
# save trajectory function in the population's local namespace
# so that the pyEval operator can access it.
pop.dvars().traj = traj.func()
simu = simulator(pop, controlledRandomMating(loci=[0, 1], alleles=[1, 1],
        subPopSize=Nt, freqFunc=traj.func()))
simu.evolve(
    preOps = [initSex()],
    ops = traj.mutators(loci=[0, 1]) + [
        stat(alleleFreq=[0, 1], begin=500, step=100),
        pyEval(r"'%4d: %.3f (exp: %.3f), %.3f (exp: %.3f)\n' % (gen, alleleFreq[0][1],"
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP import *
from simuPOP.plotter import varPlotter
pop = population(size=1000, loci=2)
simu = simulator(pop, randomMating(ops=recombinator(rates=0.01)), rep=3)
simu.evolve(
    preOps = [
        initSex(),
        initByValue([1, 2, 2, 1])
    ],
    ops = [
        stat(LD=[0, 1]),
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP import *
from simuPOP.plotter import varPlotter
pop = population(size=1000, loci=1*4)
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
    preOps = [initSex()] +
        [initByFreq([0.1*(x+1), 1-0.1*(x+1)], loci=x) for x in range(4)],
    ops = [
        stat(alleleFreq=range(4)),
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP import *
from simuPOP.plotter import varPlotter
pop = population(size=1000, loci=1*4)
simu = simulator(pop, randomMating(), rep=3)
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
    preOps = [initSex()]+
        [initByFreq([0.1*(x+1), 1-0.1*(x+1)], loci=x) for x in range(4)],
    ops = [
        stat(alleleFreq=range(4)),
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP import *
from simuPOP.plotter import scatterPlotter
import random
pop = population([500], infoFields=['x', 'y', 'anc'])
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(infoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
#
def passInfo(fields):
    'Parental fields will be passed as x1, y1, anc1, x2, y2, anc2'
    x1, y1, anc1, x2, y2, anc2 = fields
    anc = (anc1 + anc2)/2.
    x = (x1 + x2)/2 + random.normalvariate(anc - 0.5, 0.1)
    y = (y1 + y2)/2 + random.normalvariate(0, 0.1)
    return x, y, anc

simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        # random geographic location
        initInfo(random.random, infoFields=['x', 'y']),
        # anc is 0 or 1
        initInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    ops = [
        pyTagger(passInfo, infoFields=['x', 'y', 'anc']),
        scatterPlotter(['x', 'y'], stage=PreMating,            
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
from simuPOP import *
from simuPOP.plotter import histPlotter, qqPlotter, boxPlotter
import random
pop = population([500], infoFields=['x', 'y', 'anc'])
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sexSplitter())
def passInfo(fields):
    'Parental fields will be passed as x1, y1, anc1, x2, y2, anc2'
    x1, y1, anc1, x2, y2, anc2 = fields
    anc = (anc1 + anc2)/2.
    x = (x1 + x2)/2 + random.normalvariate(anc - 0.5, 0.1)
    y = (y1 + y2)/2 + random.normalvariate(0, 0.1)
    return x, y, anc

simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [
        initSex(),
        # random geographic location
        initInfo(random.random, infoFields=['x', 'y']),
        # anc is 0 or 1
        initInfo(lambda : random.randint(0, 1), infoFields='anc')
    ],
    ops = [
        pyTagger(passInfo, infoFields=['x', 'y', 'anc']),
        histPlotter(infoFields='anc', stage=PreMating,            
            subPops=[(0,0), (0,1)], col_sp=['blue', 'red'],
            saveAs='log/histPlotter.png',
            main="!'Histogram of ancestry values at generation %d' % gen",
        ),
        qqPlotter(infoFields='anc', stage=PreMating,            
            subPops=[(0,0), (0,1)], col_sp=['blue', 'red'],
            saveAs='log/qqPlotter.png',
            main="!'QQ plot of ancestry values at generation %d' % gen",
        ),
        boxPlotter(infoFields='anc', stage=PreMating,
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
from simuPOP import *
GetRNG().setSeed(12345)
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
     'label': 'Initial population',
     'allowedTypes': [types.StringType],
     'description': '''Use one of the HapMap populations as the initial
            population for this simulation. You can choose from:
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
if not os.path.isfile('getParam.png'):
    print 'Run a GUI if getParam has not been runned'
else:
    sys.argv = ['getParam.py', '--rate=[0.25]', '--rep=5', '--pop="CEU"']
    simuOpt.setOptions(gui=False)

from simuPOP import *
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

#begin_file log/reichDemo.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
import math
def demo_model(model, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    model: linear or exponential
    N0:   Initial population size.
    N1:   Ending population size.
    G0:   Length of burn-in stage.
    G1:   Length of population expansion stage.
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
# population size at generation 700
print demo_func(700)
#end_file

#begin_file log/reichStat.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
class ne(pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    #
    def calcNe(self, pop):
        Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(freq[x]/sumFreq)**2 for x in freq.keys() if x != 0])
        # save the result to the population.
        pop.dvars().ne = ne
        return True

def Ne(pop, loci):
    '''Function form of operator ne'''
    ne(loci).apply(pop)
    return pop.dvars().ne

pop = population(100, loci=[10])
InitByFreq(pop, [.2] * 5)
print Ne(pop, loci=[2, 4])
#end_file

#begin_file log/reichEvolve.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore

#begin_ignore
import math
def demo_model(model, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    model: linear or exponential
    N0:   Initial population size.
    N1:   Ending population size.
    G0:   Length of burn-in stage.
    G1:   Length of population expansion stage.
    '''
    rate = (math.log(N1) - math.log(N0))/G1
    def ins_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:
            return N1
    
    def exp_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    
    if model == 'instant':
        return ins_expansion
    elif model == 'exponential':
        return exp_expansion

class ne(pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    
    def calcNe(self, pop):
        Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(freq[x]/sumFreq)**2 for x in freq.keys() if x != 0])
        
        # save the result to the population.
        pop.dvars().ne = ne
        return True

#end_ignore

def simulate(model, N0, N1, G0, G1, spec, s, mu, k):
    '''Evolve a population using given demographic model
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
    simu = simulator(
        population(size=demo_func(0), loci=1, infoFields='fitness'),
        randomMating(subPopSize=demo_func)
    )
    simu.evolve(
        preOps = [
            initSex(),
            initByFreq(loci=0, alleleFreq=spec)
        ],
        ops=[
            kamMutator(k=k, rates=mu),
            maSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=[0], step=100),
            pyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
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
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
#!/usr/bin/env python
#
# Author:  Bo Peng
# Purpose: A real world example for simuPOP user's guide.
#
'''
Simulation the evolution of allelic spectra (number and frequencies
of alleles at a locus), under the influence of population expansion,
mutation, and natural selection.
'''
import simuOpt
simuOpt.setOptions(quiet=True, alleleType='long')
from simuPOP import *
import sys, types, os, math
options = [
    {'longarg': 'demo=',
     'default': 'instant',
     'label': 'Population growth model',
     'description': 'How does a population grow from N0 to N1.',
     'chooseOneOf': ['instant', 'exponential'],
    },
    {'longarg': 'N0=',
     'default': 10000,
     'label': 'Initial population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial population size. This size will be maintained
                till the end of burnin stage''',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'N1=',
     'default': 100000,
     'label': 'Final population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Ending population size (after population expansion)',
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
     'description': 'Number of geneartions of the population expansion stage',
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
     'label': 'Maximum allelic state',
     'allowedTypes': [types.IntType],
     'description': 'Maximum allelic state for a k-allele mutation model',
     'validate': simuOpt.valueGT(1),
    },
]

def demo_model(type, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    type: linear or exponential
    N0:   Initial population size.
    N1:   Ending population size.
    G0:   Length of burn-in stage.
    G1:   Length of population expansion stage.
    '''
    rate = (math.log(N1) - math.log(N0))/G1
    def ins_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:
            return N1
    
    def exp_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    
    if type == 'instant':
        return ins_expansion
    elif type == 'exponential':
        return exp_expansion

class ne(pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    
    def calcNe(self, pop):
        Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc][1:]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(x/sumFreq)**2 for x in freq])
        # save the result to the population.
        pop.dvars().ne = ne
        return True

def simuCDCV(model, N0, N1, G0, G1, spec, s, mu, k):
    '''Evolve a population using given demographic model
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
    simu = simulator(
        population(size=demo_func(0), loci=1, infoFields='fitness'),
        randomMating(subPopSize=demo_func)
    )
    simu.evolve(
        preOps = [
            initSex(),
            initByFreq(loci=0, alleleFreq=spec)
        ],
        ops = [
            kamMutator(rate=mu, maxAllele=k),
            maSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=0, step=100),
            pyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
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

#begin_file log/mapSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=1000, ploidy=2, loci=1, infoFields='fitness'),
    randomMating())
s1 = .1
s2 = .2
simu.evolve(
    preOps = [
        initSex(),
        initByFreq(alleleFreq=[.2, .8])
    ],
    ops = [
        stat(alleleFreq=0, genoFreq=0),
        mapSelector(loci=0, fitness={(0,0):(1-s1), (0,1):1, (1,1):(1-s2)}),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    gen=300
)
#end_file


#begin_file log/maSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=1000, ploidy=2, loci=1, infoFields='fitness'),
    randomMating())
s1 = .1
s2 = .2
simu.evolve(
    preOps = [
        initSex(),
        initByFreq(alleleFreq=[.2, .8])
    ],
    ops = [
        stat(alleleFreq=0, genoFreq=0),
        maSelector(loci=0, fitness=[1-s1, 1, 1-s2]),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    gen = 300)
#end_file


#begin_file log/mlSelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=10, ploidy=2, loci=2, 
    infoFields=['fitness', 'spare']),
    randomMating())
simu.evolve(
    preOps = [
        initSex(),
        initByFreq(alleleFreq=[.2, .8])
    ],
    ops = [ mlSelector([
         mapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}),
         mapSelector(loci=1, fitness={(0,0):1, (0,1):1, (1,1):.8}),
         ], mode = Additive),
    ],
    gen = 2
)
#end_file

#begin_file log/pySelector.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
simu = simulator(
    population(size=1000, ploidy=2, loci=3, infoFields='fitness'),
    randomMating()
)
s1 = .2
s2 = .3
# the second parameter gen can be used for varying selection pressure
def sel(arr, gen=0):
  if arr[0] == 1 and arr[1] == 1:
    return 1 - s1
  elif arr[0] == 1 and arr[1] == 2:
    return 1
  elif arr[0] == 2 and arr[1] == 1:
    return 1
  else:
    return 1 - s2

# test func
print sel([1, 1])
simu.evolve(
    preOps = [
        initSex(),
        initByFreq(alleleFreq=[.2, .8])
    ],
    ops = [
        stat(alleleFreq=0, genoFreq=0),
        pySelector(loci=[0, 1], func=sel),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=25)
    ],
    gen=100
)
#end_file

#begin_file log/mapPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(size=[2,8], ploidy=2, loci=2 )
InitByFreq(pop, [.2, .8])
MapPenetrance(pop, loci=0, 
    penetrance={(0,0):0, (0,1):1, (1,1):1})
Stat(pop, numOfAffected=1)
#end_file

#begin_file log/mlPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(1000, loci=3)
InitByFreq(pop, [0.3, 0.7])
pen = []
for loc in (0, 1, 2):
    pen.append(maPenetrance(loci=loc, wildtype=[1],
        penetrance=[0, 0.3, 0.6] ) )

# the multi-loci penetrance
MlPenetrance(pop, mode=Multiplicative, peneOps=pen)
Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#end_file

#begin_file log/pyPenetrance.py
#begin_ignore
import simuOpt
simuOpt.setOptions(quiet=True)
from simuPOP import *
GetRNG().setSeed(12345)
#end_ignore
pop = population(1000, loci=3)
InitByFreq(pop, [0.3, 0.7])
def peneFunc(geno):
    p = 1
    for l in range(len(geno)/2):
        p *= (geno[l*2]+geno[l*2+1])*0.3
    #
    return p

PyPenetrance(pop, func=peneFunc, loci=(0, 1, 2))
Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#
# You can also define a function, that returns a penetrance
# function using given parameters
def peneFunc(table):
    def func(geno):
      return table[geno[0]][geno[1]]
    #  
    return func

# then, given a table, you can do
PyPenetrance(pop, loci=(0, 1, 2),
    func=peneFunc( ((0, 0.5), (0.3, 0.8)) ) )
#end_file

