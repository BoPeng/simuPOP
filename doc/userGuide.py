# This script will genreate code/result pieces
# that will be inserted into simuPOP user's guide
# and reference manual
# 
# #file_begin file
#
# #file_end file
#
# will be used in this file so the running result can be
# separated into files specified by 'file'
#
# #PS after commands that will be executed.
#
# create directory log if not exist
import os, sys

if not os.path.isdir('log'):
  try:
    os.mkdir('log')
  except:
    print "Failed to make output directory log"
    sys.exit(1)

#

#file log/standard.log
from simuPOP import *
pop = population(10, loci=[2])
pop.locusPos(10)
pop.individual(20).setAllele(1, 0)
#end

# make sure each run generates the same output to avoid unnecessary
# documentation changes.
rng().setSeed(12345)

#file log/importSimuPOP.log
import simuOpt
simuOpt.setOptions(optimized=False, alleleType='long', quiet=False)
from simuPOP import *
#end

#file log/simpleExample.log
from simuPOP import *
pop = population(size=1000, loci=[2])
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
    preOps = [initByValue([1, 2, 2, 1])],  
    ops = [
        recombinator(rate=0.01),
        stat(LD=[0, 1]),
        pyEval(r"'%.2f\t' % LD[0][1]", step=10),
        pyOutput('\n', rep=-1, step=10)
    ],
    gen=100
)
#end

#file log/help.log
help(population.addInfoField)
#end

#file log/absIndex.log
pop = population(size=[10, 20], loci=[5, 7])
print pop.chromLocusPair(7)
print pop.absLocusIndex(1, 1)
print pop.absIndIndex(2, 1)
print pop.subPopIndPair(25)
#end

#file log/iterator.log
pop = population(size=2, loci=[5, 6])
InitByFreq(pop, [0.2, 0.3, 0.5])
for ind in pop.individuals():
    for loc in range(pop.chromBegin(1), pop.chromEnd(1)):
        print ind.allele(loc),
    print 

#end

#file log/carray.log
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
#end

#file log/genoStru.log
pop = population(size=[2, 3], ploidy=2, loci=[5, 10],
    lociPos=[range(0, 5), range(0, 20, 2)], chromNames=['Chr1', 'Chr2'],
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
#end


#file log/haplodiploid.log
pop = population(size=[2,5], ploidy=Haplodiploid, loci=[3, 5])
InitByFreq(pop, [0.3, 0.7])
Dump(pop)
#end

#file log/chromType.log
pop = population(size=6, ploidy=2, loci=[3, 3, 6, 4, 4, 4],
    chromTypes=[Autosome]*2 + [ChromosomeX, ChromosomeY] + [Customized]*2)
InitByFreq(pop, [0.3, 0.7])
Dump(pop, structure=False) # does not display genotypic structure information
#end

#file log/infoField.log
pop = population(10, loci=[20], ancGen=1,
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [initByValue([0]*20+[1]*20)],
    ops = [
        parentsTagger(),
        recombinator(rate=0.01)
    ],
    gen = 1
)
pop = simu.extract(0)
pop.indInfo('mother_idx')  # mother of all offspring
ind = pop.individual(0)
mom = pop.ancestor(ind.intInfo('mother_idx'), 1)
print ind.genotype(0)
print mom.genotype(0)
print mom.genotype(1)
#end


#file log/individual.log
pop = population([5, 4], loci=[2, 5], infoFields=['x'])
# get an individual
ind = pop.individual(3)
ind.ploidy()            # access to genotypic structure
ind.numChrom()
ind.affected()
ind.setAffected(True)   # access affection status,
ind.sex()               # sex,
ind.setInfo(4, 'x')     # and information fields
ind.info('x')
#end

#file log/individual_genotype.log
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
#end


#file log/subPop.log
pop = population(size=[3, 4, 5], ploidy=1, loci=[1], infoFields=['x'])
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
#end


#file log/subPopName.log
pop = population(size=[3, 4, 5], subPopNames=['x', 'y', 'z'])
pop.removeSubPops([1])
pop.subPopNames()
pop.subPopByName('z')
pop.splitSubPop(1, [2, 3])
pop.subPopNames()
pop.setSubPopName('z-1', 1)
pop.subPopNames()
pop.subPopByName('z')
#end

#file log/virtualSplitter.log
import random
pop = population(size=[200, 400], loci=[30], infoFields=['x'])
# assign random information fields
pop.setIndInfo([random.randint(0, 3) for x in range(pop.popSize())], 'x')
# define a virtual splitter by information field 'x'
pop.setVirtualSplitter(infoSplitter(field='x', values=[0, 1, 2, 3]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # Size of VSP 0 in subpopulation 0
pop.subPopSize([1, 0])    # Size of VSP 0 in subpopulation 1
# use a combined splitter that defines additional VSPs by sex
InitSex(pop)
pop.setSubPopName('subPop 1', 0)
pop.setVirtualSplitter(combinedSplitter([
    infoSplitter(field='x', values=[0, 1, 2, 3]),
    sexSplitter()])
)
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 4])    # VSP 4 is the first VSP defined by the sex splitter
pop.subPopSize([0, 4])    # Number of male individuals
#end


#file log/virtualSubPop.log
import random
pop = population(10, loci=[2, 3], infoFields=['Sex'])
InitSex(pop)
pop.setVirtualSplitter(sexSplitter())
# initialize male and females with different genotypes. Set initSex
# to False because this operator will by default also initialize sex.
InitByValue(pop, [[0]*5, [1]*5], subPops=([0, 0], [0, 1]), initSex=False)
# set Sex information field to 0 for all males, and 1 for all females
pop.setIndInfo([1], 'Sex', [0, 0])
pop.setIndInfo([2], 'Sex', [0, 1])
# Print individual genotypes, followed by values at information field Sex
Dump(pop, structure=False)
#end


#file log/accessIndividual.log
# create a population with two generations. The current generation has values
# 0-9 at information field x, the parental generation has values 10-19.
pop = population(size=[5, 5], loci=[2, 3], infoFields=['x'], ancGen=1)
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
#end


#file log/batchAccess.log
import random
pop = population(size=[4, 6], loci=[2], infoFields=['x'])
pop.setIndInfo([random.randint(0, 10) for x in range(10)], 'x')
pop.indInfo('x')
pop.setGenotype([0, 1, 2, 3], 0)
pop.genotype(0)
pop.setVirtualSplitter(infoSplitter(cutoff=[3], field='x'))
pop.setGenotype([0])    # clear all values
pop.setGenotype([5, 6, 7], [1, 1])
pop.indInfo('x', 1)
pop.genotype(1)
#end


#file log/popInfo.log
pop = population(10)
pop.setInfoFields(['a', 'b'])
pop.addInfoField('c')
pop.addInfoFields(['d', 'e'])
pop.infoFields()
#
cIdx = pop.infoIdx('c')
eIdx = pop.infoIdx('e')
pop.setIndInfo([1], cIdx)
for ind in pop.individuals():
    ind.setInfo(ind.info(cIdx) + 1, eIdx)

print pop.indInfo(eIdx)
#end          


#file log/ancestralPop.log
simu = simulator(population(500, loci=[1]), randomMating())
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        # start recording ancestral generations at generation 18
        setAncestralDepth(2, at=[-2]),
        stat(alleleFreq=[0], begin=-3),
        pyEval(r"'%.3f\n' % alleleFreq[0][0]", begin=-3)
    ],
    gen = 20
)
pop = simu.population(0)
# start from current generation
for i in range(pop.ancestralGens(), -1, -1):
  pop.useAncestralGen(i)
  Stat(pop, alleleFreq=[0])
  print '%d   %.3f' % (i, pop.dvars().alleleFreq[0][0])

# restore to the current generation  
pop.useAncestralGen(0)  
#end

#file log/addRemoveLoci.log
pop = population(10, loci=[3], chromNames=['chr1'])
# 1 1 1, 
pop.setGenotype([1])
# 1 1 1, 0 0 0
pop.addChrom(lociPos=[0.5, 1, 2], lociNames=['rs1', 'rs2', 'rs3'],
    chromName='chr2')
pop1 = population(10, loci=[3], chromNames=['chr3'],
    lociNames=['rs4', 'rs5', 'rs6'])
# 2 2 2,
pop1.setGenotype([2])
# 1 1 1, 0 0 0, 2 2 2
pop.addChromFrom(pop1)
# 1 1 1, 0 0 0, 2 0 2 2 0
pop.addLoci(chrom=[2, 2], pos=[1.5, 3.5], names=['rs7', 'rs8'])
# 1 1 1, 0 0 0, 2 0 2 0
pop.removeLoci([8])
Dump(pop)
#end

#file log/extract.log
import random
pop = population(size=[10, 10], loci=[5, 5],
    infoFields=['x', 'y'])
InitByValue(pop, range(10))
pop.setIndInfo([-1]*4 + [0]*3 + [-1]*3 + [2]*4 + [-1]*3 + [1]*4, 'x')
pop1 = pop.extract(field='x', loci=[1, 2, 3, 6, 7], infoFields=['x'])
Dump(pop1, structure=False)
#end

#file log/popVars.log
from pprint import pprint
pop = population(100, loci=[2])
InitByFreq(pop, [0.3, 0.7])
print pop.vars()    # No variable now
pop.dvars().myVar = 21
print pop.vars()
Stat(pop, popSize=1, alleleFreq=[0])
# pprint prints in a less messy format
pprint(pop.vars())
# print number of allele 1 at locus 0
print pop.vars()['alleleNum'][0][1]
# use the dvars() function to access dictionary keys as attributes
print pop.dvars().alleleNum[0][1]
print pop.dvars().alleleFreq[0]
#end


#file log/expression.log
simu = simulator(population(100, loci=[1]),
    randomMating(), 5)
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        stat(alleleFreq=[0]),
        terminateIf('alleleFreq[0][0] == 0. or alleleFreq[0][0] == 1.')
    ]
)
#end


#file log/savePop.log
pop = population(100, loci=[5], chromNames=['chrom1'])
pop.dvars().name = 'my population'
pop.save('sample.pop')
pop1 = LoadPopulation('sample.pop')
pop1.chromName(0)
pop1.dvars().name
#end

os.remove('sample.pop')

#file log/stageAndGen.log
simu = simulator(population(100, loci=[20]), randomMating())
simu.evolve(
    preOps = [initByFreq([0.2, 0.8])],
    ops = [
        stat(alleleFreq=[0], begin=80, step=10),
        pyEval(r"'After gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        pyEval(r"'Around gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1], stage=PrePostMating)
    ],
    postOps = [savePopulation(output='sample.pop')],
    gen=100
)
#end

os.remove('sample.pop')

#file log/dryrun.log
simu = simulator(population(100, loci=[20]), randomMating())
simu.evolve(
    preOps = [initByFreq([0.2, 0.8])],
    ops = [
        stat(alleleFreq=[0], begin=80, step=10),
        pyEval(r"'After gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        pyEval(r"'Around gen %d: alleleFreq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1], stage=PrePostMating)
    ],
    postOps = [savePopulation(output='sample.pop')],
    gen=100,
    dryrun = True
)
#end


#file log/replicate.log
simu = simulator(population(100, loci=[20]), randomMating(), 5)
simu.evolve(
    preOps = [initByFreq([0.2, 0.8])],
    ops = [
        stat(alleleFreq=[0], step=10),
        pyEval('gen', step=10, rep=0),
        pyEval(r"'\t%.2f' % alleleFreq[0][0]", step=10, rep=(0, 2, -1)),
        pyOutput('\n', step=10, rep=-1)
    ],
    gen=30,
)
#end


#file log/output.log
from simuPOP import *
simu = simulator(population(size=1000, loci=[2]), randomMating(), rep=3)
simu.evolve(
    preOps = [initByValue([1, 2, 2, 1])],  
    ops = [
        recombinator(rate=0.01),
        stat(LD=[0, 1]),
        pyEval(r"'%.2f\t' % LD[0][1]", step=20, output='>>LD.txt'),
        pyOutput('\n', rep=-1, step=20, output='>>LD.txt'),
        pyEval(r"'%.2f\t' % R2[0][1]", output='R2.txt'),
        pyEval(r"'%.2f\t' % LD[0][1]", step=20, output="!'>>LD_%d.txt' % rep"),
    ],
    gen=100
)
print open('LD.txt').read()
print open('R2.txt').read()    # Only the last write operation succeed.
print open('LD_2.txt').read()  # Each replicate writes to a different file.
#end


for file in ['LD.txt', 'LD_0.txt', 'LD_1.txt', 'LD_2.txt', 'R2.txt', 'LD_2.txt']:
    os.remove(file)

#file log/hybrid.log
def myPenetrance(geno):
    'A three-locus heterogeneity penetrance model'
    if sum(geno) < 2:
        return 0
    else:
        return sum(geno)*0.1

simu = simulator(population(1000, loci=[20]*3), randomMating())
simu.evolve(
    preOps = [initByFreq([0.8, 0.2])],
    ops = [
        pyPenetrance(func=myPenetrance, loci=[10, 30, 50]),
        stat(numOfAffected=True),
        pyEval(r"'%d: %d\n' % (gen, numOfAffected)")
    ],
    gen = 5
)
#end


#file log/pyOperator.log
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
            KamMutate(pop, maxAllele=1, rate=mu1, loci=[i])
        else:
            KamMutate(pop, maxAllele=1, rate=mu2, loci=[i])
    return True

simu = simulator(population(size=10000, loci=[2, 3]),
    randomMating())
simu.evolve(
    preOps = [ 
        initByFreq([.99, .01], loci=[0, 2, 4]),
        initByFreq([.8, .2], loci=[1, 3])],
    ops = [ 
        pyOperator(func=dynaMutator, param=(.2, 1e-2, 1e-5), stage=PreMating),
        stat(alleleFreq=range(5), step=10),
        pyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)                
#end

#file log/pyDuringMatingOperator.log
def rejectInd(off):
    'reject an individual if it off.allele(0) == 1'
    return off.allele(0) == 0

simu = simulator(population(size=100, loci=[1]),
    randomMating())
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [ 
        pyOperator(func=rejectInd, stage=DuringMating, offspringOnly=True),
    ],
    gen = 1
)
# You should see no individual with allele 1 at locus 0, ploidy 0.
simu.population(0).genotype()[:20]
#end


#file log/newOperator.log
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
                KamMutate(pop, maxAllele=1, rate=self.mu1, loci=[i])
            else:
                KamMutate(pop, maxAllele=1, rate=self.mu2, loci=[i])
        return True

simu = simulator(population(size=10000, loci=[2, 3]),
    randomMating())
simu.evolve(
    preOps = [ 
        initByFreq([.99, .01], loci=[0, 2, 4]),
        initByFreq([.8, .2], loci=[1, 3])],
    ops = [ 
        dynaMutator(cutoff=.2, mu1=1e-2, mu2=1e-5, stage=PreMating),
        stat(alleleFreq=range(5), step=10),
        pyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)          
#end


#file log/InitByFreq.log
def InitByFreq(pop, *args, **kwargs):
    initByFreq(*args, **kwargs).apply(pop)

InitByFreq(pop, [.2, .3, .5])
#end

#file log/migrSize.log
simu = simulator(
    population(size=[500, 1000], infoFields=['migrate_to']),
    randomMating())
simu.evolve(
    preOps = [initSex()],
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end

#file log/migrFixedSize.log
simu = simulator(
    population(size=[500, 1000], infoFields=['migrate_to']),
    randomMating(subPopSize=[500, 1000]))
simu.evolve(
    preOps = [initSex()],
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True, stage=PrePostMating),
        pyEval(r'"%s\n" % subPopSize', stage=PrePostMating)
    ],
    gen = 3
)
#end

#file log/demoFunc.log
def demo(gen, oldSize=[]):
    return [500 + gen*10, 1000 + gen*10]

simu = simulator(
    population(size=[500, 1000], infoFields=['migrate_to']),
    randomMating(subPopSize=demo))
simu.evolve(
    preOps = [initSex()],
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end

#file log/numOff.log
def checkNumOffspring(ms):
    '''Check the number of offspring for each family using
       information field father_idx
    '''
    simu = simulator(
        population(size=[30], infoFields=['father_idx', 'mother_idx']),
        matingScheme=ms)
    simu.evolve(
        preOps = [initSex()],
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
#end


#file log/sexMode.log
def checkSexMode(ms):
    '''Check the assignment of sex to offspring'''
    simu = simulator(
        population(size=[40]),
        matingScheme=ms)
    simu.evolve(preOps = [initSex()], ops=[], gen=1)
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
# Case 5: NumOfFamel (Specify number of female in each family)
checkSexMode(randomMating(
    numOffspring=(UniformDistribution, 4, 6),
    sexMode=(NumOfFemale, 2))
)
#end

#file log/monogamous.log
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
#end

#file log/polygamous.log
simu = simulator(population(100, infoFields=['father_idx', 'mother_idx']),
    polygamousMating(polySex=Male, polyNum=2))
simu.evolve(
    preOps = [initSex()],
    ops = [parentsTagger()],
    gen = 5
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals()][:20]
[ind.intInfo('mother_idx') for ind in pop.individuals()][:20]
#end

#file log/randomSelection.log
simu = simulator(population(100, ploidy=1, loci=[5, 5], ancGen=1,
    infoFields=['parent_idx']),
    randomSelection())
simu.evolve(
    preOps = [initByFreq([0.3, 0.7])],
    ops = [parentTagger()],
    gen = 5
)
pop = simu.extract(0)
ind = pop.individual(0)
par = pop.ancestor(ind.intInfo('parent_idx'), 1)
print ind.sex(), ind.genotype()
print par.sex(), par.genotype()
#end

#file log/alphaMating.log
simu = simulator(population(1000, loci=[5], 
    infoFields=['father_idx', 'mother_idx', 'fitness']),
    alphaMating(alphaSex=Male, alphaNum=2))
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [parentsTagger(),
        maSelector(loci=[0], fitness=[0.8, 0.8, 1]),
        stat(alleleFreq=[0]),
        pyEval(r'"%.2f\n" % alleleFreq[0][1]', step=5)
    ],
    gen = 20,
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals()][:10]
[ind.intInfo('mother_idx') for ind in pop.individuals()][:10]
#end

#file log/haplodiploidMating.log
pop = population(10, ploidy=Haplodiploid, loci=[5, 5],
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(sexSplitter())
simu = simulator(pop, haplodiploidMating())
simu.evolve(
    preOps = [initSex(),
        initByValue([0]*10, subPops=[(0, 0)], initSex=False),
        initByValue([1]*10+[2]*10, subPops=[(0, 1)], initSex=False)],
    ops = [parentsTagger(),
        dumper(structure=False, stage=PrePostMating)],
    gen = 1
)
#end

#file log/selfMating.log
pop = population(20, loci=[8])
# every chromosomes are different. :-)
for idx, ind in enumerate(pop.individuals()):
    ind.setGenotype([idx*2], 0)
    ind.setGenotype([idx*2+1], 1)

simu = simulator(pop, selfMating())
simu.evolve(
    ops = [recombinator(rate=0.1)],
    gen = 1
)
Dump(simu.population(0), width=3, structure=False, max=10)
#end

#file log/heteroMatingSP.log
pop = population(size=[1000, 1000], loci=[2],
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, heteroMating(
    [randomMating(numOffspring=2, subPop=0),
     randomMating(numOffspring=4, subPop=1)
    ])
)
simu.evolve(
    preOps = [initSex()],
    ops= [parentsTagger()],
    gen=10
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals(0)][:10]
[ind.intInfo('father_idx') for ind in pop.individuals(1)][:10]
#end

#file log/heteroMatingVSP.log
pop = population(size=[1000], loci=[2],
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(proportionSplitter([0.2, 0.8]))
simu = simulator(pop, heteroMating(
    matingSchemes = [
        selfMating(subPop=(0, 0)),
        randomMating(subPop=(0, 1))
    ])
)
simu.evolve(
    preOps = [initSex()],
    ops= [parentsTagger()],
    gen = 10
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals(0)][:15]
[ind.intInfo('mother_idx') for ind in pop.individuals(0)][:15]
#end

#file log/heteroMatingWeight.log
pop = population(size=[1000], loci=[2],
    infoFields=['mark'])
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
        randomMating(subPop=0, weight=-0.5, ops=[markOff(0)]),
        randomMating(subPop=(0, 0), weight=2, ops=[markOff(1)]),
        randomMating(subPop=(0, 1), weight=3, ops=[markOff(2)])
    ])
)
simu.evolve(
    preOps = [initSex()],
    ops= [],
    gen = 10
)
marks = list(simu.extract(0).indInfo('mark'))
marks.count(0.)
marks.count(1.)
marks.count(2.)
#end


#file log/randomMating.log
def mendelianOffspringGenerator(ops=[], *args, **kwargs):
    'An offspring generator that uses mendelianGenoTransmitter()'
    return  offspringGenerator([mendelianGenoTransmitter()] + ops, *args, **kwargs)

def randomMating(numOffspring = 1., sexMode = RandomSex, ops = [], subPopSize = [],
        subPop = (), weight = 0, selectionField = 'fitness'):
    'A basic diploid sexual random mating scheme.'
    return homoMating(
        chooser = randomParentsChooser(True, selectionField),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)
#end


#file log/sequentialSelfing.log
simu = simulator(population(100, loci=[5]*3, infoFields=['parent_idx']),
    homoMating(sequentialParentChooser(), selfingOffspringGenerator()))
simu.evolve(
    preOps = [initByFreq([0.2]*5)],
    ops = [
        parentTagger(),
        dumper(structure=False, stage=PrePostMating, max=5)],
    gen = 1
)
#end

#file log/controlledOffGenerator.log
def traj(gen):
    return [0.5 + gen * 0.01]

simu = simulator(population(1000, loci=[10]*2),
    homoMating(randomParentChooser(),
        controlledOffspringGenerator(loci=[5],
            alleles=[0], freqFunc=traj,
            ops = [selfingGenoTransmitter()]))
)

# evolve the population while keeping allele frequency 0.5
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [stat(alleleFreq=[5, 15]),
        pyEval(r'"%.2f\t%.2f\n" % (alleleFreq[5][0], alleleFreq[15][0])')],
    gen = 5
)

#end

#file log/mitochondrial.log
pop = population(10, loci=[5]*5,
    # one autosome, two sex chromosomes, and two mitochondrial chromosomes
    chromTypes=[Autosome, ChromosomeX, ChromosomeY] + [Customized]*2,
    infoFields=['father_idx', 'mother_idx'])

simu = simulator(pop, randomMating(ops=[mitochondrialGenoTransmitter()]))

simu.evolve(
    preOps=[initByFreq([0.4] + [0.2]*3)],
    ops=[
        recombinator(rate=0.1),
        parentsTagger(),
        dumper(structure=False),
    ],
    gen = 2
)
#end


#file log/sexSpecificRec.log
class sexSpecificRecombinator(pyOperator):
    def __init__(self, intensity=0, rate=0, loci=[], convMode=NoConversion,
            maleIntensity=0, maleRate=0, maleLoci=[], maleConvMode=NoConversion,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.recombinator = recombinator(intensity, rate, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = recombinator(maleIntensity, maleRate,
            maleLoci, maleConvMode)
        #
        self.initialized = False
        # Note the use of parameter formOffGenotype.
        pyOperator.__init__(self, func=self.transmitGenotype,
            stage=DuringMating, formOffGenotype=True, *args, **kwargs)
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
simu = simulator(pop, randomMating())
simu.evolve(
    preOps=[initByFreq([0.4] + [0.2]*3)],
    ops=[
        parentsTagger(),
        sexSpecificRecombinator(rate=0.1, maleRate=0),
        dumper(structure=False),
    ],
    gen = 2
)
#end


#file log/infoChooser.log
pop = population(100, loci=[10],
    infoFields=['father_idx', 'mother_idx', 'sibling'])

pop.setVirtualSplitter(sexSplitter())

def locate_sibling(pop):
    '''The population is arranged as MFMFMFMF... where MF are siblings, so the
    sibling of males are 1, 3, 5, .. and the slibling of females are 0, 2, 4, ...
    '''
    pop.setIndInfo([2*x+1 for x in range(pop.popSize()/2)], 'sibling', (0, 0))
    pop.setIndInfo([2*x for x in range(pop.popSize()/2)], 'sibling', (0, 1))

simu = simulator(pop, consanguineousMating(func=locate_sibling, infoFields=['sibling'],
    numOffspring=2, sexMode=(NumOfMale, 1)))
simu.evolve(
    preOps = [initByFreq([0.2, 0.8], sex=[Male, Female])],
    ops = [
        parentsTagger(),
        dumper(structure=False, max=6, at=[-1])
    ],
    gen = 2
)
#end


#file log/generator.log
def func():
    i = 1
    all = 0
    while i <= 5:
        all += 1./i
        i += 1
        yield all 

for i in func():
    print '%.3f' % i,

#end

#file log/pyParentsChooser.log
from random import randint
def randomChooser(pop, sp):
    males = []
    females = []
    # identify males and females in each social rank
    for rank in range(3):
        males.append([x for x in range(pop.subPopSize(sp)) \
            if pop.individual(x, sp).sex() == Male and \
                pop.individual(x, sp).info('rank') == rank])
        females.append([x for x in range(pop.subPopSize(sp)) \
            if pop.individual(x, sp).sex() == Female and \
                pop.individual(x, sp).info('rank') == rank])
    while True:
        # choose a parent randomly
        idx = randint(0, pop.subPopSize(sp) - 1)
        par = pop.individual(idx, sp)
        # then choose a spouse in the same rank randomly
        if par.sex() == Male:
            rank = par.intInfo('rank')
            yield idx, females[rank][randint(0, len(females[rank]) - 1)]
        else:
            rank = par.intInfo('rank')
            yield males[rank][randint(0, len(males[rank]) - 1)], idx

def setRank(pop, dad, mom, off):
    'The rank of offspring can increase or drop to zero randomly'
    off.setInfo((dad.info('rank') + randint(-1, 1)) % 3, 'rank')

pop = population(size=[1000, 2000], loci=[1], infoFields=['rank'])
pop.setIndInfo([randint(0, 2) for x in range(pop.popSize())], 'rank')

simu = simulator(pop, homoMating(
    pyParentsChooser(randomChooser),
    mendelianOffspringGenerator()))
simu.evolve(
    preOps = [initSex()],
    ops = [],
    gen = 5
)    

#end

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

sys.path.append('log')

try:
    import myParentsChooser
except:
    os.chdir('log')
    os.system('python setup.py build_ext --swig-opts="-O -templatereduce -shadow -c++ -keyword -nodefaultctor" install --install-purelib="." --install-platlib="."')
    os.chdir('..')

#file log/cppParentChooser.log
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

pop = population(100, loci=[1])
simu = simulator(pop,
    homoMating(pyParentsChooser(parentsChooser), mendelianOffspringGenerator())
)
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [],
    gen = 100
)

#end


#file log/simuGen.log
simu = simulator(population(50, loci=[10], ploidy=1),
    randomSelection(), rep=3)
simu.evolve(ops = [], gen = 5)
simu.gen()
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        stat(alleleFreq=[5]),
        ifElse('alleleNum[5][0] == 0',
            pyEval(r"'Allele 0 is lost in rep %d at gen %d\n' % (rep, gen)")),
        ifElse('alleleNum[5][0] == 50',
            pyEval(r"'Allele 0 is fixed in rep %d at gen %d\n' % (rep, gen)")),
        terminateIf('alleleNum[5][0] == 0 or alleleNum[5][0] == 50'),
    ],
)
simu.gen()
#end


#file log/twoStage.log
# First stage: use the standard random mating scheme, do not use any
# information field for efficiency considerations.
simu = simulator(population(500, loci=[10]), randomMating())
simu.evolve(preOps = [initByFreq([0.5, 0.5])],
    ops = [], gen = 50)
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
#end

#file log/changeStru.log
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
    preOps = [initSex()],
    ops = [pyOperator(func=mutator, param=(10000, 2e-6))],
    gen = 200
)
for pop in simu.populations():
    print pop.totNumLoci(), pop.lociPos()

#end


#file log/simuFunc.log
simu = simulator(population(100, loci=[5, 10], infoFields=['x']),
    randomMating(), rep=5)
simu.evolve(preOps=[initByFreq([0.4, 0.6])],
    ops=[], gen=10)
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
#end

os.remove('sample.sim')


#file log/initSex.log
pop = population(size=[1000, 1000])
InitSex(pop, maleFreq=0.3, subPops=0)
InitSex(pop, sex=[Male, Female, Female], subPops=1)
Stat(pop, numOfMale=True)
print pop.dvars(0).numOfMale
print pop.dvars(1).numOfMale
#end

#file log/initByFreq.log
pop = population(size=[2, 3], loci=[5, 7])
InitByFreq(pop, alleleFreq=[[.2, .8], [.8, .2]])
Dump(pop, structure=False)
#end

#file log/initByFreqIdenticalInds.log
pop = population(size=[2, 3], loci=[5, 7])
InitByFreq(pop, alleleFreq=[.2, .8], identicalInds=True)
Dump(pop, structure=False)
#end

#file log/initByValue.log
pop = population(size=[2, 3], loci=[5, 7])
InitByValue(pop, [1]*5 + [2]*7 + [3]*5 +[4]*7)
Dump(pop, structure=False)
#end

#file log/initByValueProp.log
pop = population(size=[6, 8], loci=[5, 7])
pop.setVirtualSplitter(sexSplitter())
# initialize sex and the first two loci
InitByValue(pop, loci=range(5), value=range(10))
# initialize all males
InitByValue(pop, loci=range(5, 12), value=[2]*7,
    subPops=[(0, 0), (1, 0)], initSex=False)
# initialize females by proportion
InitByValue(pop, loci=range(5, 12), ploidy=1, value=[[3]*7, [4]*7],
    initSex=False, subPops=[(0, 1), (1, 1)], proportions=[0.4, 0.6])
Dump(pop, structure=False)
#end

#file log/dumper.log
pop = population(size=[10, 10], loci=[20, 30], infoFields=['gen'],
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

#end

#file log/savePopulation.log
simu = simulator(population(100, loci=[2]),
    randomMating(), rep=5)
simu.evolve(
    preOps = [initByFreq([0.2, 0.8])],
    ops = [
        savePopulation(output="!'snapshot_%d_%d.pop' % (rep, gen)",
            step = 10),
        ],
    gen = 50
)
#end

for rep in range(5):
    for gen in range(0, 50, 10):
        os.remove('snapshot_%d_%d.pop' % (rep, gen))

#file log/setAncDepth.log
simu = simulator(population(100, infoFields=['father_idx', 'mother_idx']),
    randomMating(), rep=5)
simu.evolve(
    preOps = [initByFreq([0.3, 0.7])],
    ops = [
        setAncestralDepth(2, at=-2),
        parentsTagger(begin=-2)
    ],
    gen = 100
)
pop = simu.population(3)
print pop.ancestralGens()
print pop.ancestor(10, 1).info('father_idx')
#end


#file log/ifElse.log
simu = simulator(
    population(size=1000, loci=[1]),
    randomMating(), rep=4)
simu.evolve(
    preOps = [
        initByFreq([0.5, 0.5]),
        pyExec('below40, above60 = 0, 0')
    ],
    ops = [
        stat(alleleFreq=[0]),
        ifElse('alleleFreq[0][1] < 0.4',
            pyExec('below40 += 1')),
        ifElse('alleleFreq[0][1] > 0.6',
            pyExec('above60 += 1')),
        ifElse('alleleFreq[0][1] == 0 or alleleFreq[0][1] == 1',
            pyExec('stoppedAt = gen')),
        terminateIf('alleleFreq[0][1] == 0 or alleleFreq[0][1] == 1')
    ]
)
for pop in simu.populations():
    print 'Overall: %4d, below 40%%: %4d, above 60%%: %4d' % \
        (pop.dvars().stoppedAt, pop.dvars().below40, pop.dvars().above60)

#end

#file log/debug.log
simu = simulator(population(100, loci=[1]), randomMating(), rep=5)
simu.evolve(
    preOps = [initByFreq([0.1, 0.9])],
    ops = [
        stat(alleleFreq=[0]),
        ifElse('alleleNum[0][0] == 0',
            turnOnDebug(DBG_MUTATOR),
            turnOffDebug(DBG_MUTATOR)),
        ifElse('alleleNum[0][0] == 0',
            pointMutator(loci=[0], toAllele=0, inds=[0])),
    ],
    gen = 100
)
#end

#file log/pause.log
simu = simulator(population(100), randomMating(), rep=10)
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [pause(stopOnKeyStroke=str(x), rep=x) for x in range(10)],
    gen = 100
)
#end

#file log/ticToc.log
simu = simulator(population(10000, loci=[100]*5), randomMating(), rep=2)
simu.evolve(
    preOps = [initByFreq([0.1, 0.9])],
    ops = [
        stat(alleleFreq=[0]),
        ticToc(step=50, rep=-1),
    ],
    gen = 101
)
#end

#file log/pyExec.log
simu = simulator(population(100, loci=[1]),
    randomMating(), rep=2)
simu.evolve(
    preOps = [
        initByFreq([0.2, 0.8]),
        pyExec('traj=[]')
    ],
    ops = [
        stat(alleleFreq=[0]),
        pyExec('traj.append(alleleFreq[0][1])'),
    ],
    gen=5
)
# print trajectory
print ', '.join(['%.3f' % x for x in simu.dvars(0).traj])
#end

#file log/pyEval.log
simu = simulator(population(1000, loci=[1],
    infoFields=['mother_idx', 'father_idx']),
    randomMating())
simu.evolve(
    preOps = [initSex()],
    ops = [
        stat(alleleFreq=[0]),
        parentsTagger(),
        pyEval(r'"gen %d, #father %d, #mother %d\n"' \
            ' % (gen, numFather, numMother)',
            stmts="numFather = len(set(pop.indInfo('father_idx')))\n"
                "numMother = len(set(pop.indInfo('mother_idx')))",
            exposePop='pop')
    ],
    gen=3
)
#end

#file log/infoEval.log
import random
pop = population(20, loci=[1], infoFields=['a'])
pop.setVirtualSplitter(infoSplitter('a', cutoff=[3]))
InitByFreq(pop, [0.2, 0.8])
pop.setIndInfo([random.randint(2, 5) for x in range(20)], 'a')
InfoEval(pop, 'a', subPops=[(0, 0)]);print
InfoEval(pop, 'ind.allele(0, 0)', exposeInd='ind');print
# use population variables
pop.dvars().b = 5
InfoEval(pop, '"%d " % (a+b)', usePopVars=True);print
#end

#file log/infoExec.log
pop = population(100, loci=[1], infoFields=['a', 'b', 'c'])
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
#end

#file log/migrateByProb.log
simu = simulator(
    population(size=[1000]*3, infoFields=['migrate_to']),
    randomMating())
simu.evolve(
    preOps = [initSex()],
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
#end

#file log/migrateByPropAndCount.log
simu = simulator(
    population(size=[1000]*3, infoFields=['migrate_to']),
    randomMating())
simu.evolve(
    preOps = [initSex()],
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
#end

#file log/migrateVSP.log
pop = population(size=[1000]*2, infoFields=['migrate_to'])
pop.setVirtualSplitter(sexSplitter())
simu = simulator(pop, randomMating())
simu.evolve(
    # 500 males and 500 females
    preOps = [initSex(sex=[Male, Female])],
    ops = [
        migrator(rate=[
            [0, 0.10],
            [0, 0.05],
            ],
            mode = ByProportion,
            subPops=[(0, 0), (0, 1)]),
        stat(popSize=True, numOfMale=True, stage=PrePostMating),
        pyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMale'], subPopSize[0], "
            "subPop[1]['numOfMale'], subPopSize[1])", stage=PrePostMating),
    ],
    gen = 2
)   
#end

#file log/manualMigration.log
pop = population([10]*2, infoFields=['migrate_to'])
pop.setIndInfo([0, 1, 2, 3]*5, 'migrate_to')
Migrate(pop, mode=ByIndInfo)
pop.subPopSizes()
#end

#file log/getParam.log
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

pars = simuOpt.simuOpt(options, 'A demo simulation')
print pars.usage()
# You can manually feed parameters...
pars.processArgs(['--rep=10'])
pars.rep
#beginignore
oldArg = sys.argv[-1]
sys.argv.pop()
if not os.path.isfile('getParam.png'):
    print 'Run a GUI if getParam has not been runned'
else:
    sys.argv.extend(['--rate=[0.25]', '--rep=5', '--pop="CEU"'])
    simuOpt.setOptions(gui=False)
    pars.processArgs(sys.argv)

#endignore
# but simuOpt.getParam is the easiest to used
if not pars.getParam():
    sys.exit(1)

#beginignore
sys.argv[1] = oldArg
#endignore
# save prameters to a configuration file
pars.saveConfig('sample.cfg')
# post-process parameters
pars.rate
pars.rep
pars.rate = pars.rate * pars.rep
# extract parameters as a dictionary or a list
pars.asDict()
pars.asList()
#end


################################################

#file log/splitMerge.log
pop = population(1000, loci=[1], infoFields=['migrate_to'])
simu = simulator(pop, randomSelection())
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.2, 0.8], at = 3),
        splitSubPop(1, proportions=[0.4, 0.6], at = 5),
        mergeSubPops([0, 2], at = 7),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize'),
    ],
    gen = 10
)
#end

#file log/splitMigration.log
from simuPOP import *
pop = population(1000, loci=[1], infoFields=['migrate_to'])
simu = simulator(pop, randomSelection())
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.2, 0.3, 0.5], at = 3),
        migrator(rate = [0.2], subPops=0, toSubPops=1, 
            begin = 3, end = 4),
        migrator(rate = [
            [0, 0.2, 0.4],
            [0, 0,   0.1],
            [0.1, 0.1, 0]],
            begin = 4),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize'),
    ],
    gen = 10
)
#end

#file log/splitMigration2.log
from simuPOP import *
pop = population(1000, loci=[1], infoFields=['migrate_to'])
def popSize(gen, oldSize=[]):
    if gen < 3:
        return [1000]
    elif gen < 5:
        return [400, 500]
    else:
        return [300, 400, 600]

simu = simulator(pop, randomSelection(subPopSize=popSize))
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.3, 0.7], at = 3),
        migrator(rate = [0.2], subPops=0, toSubPops=1, 
            begin = 3, end = 4),
        splitSubPop(0, proportions=[0.3, 0.7], at = 5),
        migrator(rate = [
            [0, 0.2, 0.4],
            [0, 0,   0.1],
            [0.1, 0.1, 0]],
            begin = 5),
        stat(popSize=True, stage=PreMating),
        pyEval(r'"From %s\t" % subPopSize', stage=PreMating),
        stat(popSize=True),
        pyEval(r'"to: %s\n" % subPopSize'),
    ],
    gen = 10
)
#end

#file log/splitMigration3.log
from simuPOP import *
pop = population(1000, loci=[1], infoFields=['migrate_to'])
def popSize(gen, oldSize=[]):
    return [x*2 for x in oldSize]

simu = simulator(pop, randomSelection(subPopSize=popSize))
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.3, 0.7], at = 3),
        migrator(rate = [0.2], subPops=0, toSubPops=1, 
            begin = 3, end = 4),
        splitSubPop(0, proportions=[0.3, 0.7], at = 5),
        migrator(rate = [
            [0, 0.2, 0.4],
            [0, 0,   0.1],
            [0.1, 0.1, 0]],
            begin = 5),
        stat(popSize=True, stage=PrePostMating),
        pyEval(r'"From %s\t" % subPopSize', stage=PreMating),
        pyEval(r'"to: %s\n" % subPopSize'),
    ],
    gen = 10
)
#end

#file log/migrator.log
from simuUtil import *

#number of cities
nc = 6

# how to change subpop size?
def changeSPSize(gen, oldSize):
    size = [0]*len(oldSize)
    for i in range(0, len(size)):
        size[i] = oldSize[i]*1.2
    if size[i] > 1000:
        size[i] /= 2
    return size

# migration between subpopulaitons
rates = []
for i in range(nc):
    rates.append([0.]*nc)

for i in range(1, nc-1):
    rates[i][i+1]=0.05
    rates[i][i-1]=0.05

#
rates[0][1] = 0.1
rates[nc-1][nc-2] = 0.1

# print rates
print rates
migr = migrator(rate=rates, mode=ByProbability)

# initially, we need to set everyone to middle subpop
initMigr = migrator(rate=[[1]], mode=ByProportion,
    subPops=[0], toSubPops=[nc/2])

pop = population(size=[500]*nc, infoFields=['migrate_to'])

# the new popsize relies on a variable newSPSize
# which is calculated from subPopSize bu newSize operator
simu = simulator(pop,
    randomMating(subPopSize=changeSPSize) )

# evolve!
simu.evolve(
    preOps = [initSex(), initMigr ],
    ops = [
        migr,
        stat(popSize=True),
        pyEval('list(subPopSize)'),
        pyOutput('\n', rep=-1)
    ],
    gen=10
)

#end


#file log/kamMutator.log
pop = population(size=[200, 300], ploidy=2, loci=[5, 10],
    lociPos=[range(0, 5), range(0, 20, 2)],
    alleleNames=['A', 'C', 'T', 'G'])
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
    preOps = [initByFreq([.8, .2])],
    ops = [
        stat(alleleFreq=[0, 1], Fst=[1], step=10),
        kamMutator(rate=0.001, rep=1),
        kamMutator(rate=0.0001, rep=2)
    ],
    gen=10
)
#end
 
#file log/smmMutator.log
simu = simulator(population(size=300, loci=[3, 5]), randomMating())
simu.evolve(
    preOps = [initByFreq( [.2, .3, .5])],
    ops = [
        smmMutator(rate=1,  incProb=.8),
    ],
    gen=1
)
#end

#file log/gsmMutator.log
import random
simu = simulator(population(size=300, loci=[3, 5]), randomMating())
simu.evolve(
    preOps = [initByFreq( [.2, .3, .5])],
    ops = [
        gsmMutator(rate=1, p=.8, incProb=.8),
    ],
    gen=1
)

def rndInt():
  return random.randrange(3, 6)

simu.evolve(
    preOps = [initByFreq( [.2, .3, .5])],
    ops = [
        gsmMutator(rate=1, func=rndInt, incProb=.8),
    ],
    gen=1
)
#end


#file log/pyMutator.log
import random
simu = simulator(population(size=300, loci=[3, 5]), randomMating())

def mutateTo(allele):
  return allele + random.randrange(3, 6)

simu.evolve(
    preOps = [initByValue([50]*3 + [100]*5)],
    ops = [
        pyMutator(rate=0.001, func=mutateTo),
    ],
    gen=1
)
#end

#### #file log/recombinator.log
simu = simulator(population(4, loci=[4, 5, 6], 
    infoFields=['father_idx', 'mother_idx']),
    randomMating())
simu.evolve(
    preOps = [initByFreq([.2, .2, .4, .2]), dumper(structure=False) ],
    ops = [parentsTagger()],
    postOps = [ dumper(structure=False)],
    gen=1
)
simu.evolve(
    ops = [
        parentsTagger(),
        recombinator(rate=[1, 1, 1], loci=[2, 6, 10])
    ],
    postOps = [dumper(structure=False)],
    gen=1
)
#end




#file log/mapSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
    randomMating())
s1 = .1
s2 = .2
simu.evolve(
    preOps=[initByFreq(alleleFreq=[.2, .8])],
    ops = [
        stat(alleleFreq=[0], genoFreq=[0]),
        mapSelector(loci=0, fitness={'0-0':(1-s1), '0-1':1, '1-1':(1-s2)}),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    gen=300
)
#end

 
#file log/maSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
    randomMating())
s1 = .1
s2 = .2
simu.evolve(
    preOps=[initByFreq(alleleFreq=[.2, .8])],
    ops = [
        stat(alleleFreq=[0], genoFreq=[0]),
        maSelector(loci=0, fitness=[1-s1, 1, 1-s2]),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    gen = 300)
#end

 
#file log/mlSelector.log
simu = simulator(
    population(size=10, ploidy=2, loci=[2], 
    infoFields=['fitness', 'spare']),
    randomMating())
simu.evolve(
    [ mlSelector([
         mapSelector(loci=0, fitness={'0-0':1, '0-1':1, '1-1':.8}),
         mapSelector(loci=1, fitness={'0-0':1, '0-1':1, '1-1':.8}),
         ], mode = Additive),
    ],
    preOps = [
        initByFreq(alleleFreq=[.2, .8])
    ],
    gen = 2
)
#end



#file log/pySelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[3], infoFields=['fitness'] ),
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
    preOps=[initByFreq(alleleFreq=[.2, .8])],
    ops = [
        stat(alleleFreq=[0], genoFreq=[0]),
        pySelector(loci=[0, 1], func=sel),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=25)
    ],
    gen=100
)
#end


#file log/mapPenetrance.log
pop = population(size=[2,8], ploidy=2, loci=[2] )
InitByFreq(pop, [.2, .8])
MapPenetrance(pop, loci=0, 
    penetrance={'0-0':0, '0-1':1, '1-1':1})
Stat(pop, numOfAffected=1)
#end


#file log/mlPenetrance.log
pop = population(1000, loci=[3])
InitByFreq(pop, [0.3, 0.7])
pen = []
for loc in (0, 1, 2):
    pen.append(maPenetrance(loci=loc, wildtype=[1],
        penetrance=[0, 0.3, 0.6] ) )

# the multi-loci penetrance
MlPenetrance(pop, mode=Multiplicative, peneOps=pen)
Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#end

#file log/pyPenetrance.log
pop = population(1000, loci=[3])
InitByFreq(pop, [0.3, 0.7])
def peneFunc(geno):
    p = 1
    for l in range(len(geno)/2):
        p *= (geno[l*2]+geno[l*2+1])*0.3
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
    return func

# then, given a table, you can do
PyPenetrance(pop, loci=(0, 1, 2),
    func=peneFunc( ((0, 0.5), (0.3, 0.8)) ) )
#end





## 
## 
#file log/reichParam.py
initSize =  10000            # initial population size
finalSize = 1000000          # final population size
burnin = 500                 # evolve with constant population size
endGen = 1000                # last generation
mu = 3.2e-5                  # mutation rate
C_f0 = 0.2                   # initial allelic frequency of *c*ommon disease
R_f0 = 0.001                 # initial allelic frequency of *r*are disease
max_allele = 255             # allele range 1-255 (1 for wildtype)
C_s = 0.0001                 # selection on common disease
R_s = 0.9                    # selection on rare disease
psName = 'lin_exp'           # filename of saved figures 

# allele spectrum
C_f = [1-C_f0] + [x*C_f0 for x in [0.9, 0.02, 0.02, 0.02, 0.02, 0.02]]
R_f = [1-R_f0] + [x*R_f0 for x in [0.9, 0.02, 0.02, 0.02, 0.02, 0.02]]
#end

#file log/reichSimulator.py
from simuOpt import setOptions
setOptions(alleleType='long')
from simuPOP import *

# instantaneous population growth
def ins_exp(gen, oldSize=[]):
    if gen < burnin:
        return [initSize]
    else:
        return [finalSize]

def simulate(incScenario):
    simu = simulator(                                        # create a simulator
        population(size=incScenario(0), loci=[1,1],
            infoFields=['fitness']),                         # inital population
        randomMating(subPopSize=incScenario)           # random mating
    )

#simulate(ins_exp)
#end

#file log/reichMutSel.py
def simulate(incScenario):
    simu = simulator(                       # create a simulator
        population(size=incScenario(0), loci=[1,1],
            infoFields=['fitness']),        # inital population
        randomMating(subPopSize=incScenario)
    )
    simu.evolve(                            # start evolution
        preOps=[                            # operators that will be applied before evolution
            # initialize locus 0 (for common disease)
            initByFreq(loci=[0], alleleFreq=C_f),
            # initialize locus 1 (for rare disease)
            initByFreq(loci=[1], alleleFreq=R_f),
        ],
        ops=[                               # operators that will be applied at each gen
            # mutate: k-alleles mutation model
            kamMutator(rate=mu, maxAllele=max_allele),
            # selection on common and rare disease,
            mlSelector([                # multiple loci - multiplicative model
                maSelector(loci=0, fitness=[1,1,1-C_s], wildtype=[0]),
                maSelector(loci=1, fitness=[1,1,1-R_s], wildtype=[0])
            ], mode=SEL_Multiplicative),
        ],
        gen = endGen
    )

#simulate(ins_exp)
#end

#file log/reichStat.py
def ne(pop):
    ' calculate effective number of alleles '
    Stat(pop, alleleFreq=[0,1])
    f0 = [0, 0]
    ne = [0, 0]
    for i in range(2):
        freq = pop.dvars().alleleFreq[i][1:]
        f0[i] = 1 - pop.dvars().alleleFreq[i][0]
        if f0[i] == 0:
            ne[i] = 0
        else:
            ne[i] = 1. / sum([(x/f0[i])**2 for x in freq])
    print '%d\t%.3f\t%.3f\t%.3f\t%.3f' % (pop.gen(), f0[0], f0[1], ne[0], ne[1])
    return True

#end

#file log/reichOpt.py
options = [
    {'arg': 'h',
     'longarg': 'help',
     'default': False, 
     'description': 'Print this usage message.',
     'allowedTypes': [types.NoneType, type(True)],
     'jump': -1                    # if -h is specified, ignore any other parameters.
    },
    {'longarg': 'initSize=',
     'default': 10000,
     'label': 'Initial population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial population size. This size will be maintained
                till the end of burnin stage''',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'finalSize=',
     'default': 1000000,
     'label': 'Final population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Ending population size (after expansion.',
     'validate': simuOpt.valueGT(0)
    }, 
    {'longarg': 'burnin=',
     'default': 500,
     'label': 'Length of burn-in stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of generations of the burn in stage.',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'endGen=',
     'default': 1000,
     'label': 'Last generation',
     'allowedTypes': [types.IntType],
     'description': 'Ending generation, should be greater than burnin.',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'growth=',
     'default': 'instant',
     'label': 'Population growth model',
     'description': '''How population is grown from initSize to finalSize.
                Choose between instant, linear and exponential''',
     'chooseOneOf': ['linear', 'instant'],
    },
    {'longarg': 'name=',
     'default': 'cdcv',
     'allowedTypes': [types.StringType],
     'label': 'Name of the simulation',
     'description': 'Base name for configuration (.cfg) log file (.log) and figures (.eps)'
    },
]

def getOptions(details=__doc__):
    # get all parameters, __doc__ is used for help info
    allParam = simuOpt.getParam(options, 
        'This program simulates the evolution of a common and a rare direse\n' +
        'and observe the evolution of allelic spectra\n', details)
    #
    # when user click cancel ...
    if len(allParam) == 0:
        sys.exit(1)
    # -h or --help
    if allParam[0]:    
        print simuOpt.usage(options, __doc__)
        sys.exit(0)
    # automatically save configurations
    name = allParam[-1]
    if not os.path.isdir(name):
        os.makedirs(name)
    simuOpt.saveConfig(options, os.path.join(name, name+'.cfg'), allParam)
    # return the rest of the parameters
    return allParam[1:-1]

#
# IGNORED
# 

if __name__ == '__main__':
    # get parameters
    (initSize, finalSize, burnin, endGen, growth) = getOptions()
    # 
    from simuPOP import *
    #
    if initSize > finalSize:
        print 'Initial size should be greater than final size'
        sys.exit(1)
    if burnin > endGen:
        print 'Burnin gen should be less than ending gen'
        sys.exit(1)
    if growth == 'linear':
        simulate(lin_exp)
    else:
        simulate(ins_exp)

#end

#file log/reich.py
#!/usr/bin/env python

'''
simulation for Reich(2001):
     On the allelic spectrum of human disease

'''

import simuOpt
simuOpt.setOptions(alleleType='long')

import sys, types, os

options = [
    {'arg': 'h',
     'longarg': 'help',
     'default': False, 
     'description': 'Print this usage message.',
     'allowedTypes': [types.NoneType, type(True)],
     'jump': -1                    # if -h is specified, ignore any other parameters.
    },
    {'longarg': 'initSize=',
     'default': 10000,
     'label': 'Initial population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial population size. This size will be maintained
                till the end of burnin stage''',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'finalSize=',
     'default': 1000000,
     'label': 'Final population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Ending population size (after expansion.',
     'validate': simuOpt.valueGT(0)
    }, 
    {'longarg': 'burnin=',
     'default': 500,
     'label': 'Length of burn-in stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of generations of the burn in stage.',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'endGen=',
     'default': 1000,
     'label': 'Last generation',
     'allowedTypes': [types.IntType],
     'description': 'Ending generation, should be greater than burnin.',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'growth=',
     'default': 'instant',
     'label': 'Population growth model',
     'description': '''How population is grown from initSize to finalSize.
                Choose between instant, linear and exponential''',
     'chooseOneOf': ['linear', 'instant'],
    },
    {'longarg': 'name=',
     'default': 'cdcv',
     'allowedTypes': [types.StringType],
     'label': 'Name of the simulation',
     'description': 'Base name for configuration (.cfg) log file (.log) and figures (.eps)'
    },
]

def getOptions(details=__doc__):
    # get all parameters, __doc__ is used for help info
    allParam = simuOpt.getParam(options, 
      'This program simulates the evolution of a common and a rare direse\n' +
        'and observe the evolution of allelic spectra\n', details)
    #
    # when user click cancel ...
    if len(allParam) == 0:
        sys.exit(1)
    # -h or --help
    if allParam[0]:    
        print simuOpt.usage(options, __doc__)
        sys.exit(0)
    # automatically save configurations
    name = allParam[-1]
    if not os.path.isdir(name):
        os.makedirs(name)
    simuOpt.saveConfig(options, os.path.join(name, name+'.cfg'), allParam)
    # return the rest of the parameters
    return allParam[1:-1]


# these can be put as options as well.
mu = 3.2e-5                  # mutation rate
C_f0 = 0.2                   # initial allelic frequency of *c*ommon disease
R_f0 = 0.001                 # initial allelic frequency of *r*are disease
max_allele = 255             # allele range 1-255 (1 for wildtype)
C_s = 0.0001                 # selection on common disease
R_s = 0.9                    # selection on rare disease

C_f = [1-C_f0] + [x*C_f0 for x in [0.9, 0.02, 0.02, 0.02, 0.02, 0.02]]
R_f = [1-R_f0] + [x*R_f0 for x in [0.9, 0.02, 0.02, 0.02, 0.02, 0.02]]

# instantaneous population growth
def ins_exp(gen, oldSize=[]):
    if gen < burnin:
        return [initSize]
    else:
        return [finalSize]

# linear growth after burn-in
def lin_exp(gen, oldSize=[]):
    if gen < burnin:
        return [initSize]
    elif gen % 10 != 0:
        return oldSize
    else:
        incSize = (finalSize-initSize)/(endGen-burnin)
        return [oldSize[0]+10*incSize]

def ne(pop):
    ' calculate effective number of alleles '
    Stat(pop, alleleFreq=[0,1])
    f0 = [0, 0]
    ne = [0, 0]
    for i in range(2):
        freq = pop.dvars().alleleFreq[i][1:]
        f0[i] = 1 - pop.dvars().alleleFreq[i][0]
        if f0[i] == 0:
            ne[i] = 0
        else:
            ne[i] = 1. / sum([(x/f0[i])**2 for x in freq])
    print '%d\t%.3f\t%.3f\t%.3f\t%.3f' % (pop.gen(), f0[0], f0[1], ne[0], ne[1])
    return True

def simulate(incScenario):
    simu = simulator(                                        # create a simulator
        population(subPop=incScenario(0), loci=[1,1],
            infoFields=['fitness']),                         # inital population
        randomMating(subPopSizeFunc=incScenario)           # random mating
    )
    simu.evolve(                            # start evolution
        preOps=[                            # operators that will be applied before evolution
            # initialize locus 0 (for common disease)
            initByFreq(atLoci=[0], alleleFreq=C_f),
            # initialize locus 1 (for rare disease)
            initByFreq(atLoci=[1], alleleFreq=R_f),
        ],
        ops=[                               # operators that will be applied at each gen
            # mutate: k-alleles mutation model
            kamMutator(rate=mu, maxAllele=max_allele),
            # selection on common and rare disease,
            mlSelector([                # multiple loci - multiplicative model
                maSelector(loci=0, fitness=[1,1,1-C_s], wildtype=[0]),
                maSelector(loci=1, fitness=[1,1,1-R_s], wildtype=[0])
            ], mode=SEL_Multiplicative),
            # report generation and popsize and total disease allele frequency.
            pyOperator(func=ne, step=5),
            # monitor time
            ticToc(step=100),
            # pause at any user key input (for presentation purpose)
            pause(stopOnKeyStroke=1)
        ],
        end=endGen
    )


if __name__ == '__main__':
    # get parameters
    (initSize, finalSize, burnin, endGen, growth) = getOptions()
    # 
    from simuPOP import *
    #
    if initSize > finalSize:
        print 'Initial size should be greater than final size'
        sys.exit(1)
    if burnin > endGen:
        print 'Burnin gen should be less than ending gen'
        sys.exit(1)
    if growth == 'linear':
        simulate(lin_exp)
    else:
        simulate(ins_exp)

#end



## 
## #file log/expr.log
## simu = simulator(population(10), randomMating(), rep=2)
## # evaluate an expression in different areas
## print simu.vars(0)
## print simu.population(0).evaluate("gen+1")
## # a statement (no return value)
## simu.population(0).execute("myRep=2+rep*rep")
## simu.population(1).execute("myRep=2*rep")
## print simu.vars(0)
## simu.evolve(
##     ops=[ pyExec("myRep=2+rep*rep") ],
##     gen=1)
## print simu.vars(0)
## #end
## 
## 
## # Note that I can not use pop now, since it is
## # obtained from simu.population(0) which is invalid now.
#### 
## #file log/operatoroutputexpr.log
## outfile="'>>a'+str(rep)+'.txt'"
## simu.evolve(
##   ops = [
##     stat(alleleFreq=[0]),
##     pyEval('alleleFreq[0][0]', outputExpr=outfile)
##   ],
##   gen=1
## )
## print open("a0.txt").read()
## print open("a1.txt").read()
## #end
## 
## #file log/varPlotter.log
## from simuUtil import *
## from simuRPy import *
## 
## simu = simulator(
##     population(size=[50, 50, 100], ploidy=2, loci=[3, 4], infoFields=['migrate_to']),
##     randomMating(), rep=4)
## 
## # migrate
## migr = migrator([[0, .2, .1], [.25, 0, .1], [.1, .2, 0]],
##     mode=ByProbability)
## # and count the size of subpopulations
## stat = stat(popSize=1, stage=PreMating)
## # plot subPopSize. 
## simu.evolve(
##     ops = [
##         migr, 
##         stat,
##         varPlotter('subPopSize', numRep=4, byRep=1, 
##             varDim=3, win=10, title='subPop size', saveAs='log/simuDemo')
##     ],
##     gen=30
## )
## 
## #end
## #PS /usr/bin/convert log/simuDemo16.eps log/simuDemo16.png
## #PS /bin/rm -f log/simuDemo*.eps
## 
## 
## #file log/expcomplex.log
## #
## 
## numSubPop = 100         # number of archipelagos
## numFamilies = 10        # real simulation uses 1000
## numOffspring = 4     # kind of family size
## numReplicate = 1
## loci = [20]*20            # 400 loci on 20 chromosomes
## endGen = 10                 # should be at leat 1000
## maxAllele = 30
## mutationRate = 0.001
## recombinationRate = 0.02
## 
## popSize = numFamilies*numOffspring*numSubPop
## subPopSize = [numFamilies*numOffspring]*numSubPop
## 
## # intializer
## init = initByFreq( alleleFreq=[1./maxAllele]*maxAllele )
## 
## # migration: island model
## #     by proportion, .1 to all others
## #
## migrRate = .1
## # rate[i->i] will be ignored so we can do the following
## migrRates = [[migrRate/(numSubPop-1)]*numSubPop]*numSubPop 
## migrMode = ByProbability
## #
## migrate = migrator(migrRates, mode=migrMode)
## 
## # mutation
## mutate = kamMutator(rate=mutationRate, maxAllele=maxAllele)
## 
## # recombination
## recombine = recombinator( rate = recombinationRate )
## 
## # create a simulator 
## simu = simulator(
##     population(size=subPopSize, ploidy=2, loci=loci),
##     randomMating(numOffspring = numOffspring,
##                          subPopSize=subPopSize) )
## #
## # evolve
## simu.evolve([
##     migrate, 
##     recombine, 
##     mutate,
##     pyEval(r"gen", rep=0),    # report progress
##     endl(rep=-1)
##     ],
##     preOps=[init],
##     gen=endGen)
## 
## 
## #end


#### 
