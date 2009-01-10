#
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
InitByValue(pop, [[0]*5, [1]*5], subPop=([0, 0], [0, 1]), initSex=False)
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

#end

#file log/usePyOperator.log
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
        initByValue([0]*10, subPop=[(0, 0)], initSex=False),
        initByValue([1]*10+[2]*10, subPop=[(0, 1)], initSex=False)],
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

#file log/cppParentChooser.log
# The class myParentsChooser is defined in module myParentsChooser
from myParentsChooser import myParentsChooser

def parentsChooser(pop, sp):
    'How to call a C++ level parents chooser.'
    # create an object with needed information (such as x, y) ...
    pc = myParentsChooser(
        [x.info('x') for x in pop.individuals() if x.sex() == Male],
        [x.info('y') for x in pop.individuals() if x.sex() == Male],
        [x.info('x') for x in pop.individuals() if x.sex() == Female],
        [x.info('y') for x in pop.individuals() if x.sex() == Female])
    while True:
        # return indexes of parents repeatedly
        yield pc.chooseParents()


pop = population(100, loci=[1], infoFields=['x', 'y'])
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
        terminateIf('alleleNum[5][0] == 0 or alleleNum[5][0] == 50')
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
simu.addInfoFields(['father_idx', 'mother_idx'])
simu.setAncestralDepth(1)
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
[ind.intInfo('father_idx')  for ind in sample.individuals()]
#end

################################################################################
#

#file log/simulatorCloneSaveLoad.log
simu = simulator(population(100, loci=[5, 10], infoFields=['x']),
    randomMating(), rep=5)
simu.evolve(preOps=[initByFreq([0.4, 0.6])],
    ops=[], gen=10)
# clone
cloned = simu.clone()
# save and load
simu.save("sample.sim")
loaded = LoadSimulator("sample.sim", randomMating())
# 
simu.numRep()
loaded.numRep()
for rep in range(5):
    assert cloned.population(rep) == loaded.population(rep)

# continue to evolve
simu.evolve(ops=[], gen=10)
simu.gen()
#end

os.remove('sample.sim')

#file log/splitAndMerge.log
from simuPOP import *
pop = population(1000, loci=[1], infoFields=['migrate_to'])
simu = simulator(pop, randomSelection())
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.2, 0.8], at = 3),
        splitSubPop(1, proportions=[0.4, 0.6], at = 5),
        mergeSubPops([0,2], at = 7),
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
        migrator(rate = [0.2], fromSubPop=[0], toSubPop=[1], 
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
        migrator(rate = [0.2], fromSubPop=[0], toSubPop=[1], 
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
        migrator(rate = [0.2], fromSubPop=[0], toSubPop=[1], 
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


#file log/saveQTDT.log
def SaveQTDT(pop, output='', loci=[], 
        fields=[], combine=None, shift=1, **kwargs):
    """ save population in Merlin/QTDT format. The population must have
        pedindex, father_idx and mother_idx information fields.
         
        pop: population to be saved. If pop is a filename, it will be loaded.

        output: base filename. 

        loci: loci to output

        fields: information fields to output

        combine: an optional function to combine two alleles of a diploid 
            individual.

        shift: if combine is not given, output two alleles directly, adding
            this value (default to 1).
    """
    if type(pop) == type(''):
        pop = LoadPopulation(pop)
    if output != '':
        file = output
    else:
        raise exceptions.ValueError, "Please specify output"
    # open data file and pedigree file to write.
    try:
        datOut = open(file + ".dat", "w")
        mapOut = open(file + ".map", "w")
        pedOut = open(file + ".ped", "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + " to write."
    if loci == []:
        loci = range(0, pop.totNumLoci())
    # write dat file
    # 
    if 'affection' in fields:
        outputAffectation = True
        fields.remove('affection')
        print >> datOut, 'A\taffection'
    else:
        outputAffectation = False
    for f in fields:
        print >> datOut, 'T\t%s' % f
    for marker in loci:
        print >> datOut, 'M\t%s' % pop.locusName(marker)
    datOut.close()
    # write map file
    print >> mapOut, 'CHROMOSOME MARKER POSITION'
    for marker in loci:
        print >> mapOut, '%d\t%s\t%f' % (pop.chromLocusPair(marker)[0] + 1, 
            pop.locusName(marker), pop.locusPos(marker))
    mapOut.close()
    # write ped file
    def sexCode(ind):
        if ind.sex() == Male:
            return 1
        else:
            return 2
    # disease status: in linkage affected is 2, unaffected is 1
    def affectedCode(ind):
        if ind.affected():
            return 'a'
        else:
            return 'u'
    #
    pldy = pop.ploidy()
    def writeInd(ind, famID, id, fa, mo):
        print >> pedOut, '%d %d %d %d %d' % (famID, id, fa, mo, sexCode(ind)),
        if outputAffectation:
            print >> pedOut, affectedCode(ind),
        for f in fields:
            print >> pedOut, '%.3f' % ind.info(f),
        for marker in loci:
            for p in range(pldy):
                print >> pedOut, "%d" % (ind.allele(marker, p) + shift), 
        print >> pedOut
    # number of pedigrees
    # get unique pedgree id numbers
    from sets import Set
    peds = Set(pop.indInfo('pedindex', False))
    # do not count peds -1
    peds.discard(-1)
    #
    newPedIdx = 1
    #
    for ped in peds:
        id = 1
        # -1 means no parents
        pastmap = {-1:0}
        # go from generation 2, 1, 0 (for example)
        for anc in range(pop.ancestralDepth(), -1, -1):
            newmap = {-1:0}
            pop.useAncestralPop(anc)
            # find all individual in this pedigree
            for i in range(pop.popSize()):
                ind = pop.individual(i)
                if ind.info('pedindex') == ped:
                    dad = int(ind.info('father_idx'))
                    mom = int(ind.info('mother_idx'))
                    if dad == mom and dad != -1:
                        print ("Something wrong with pedigree %d, father and mother " + \
                            "idx are the same: %s") % (ped, dad)
                    writeInd(ind, newPedIdx, id, pastmap.setdefault(dad, 0), \
                        pastmap.setdefault(mom, 0))
                    newmap[i] = id
                    id += 1
            pastmap = newmap
        newPedIdx += 1
    pedOut.close()

#end

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
                maSelector(locus=0, fitness=[1,1,1-C_s], wildtype=[0]),
                maSelector(locus=1, fitness=[1,1,1-R_s], wildtype=[0])
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
                maSelector(locus=0, fitness=[1,1,1-C_s], wildtype=[0]),
                maSelector(locus=1, fitness=[1,1,1-R_s], wildtype=[0])
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



## #file log/popAndOperator.log
## pop = population(size=[2, 3], ploidy=2, loci=[5, 10],
##     lociPos=[range(0, 5), range(0, 20, 2)],
##     alleleNames=['A', 'C', 'T', 'G'])
## simu = simulator(pop, randomMating(), rep=3)
## simu.evolve(
##     preOps = [ initByFreq([.8, .2])],
##     ops = [
##         stat(alleleFreq=[0, 1], Fst=[1], step=10),
##         kamMutator(rate=0.001, rep=1),
##         kamMutator(rate=0.0001, rep=2)
##     ],
##     gen=10
## )
## #end
## 
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
## 
## #file log/initByFreq.log
## simu = simulator( population(size=[2, 3], loci=[5, 7]),
##     randomMating(), rep=1)
## simu.evolve([
##     initByFreq(alleleFreq=[ [.2, .8], [.8, .2]]),
##     dumper(structure=False)
##   ], 
##   gen=1
##   )
## #end
## 
## #file log/initByValue.log
## simu.evolve([
##     initByValue([1]*5 + [2]*7 + [3]*5 +[4]*7),
##     dumper(structure=False)],
##     gen=1)
## #end
## 
## #file log/pyInit.log
## def initAllele(ind, p, sp):
##   return sp + ind + p
## 
## simu.evolve([
##     pyInit(func=initAllele),
##     dumper(structure=False, dispWidth=2)],
##     gen=1)
## #end
## 
## 
## #file log/kamMutator.log
## simu = simulator(population(size=5, loci=[3, 5]), randomMating())
## simu.evolve(
##     preOps = [initSex()],
##     ops = [
##         kamMutator( rate=[.2, .6, .5], atLoci=[0, 2, 6], maxAllele=9),
##         dumper(structure=False)],
##     gen=1
## )
## #end
## 
## #file log/smmMutator.log
## simu = simulator(population(size=3, loci=[3, 5]), randomMating())
## simu.evolve(
##     preOps = [initByFreq( [.2, .3, .5])],
##     ops = [
##         smmMutator(rate=1,  incProb=.8),
##         dumper(structure=False, stage=PrePostMating)
##     ],
##     gen=1
## )
## #end
## 
## #file log/gsmMutator.log
## simu = simulator(population(size=3, loci=[3, 5]), randomMating())
## simu.evolve(
##     preOps = [initByFreq( [.2, .3, .5])],
##     ops = [
##         gsmMutator(rate=1, p=.8, incProb=.8),
##         dumper(structure=False, stage=PrePostMating)
##     ],
##     gen=1
## )
## 
## import random
## def rndInt():
##   return random.randrange(3, 6)
## 
## simu.evolve(
##     preOps = [initByFreq( [.2, .3, .5])],
##     ops = [
##         gsmMutator(rate=1, func=rndInt, incProb=.8),
##         dumper(structure=False, stage=PrePostMating)
##     ],
##     gen=1
## )
## 
## #end
## 
## 
## #file log/recombinator.log
## simu = simulator(population(4, loci=[4, 5, 6], 
##     infoFields=['father_idx', 'mother_idx']),
##     randomMating())
## simu.evolve(
##     preOps = [initByFreq([.2, .2, .4, .2]), dumper(structure=False) ],
##     ops = [parentsTagger()],
##     postOps = [ dumper(structure=False)],
##     gen=1
## )
## simu.evolve(
##     ops = [
##         parentsTagger(),
##         recombinator(rate=[1, 1, 1], loci=[2, 6, 10])
##     ],
##     postOps = [dumper(structure=False)],
##     gen=1
## )
## #end
## 
## 
## #file log/basicSelector.log
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
##     randomMating())
## s1 = .1
## s2 = .2
## simu.evolve(
##     preOps=[initByFreq(alleleFreq=[.2, .8])],
##     ops = [
##         stat( alleleFreq=[0], genoFreq=[0]),
##         mapSelector(locus=0, fitness={'0-0':(1-s1), '0-1':1, '1-1':(1-s2)}),
##         pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
##     ],
##     gen=300
## )
## #end
## 
## #file log/pySelector.log
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[3], infoFields=['fitness'] ),
##     randomMating()
## )
## 
## s1 = .2
## s2 = .3
## # the second parameter gen can be used for varying selection pressure
## def sel(arr, gen=0):
##   if arr[0] == 1 and arr[1] == 1:
##     return 1 - s1
##   elif arr[0] == 1 and arr[1] == 2:
##     return 1
##   elif arr[0] == 2 and arr[1] == 1:
##     return 1
##   else:
##     return 1 - s2
## 
## # test func
## print sel([1, 1])
## 
## simu.evolve(
##     preOps=[initByFreq(alleleFreq=[.2, .8])],
##     ops = [
##         stat( alleleFreq=[0], genoFreq=[0]),
##         pySelector(loci=[0, 1], func=sel),
##         pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=25)
##     ],
##     gen=100
## )
## #end
## 
## #file log/pySubset.log
## simu = simulator(population(size=[2, 3], loci=[3, 4], infoFields=['fitness']),
##     randomMating())
## simu.evolve(
##     preOps = [initByFreq([.3, .5, .2])],
##     ops = [
##         pySubset( [1, -1, -1, 1, -1] ),
##         dumper(structure=False, stage=PrePostMating)
##     ],
##     gen=1
## )
## #end
## 
## 
## 
## #turnOnDebug(DBG_ALL)
## #turnOnDebug(DBG_SIMULATOR)
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
##     mode=MigrByProbability)
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
## #file log/rng.log
## print AvailableRNGs()
## print rng().name()
## SetRNG("taus2", seed=10)
## print rng().name()
## #end
## 
## #file log/rngrand.log
## r=rng()
## #help(RNG)
## for n in range(1, 10):
##   print r.randBinomial(10, .7),
## #end
## 
## 
## 
## 
## #file log/extgenostru.log
## pop = population(1, loci=[2, 3, 4])
## print pop.numLoci(1)
## print pop.locusPos(2)
## dis = pop.lociPos()
## print dis
## print pop.locusPos(2)
## print pop.lociPos()
## #end
## 
## #file log/extgenotype.log
## pop = population(1, loci=[2, 3, 4])
## InitByFreq(pop, [.2, .8])
## Dump(pop, alleleOnly=1)
## ind = pop.individual(0)
## print ind.allele(1, 1)
## ind.setAllele(3, 1, 1)
## Dump(pop, alleleOnly=1)
## a = ind.genotype()
## print a
## a = ind.genotype(1)
## print a
## a = ind.genotype(1, 2)
## print a
## a[2]=4
## # the allele on the third chromosome has been changed
## Dump(pop, alleleOnly=1)
## #end
## 
## #file log/extother.log
## pop = population(1, loci=[2, 3, 4])
## ind = pop.individual(1)
## print ind.sex()
## print ind.sexChar()
## ind.setSex(Female)
## ind.setAffected(True)
## print ind.tag()
## ind.setTag([1, 2])
## Dump(pop)
## #end
## 
## #file log/extsimu.log
## pop = population(1, loci=[2, 3, 4])
## simu = simulator(pop, randomMating(), rep=3)
## pop1 = simu.population(1)
## ind1 = pop1.individual(0)
## ind1.setAllele(3, 0)
## Dump(pop1)
## #end
## 
## #file log/extoperator.log
## pop = population(1, loci=[2, 3, 4])
## simu = simulator(pop, randomMating(), rep=3)
## simu.evolve(
##     ops = [ pyEval(stmts="pop=simu.population(rep)")],
##     gen=1)
## #end
#### 
## #file log/expLD.log
## #
## # this is an example of observing decay of LD
## from simuUtil import *
## from simuRPy import *
## 
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[2]),
##     randomMating(),
##     rep=4)
## 
## # see the change of allele/genotype/heplotype numbers as
## # the result of genetic drift.
## init = initByValue([1, 2, 2, 1])
## count = stat(LD=[0, 1])
## recombine = recombinator( rate=0.1 )
## simu.evolve(
##     ops = [
##         recombine,
##         count,
##         pyEval(r'"%.4f\t" % LD[0][1]'),
##    #varPlotter(expr='LD[0][1]', title='Linkage disequilibrium',
##    #  numRep = 4, ytitle='LD', saveAs='LD')
##    ],
##    preOps=[init],
##    gen=10
## )
## 
## #end
## 
## 
## 
## #file log/expcomplex.log
## #
## 
## numSubPop = 100     # number of archipelagos
## numFamilies = 10    # real simulation uses 1000
## numOffspring = 4   # kind of family size
## numReplicate = 1
## loci = [20]*20      # 400 loci on 20 chromosomes
## endGen = 10         # should be at leat 1000
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
## #   by proportion, .1 to all others
## #
## migrRate = .1
## # rate[i->i] will be ignored so we can do the following
## migrRates = [[migrRate/(numSubPop-1)]*numSubPop]*numSubPop 
## migrMode  = MigrByProbability
## #
## migrate =  migrator(migrRates, mode=migrMode)
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
##                        subPopSize=subPopSize) )
## #
## # evolve
## simu.evolve(
##     preOps=[init],
##     ops = [
##         migrate, 
##         recombine, 
##         mutate,
##         pyEval(r"gen", rep=0),  # report progress
##         pyOutput('\n', rep=-1)
##     ],
##     gen=endGen
## )
## 
## 
## #end
## 
## #file log/expmigration.log
## # this is an example of complex population size change.
## # for endl and tab
## from simuUtil import *
## 
## #number of cities
## nc = 6
## 
## # how to change subpop size?
## def changeSPSize(gen, oldSize=[]):
##   size = [0]*len(oldSize)
##   for i in range(0, len(size)):
##     size[i] = oldSize[i]*1.2
##   if size[i] > 1000:
##     size[i] /= 2
##   return size
## 
## # migration between subpopulaitons
## rates = []
## for i in range(nc):
##     rates.append([0.]*nc)
## 
## for i in range(1, nc-1):
##   rates[i][i+1]=0.05
##   rates[i][i-1]=0.05
## 
## rates[0][1] = 0.1
## rates[nc-1][nc-2] = 0.1
## 
## # print rates
## print rates
## migr = migrator(rate=rates, mode=MigrByProbability)
## 
## # initially, we need to set everyone to middle subpop
## initMigr = migrator(rate=[[1]], mode=MigrByProportion,
##        fromSubPop=[0], toSubPop=[nc/2])
## 
## pop = population(size=500)
## 
## # the new popsize relies on a variable newSPSize
## # which is calculated from subPopSize bu newSize operator
## simu = simulator(pop,
##     randomMating(subPopSize=changeSPSize) )
## 
## # evolve!
## simu.evolve(
##     ops = [migr, stat(popSize=True),
##        pyEval('list(subPopSize)'), 
##        pyOutput('\n')
##     ],
##     preOps = [ initMigr ],
##     gen=10
## )
## 
## #end
## 
## # need reich.py
## 
## #file log/absIndex.log
## pop = population(size=[20, 30], loci=[5, 6])
## print pop.chromLocusPair(7)
## print pop.absLocusIndex(1, 1)
## print pop.absIndIndex(10, 1)
## print pop.subPopIndPair(40)
## #end
## 
## 
## #file log/popAndOperator.log
## simu = simulator(pop, randomMating(), rep=3)
## simu.evolve(
##     preOps = [ initByFreq([.8, .2])],
##     ops = [
##         stat(alleleFreq=[0, 1], Fst=[1], step=10),
##         kamMutator(rate=0.001, rep=1),
##         kamMutator(rate=0.0001, rep=2)
##     ],
##     gen=10
## )
## #end
## 
## #file log/genotype.log
## pop = population(size=[3, 2], loci=[2])
## # single allele access
## for ind in pop.individuals(1):
##     for marker in range(pop.totNumLoci()):
##         ind.setAllele(marker % 2, marker, 0)
##         ind.setAllele(marker % 2, marker, 1)
##         print '%d %d ' % (ind.allele(marker, 0), ind.allele(marker, 1))
## 
## # batch access
## ind = pop.individual(4)
## geno = ind.genotype()
## print geno
## geno[2] = 3
## print ind.genotype()
## # direct modification of the underlying genotype
## geno[2:4] = [3, 4]
## print ind.genotype()
## # set genotype
## ind.setGenotype([2, 1])
## print geno
## # print genotypes of all individuals in the second subpopulation.
## print pop.genotype(1)
## #end
## 
## 
## #file log/localNamespace.log
## pop = population(size=[1000, 2000], loci=[1])
## InitByFreq(pop, [0.2, 0.8])
## Stat(pop, popSize=1, alleleFreq=[0])
## print pop.evaluate('alleleNum[0][0] + alleleNum[0][1]')
## pop.execute('newPopSize=int(popSize*1.5)')
## ListVars(pop.vars(), level=1, useWxPython=False)
## # this variable is 'local' to the population and is
## # not available in the main namespace
## newPopSize
## #
## simu = simulator(population(10), randomMating(), rep=2)
## # evaluate an expression in different areas
## print simu.vars(1)
## # a statement (no return value)
## simu.population(0).execute("myRep=2+rep*rep")
## simu.population(1).execute("myRep=2*rep")
## print simu.vars(0)
## #end
## 
## 
## #file log/pyEval.log
## simu = simulator(population(100, loci=[1]),
##     randomMating(), rep=2)
## simu.evolve(
##     preOps = [initByFreq([0.2, 0.8])],
##     ops = [ stat(alleleFreq=[0]),
##         pyExec('myNum = alleleNum[0][0] * 2'),
##         pyEval(r'"gen %d, rep %d, num %d, myNum %d\n"' \
##             ' % (gen, rep, alleleNum[0][0], myNum)')
##         ],
##     gen=3
## )
## #end
## 
## 
## #turnOnDebug(DBG_SIMULATOR)
## #turnOnDebug(DBG_UTILITY)
## 
## 
## # Note that I can not use pop now, since it is
## # obtained from simu.population(0) which is invalid now.
#### 
## 
## #file log/generator_random.log
## from random import randint
## 
## def randomChooser(pop, sp):
##     males = [x for x in range(pop.subPopSize(sp)) \
##         if pop.individual(x, sp).sex() == Male \
##             and pop.individual(x, sp).info('age') > 30]
##     females = [x for x in range(pop.subPopSize(sp)) \
##         if pop.individual(x, sp).sex() == Female \
##             and pop.individual(x, sp).info('age') > 30]
##     nm = len(males)
##     nf = len(females)
##     while True:
##         yield males[randint(0, nm-1)], females[randint(0, nf-1)]
## 
## pop = population(size=[1000, 200], loci=[1], infoFields=['age'])
## # this will initialize sex randomly
## InitByFreq(pop, [0.2, 0.8])
## for ind in pop.individuals():
##     ind.setInfo(randint(0, 60), 'age')
## 
## rc1 = randomChooser(pop, 0)
## for i in range(5):
##     print rc1.next(),
## 
## rc2 = randomChooser(pop, 1)
## for i in range(5):
##     print rc2.next(),
## 
## #end
## 
## 
## #file log/ifElse.log
## from simuRPy import *
## from simuUtil import *
## numRep=4
## popSize=100
## endGen=50
## 
## simu = simulator(population(size=popSize, loci=[1]),
##     randomMating(), rep=numRep)
## simu.evolve(
##     preOps = [ initByValue([1, 1])],
##     ops = [
##         # penetrance, additve penetrance
##         maPenetrance(locus=0, wildtype=[1], penetrance=[0, 0.5, 1]),
##         # count number of affected
##         stat(numOfAffected=True),
##         # introduce disease if no one is affected
##         ifElse(cond='numOfAffected==0',
##             ifOp=kamMutator(rate=0.01, maxAllele=2)),
##         # expose affected status
##         pyExec('pop.exposeAffectedness()', exposePop=True),
##         # plot affected status
##         varPlotter(expr='affected', plotType="image", byRep=1, update=endGen, 
##             varDim=popSize, win=endGen, numRep=numRep,
##             title='affected status', saveAs="ifElse")
##     ],
##     gen=endGen,
##     dryrun=False
## )
## #end
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
## migrMode = MigrByProbability
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
## 
## #file log/expmigration.log
## # this is an example of complex population size change.
## # for endl and tab
## from simuUtil import *
## 
## #number of cities
## nc = 6
## 
## # how to change subpop size?
## def changeSPSize(gen, oldSize):
##     size = [0]*len(oldSize)
##     for i in range(0, len(size)):
##         size[i] = oldSize[i]*1.2
##     if size[i] > 1000:
##         size[i] /= 2
##     return size
## 
## # migration between subpopulaitons
## rates = []
## for i in range(nc):
##     rates.append([0.]*nc)
## #
## for i in range(1, nc-1):
##     rates[i][i+1]=0.05
##     rates[i][i-1]=0.05
## 
## #
## rates[0][1] = 0.1
## rates[nc-1][nc-2] = 0.1
## 
## # print rates
## print rates
## migr = migrator(rate=rates, mode=MigrByProbability)
## 
## # initially, we need to set everyone to middle subpop
## initMigr = migrator(rate=[[1]], mode=MigrByProportion,
##     fromSubPop=[0], toSubPop=[nc/2])
## 
## pop = population(size=500)
## 
## # the new popsize relies on a variable newSPSize
## # which is calculated from subPopSize bu newSize operator
## simu = simulator(pop,
##     randomMating(subPopSize=changeSPSize) )
## 
## # evolve!
## simu.evolve( [
##     migr, stat(popSize=True),
##     pyEval('list(subPopSize)'), endl()],
##     preOps = [ initMigr ], gen=10
## )
## 
## #end
#### 
## #file log/src_genoStruTrait.log
## # create a population, most parameters have default values
## pop = population(size=[2,3], ploidy=2, loci=[5, 10],
##     lociPos=[range(0, 5), range(0, 20, 2)],
##     alleleNames=['A', 'C', 'T', 'G'],
##     maxAllele=3)
## print pop.popSize()
## print pop.ploidy()
## print pop.ploidyName()
## print pop.numChrom()
## print pop.locusPos(2)
## print pop.alleleName(1)
## # get the fourth individual of the population
## ind = pop.individual(3)
## # access genotypic structure info
## print ind.ploidy()
## print ind.numChrom()
## print ind.numLoci(0)
## print ind.genoSize()
## # and from simulator level
## simu = simulator(pop, randomMating(), rep=3)
## print simu.numChrom()
## #end
## 
## #file log/src_population.log
## # use of population function
## # a Wright-Fisher population
## WF = population(size=100, ploidy=1, loci=[1])
## 
## # a diploid population of size 10
## # there are two chromosomes with 5 and 7 loci respectively
## pop = population(size=[2, 8], ploidy=2, loci=[5, 7])
## 
## # a population with SNP markers (with names A, C, T, G)
## # range() are python functions
## pop = population(size=[2,3], ploidy=2, loci=[5, 10],
##     lociPos=[range(0, 5), range(0, 20, 2)],
##     alleleNames=['A', 'C', 'T', 'G'],
##     maxAllele=3)
## 
## #
## # population structure functions
## print pop.popSize()
## print pop.numSubPop()
## print pop.subPopSize(0)
## print pop.subPopSizes()
## print pop.subPopBegin(1)
## print pop.subPopEnd(1)
## print pop.subPopIndPair(3)
## print pop.absIndIndex(1, 1)
## 
## #
## # save and load population
## # save it in various formats, default format is "txt"
## pop = population(1000, loci=[2, 5, 10])
## pop.save("sample.pop")
## 
## # load it in another population
## pop1 = LoadPopulation("sample.pop")
## #end
## 
## 
## 
## #file log/src_operator.log
## simu = simulator(population(1, loci=[3]), randomSelection(), rep=2)
## op1 = pyOutput("a", begin=5, end=20, step=3)
## op2 = pyOutput("a", begin=-5, end=-1, step=2)
## op3 = pyOutput("a", at=[2, 5, 10])
## op4 = pyOutput("a", at=[-10, -5, -1])
## simu.evolve(
##     ops = [ pyEval(r"str(gen)+'\n'", begin=5, end=-1, step=2)],
##     gen=10
## )
## #
## #
## # parameter output 
## simu = simulator(population(100, loci=[3]), randomMating(), rep=2)
## simu.evolve(
##     preOps=[
##         initByFreq([0.2, 0.8], rep=0),
##         initByFreq([0.5, 0.5], rep=1) ],
##     ops = [
##         stat(alleleFreq=[0]),
##         pyEval('alleleFreq[0][0]', output='a.txt')
##     ],
##     gen=1
## )
## # only from rep 1
## print open('a.txt').read()
## 
## simu.evolve(
##     ops = [
##         stat(alleleFreq=[0]),
##         pyEval('alleleFreq[0][0]', output='>>a.txt')
##     ],
##     gen=1)
## # from both rep0 and rep1
## print open("a.txt").read()
## 
## outfile='>>>a.txt'
## simu.evolve(
##     ops = [
##         stat(alleleFreq=[0]),
##         pyEval('alleleFreq[0][0]', output=outfile),
##         pyOutput("\t", output=outfile),
##         pyOutput("\n", output=outfile, rep=0)
##     ],
##     gen=1
## )
## print open("a.txt").read()
## #
## # Output expression
## outfile="'>>a'+str(rep)+'.txt'"
## simu.evolve(
##     ops = [
##         stat(alleleFreq=[0]),
##         pyEval('alleleFreq[0][0]', outputExpr=outfile)
##     ],
##     gen=1
## )
## print open("a0.txt").read()
## print open("a1.txt").read()
## #end
## os.remove('a.txt')
## 
## 
## 
## #file log/src_initByFreq.log
## simu = simulator( 
##     population(size=[2, 3], loci=[5, 7], maxAllele=1),
##     randomMating(), rep=1)
## simu.evolve([
##     initByFreq(alleleFreq=[ [.2, .8], [.8, .2]]),
##     dumper(structure=False)
##   ],
##   gen=1)
## #end
## 
## #file log/src_initByValue.log
## simu = simulator(
##     population(size=[2, 3], loci=[5, 7], maxAllele=9),
##     randomMating(), rep=1)
## simu.evolve([
##     initByValue([1]*5 + [2]*7 + [3]*5 +[4]*7),
##     dumper(structure=False)],
##     gen=1)
## #end
## 
## #file log/src_pyInit.log
## def initAllele(ind, p, sp):
##   return sp + ind + p
## 
## simu = simulator( 
##     population(size=[2, 3], loci=[5, 7]),
##     randomMating(), rep=1)
## simu.evolve([
##     pyInit(func=initAllele),
##     dumper(structure=False, dispWidth=2)],
##     gen=1)
## #end
## 
## 
## #file log/src_mating.log
## # arbitrary demographic model
## def lin_inc(gen, oldsize=[]):
##     return [10+gen]*5
## 
## simu = simulator(
##     population(size=[50]*5, loci=[1]),
##     randomMating(subPopSize=lin_inc)
## )
## simu.evolve(
##     preOps = [initSex()],
##     ops = [
##         stat(popSize=True),
##         pyEval(r'"%d %d\n"%(gen, subPop[0]["popSize"])'),
##     ],
##     gen = 5
## )
## 
## #
## # control the number of offspring per mating event
## # famSizes is only defined when DBG_MATING is defined
## TurnOnDebug(DBG_MATING)
## simu = simulator(population(50, loci=[1]),
##     randomMating(numOffspring=(UniformDistribution, 2, 5)))
## simu.evolve(ops=[], gen=1)
## print simu.population(0).dvars().famSizes
## TurnOffDebug(DBG_MATING)
## #end
## 
## 
## #file log/src_mapSelector.log
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
##     randomMating())
## s1 = .1
## s2 = .2
## simu.evolve([
##     stat( alleleFreq=[0], genoFreq=[0]),
##     mapSelector(locus=0, fitness={'0-0':(1-s1), '0-1':1, '1-1':(1-s2)}),
##     pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2, .8])],
##     gen = 300)
## #end
## 
## 
## #file log/src_maSelector.log
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
##     randomMating())
## s1 = .1
## s2 = .2
## simu.evolve(
##     preOps=[initByFreq(alleleFreq=[.2, .8])],
##     ops = [
##         stat( alleleFreq=[0], genoFreq=[0]),
##         maSelector(locus=0, fitness=[1-s1, 1, 1-s2]),
##         pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
##     ],
##     gen = 300)
## #end
## 
## #file log/src_pySelector.log
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[3], infoFields=['fitness'] ),
##     randomMating())
## 
## s1 = .2
## s2 = .3
## # the second parameter gen can be used for varying selection pressure
## def sel(arr, gen=0):
##   if arr[0] == 1 and arr[1] == 1:
##     return 1 - s1
##   elif arr[0] == 1 and arr[1] == 2:
##     return 1
##   elif arr[0] == 2 and arr[1] == 1:
##     return 1
##   else:
##     return 1 - s2
## 
## # test func
## print sel([1, 1])
## 
## simu.evolve([
##     stat( alleleFreq=[0], genoFreq=[0]),
##     pySelector(loci=[0, 1], func=sel),
##     pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=25)
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2, .8])],
##     gen = 100)
## #end
## 
## 
## #file log/src_mlPenetrance.log
## pop = population(1000, loci=[3])
## InitByFreq(pop, [0.3, 0.7])
## pen = []
## for loc in (0, 1, 2):
##     pen.append( maPenetrance(locus=loc, wildtype=[1],
##         penetrance=[0, 0.3, 0.6] ) )
## 
## # the multi-loci penetrance
## MlPenetrance(pop, mode=PEN_Multiplicative, peneOps=pen)
## Stat(pop, numOfAffected=True)
## print pop.dvars().numOfAffected
## #end
## 
## 
## #file log/src_pyPenetrance.log
## pop = population(1000, loci=[3])
## InitByFreq(pop, [0.3, 0.7])
## def peneFunc(geno):
##     p = 1
##     for l in range(len(geno)/2):
##         p *= (geno[l*2]+geno[l*2+1])*0.3
##     return p
## 
## PyPenetrance(pop, func=peneFunc, loci=(0, 1, 2))
## Stat(pop, numOfAffected=True)
## print pop.dvars().numOfAffected
## #
## # You can also define a function, that returns a penetrance
## # function using given parameters
## def peneFunc(table):
##     def func(geno):
##       return table[geno[0]][geno[1]]
##     return func
## 
## # then, given a table, you can do
## PyPenetrance(pop, loci=(0, 1, 2),
##     func=peneFunc( ((0, 0.5), (0.3, 0.8)) ) )
## #end
## 
## 
## 
## #turnOnDebug(DBG_ALL)
## #turnOnDebug(DBG_SIMULATOR)
## 
## #file log/src_ifElse.log
## simu = simulator(
##     population(size=1000, loci=[1]),
##     randomMating(), rep=4)
## simu.evolve(
##   preOps = [ initByValue([1, 1])],
##   ops = [
##     # penetrance, additve penetrance
##     maPenetrance(locus=0, wildtype=[1], penetrance=[0, 0.5, 1]),
##     # count number of affected
##     stat(numOfAffected=True),
##     # introduce disease if no one is affected
##     ifElse(cond='numOfAffected==0',
##       ifOp=kamMutator(rate=0.01, maxAllele=2)),
##     ifElse(cond='numOfAffected==0',
##         ifOp=pyEval(r'"No affected at gen %d\n" % gen'))
##   ],
##   gen = 50
## )
## #end
## 
## #file log/src_noneOp.log
## # this may be set from command line option
## savePop = False
## # then, saveOp is defined accordingly
## if savePop:
##     saveOp = savePopulation(output='a.txt')
## else:
##     saveOp = noneOp()
## 
## simu = simulator(population(10), randomMating())
## simu.evolve(preOps=[initSex()],
##     ops=[saveOp], gen=1)
## #end
## 
## 
## #file log/src_rng.log
## print AvailableRNGs()
## print rng().name()
## SetRNG("taus2", seed=10)
## print rng().name()
## #end
## 
## #file log/src_rngrand.log
## r=rng()
## #help(RNG)
## for n in range(1, 10):
##     print r.randBinomial(10, .7),
## #end
## 
## 
## #file log/src_splitSubPop.log
## pop = population(size=10, loci=[2, 6], infoFields=['migrate_to'])
## InitByFreq(pop, [.2, .4, .4])
## SplitSubPop(pop, which=0, sizes=[2, 8], randomize=False)
## print pop.subPopSizes()
## #end
## 
## 
## #file log/src_mlSelector.log
## simu = simulator(
##     population(size=10, ploidy=2, loci=[2], 
##     infoFields=['fitness', 'spare']),
##     randomMating())
## simu.evolve(
##     [ mlSelector([
##          mapSelector(locus=0, fitness={'0-0':1, '0-1':1, '1-1':.8}),
##          mapSelector(locus=1, fitness={'0-0':1, '0-1':1, '1-1':.8}),
##          ], mode=SEL_Additive),
##     ],
##     preOps = [
##         initByFreq(alleleFreq=[.2, .8])
##     ],
##     gen = 2
## )
## #end
## 
## 
## #file log/src_mapPenetrance.log
## pop = population(size=[2,8], ploidy=2, loci=[2] )
## InitByFreq(pop, [.2, .8])
## MapPenetrance(pop, locus=0, 
##     penetrance={'0-0':0, '0-1':1, '1-1':1})
## Stat(pop, numOfAffected=1)
## #end
## 
## 
